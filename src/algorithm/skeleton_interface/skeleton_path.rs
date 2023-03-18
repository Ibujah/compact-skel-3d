use anyhow::Result;
use nalgebra::{MatrixXx1, MatrixXx3, Vector3};

use super::skeleton_interface::{IterPartialEdge, SkeletonInterface3D};
use crate::geometry::geometry_operations;

#[derive(Copy, Clone)]
pub enum PathPart {
    PartialNode(usize),
    PartialEdge(usize),
}

pub enum State {
    Computing,
    Closed,
}

pub struct SkeletonPath {
    components: Vec<PathPart>,
    opt_ind_pedge_last: Option<usize>,
}

impl SkeletonPath {
    pub fn new(ind_pedge: usize) -> SkeletonPath {
        SkeletonPath {
            components: Vec::new(),
            opt_ind_pedge_last: Some(ind_pedge),
        }
    }

    pub fn mesh_path(&self, skeleton_interface: &SkeletonInterface3D) -> Vec<usize> {
        let mut path = Vec::new();
        for ind1 in 0..self.components.len() {
            let ind2 = (ind1 + 1) % self.components.len();
            match (self.components[ind1], self.components[ind2]) {
                (PathPart::PartialNode(ind_pnode), PathPart::PartialNode(_)) => {
                    path.push(
                        skeleton_interface
                            .get_partial_node_uncheck(ind_pnode)
                            .corner(),
                    );
                }
                (_, _) => (),
            }
        }
        path
    }

    pub fn append_last(&mut self, skeleton_interface: &mut SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_last) = self.opt_ind_pedge_last {
            let edge = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_last)
                .edge();
            if !edge.is_computed() {
                let ind_edge = edge.ind();
                skeleton_interface.propagate_edge(ind_edge)?;
            }
            self.components.push(PathPart::PartialEdge(ind_pedge_last));

            let ind_pedge_next = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_last)
                .partial_edge_next()
                .ok_or(anyhow::Error::msg("No next partial edge"))?
                .ind();

            self.opt_ind_pedge_last = Some(ind_pedge_next);

            self.check_loop()
        } else {
            Ok(State::Closed)
        }
    }

    fn check_loop(&mut self) -> Result<State> {
        if let Some(ind_pedge_last) = self.opt_ind_pedge_last {
            let part_first = self.components.first().unwrap();
            let looped = match part_first {
                &PathPart::PartialNode(_) => todo!(),
                &PathPart::PartialEdge(ind_pedge) => ind_pedge == ind_pedge_last,
            };

            if looped {
                self.opt_ind_pedge_last = None
            }
        }

        if self.opt_ind_pedge_last.is_none() {
            Ok(State::Closed)
        } else {
            Ok(State::Computing)
        }
    }

    pub fn rotate_last(&mut self, skeleton_interface: &SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_last) = self.opt_ind_pedge_last {
            let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge_last);
            let pedge_neigh = pedge.partial_edge_neighbor();

            let alve = pedge_neigh.partial_alveola().alveola();
            let pedge_next = if alve.is_full() {
                pedge_neigh
                    .partial_edge_next()
                    .ok_or(anyhow::Error::msg("No next partial edge"))?
            } else {
                pedge_neigh.partial_edge_opposite()
            };
            let ind_pnode = pedge
                .partial_node_first()
                .ok_or(anyhow::Error::msg("No first node"))?
                .ind();
            if let Some(&plast) = self.components.last() {
                if let PathPart::PartialNode(nod) = plast {
                    if nod != ind_pnode {
                        self.components.push(PathPart::PartialNode(ind_pnode));
                    }
                } else {
                    self.components.push(PathPart::PartialNode(ind_pnode));
                }
            } else {
                self.components.push(PathPart::PartialNode(ind_pnode));
            }
            let ind_pedge_new = pedge_next.ind();
            self.opt_ind_pedge_last = Some(ind_pedge_new);
            self.check_loop()
        } else {
            Ok(State::Closed)
        }
    }

    pub fn ind_last_partial_edge(&self) -> Option<usize> {
        self.opt_ind_pedge_last
    }

    pub fn last_partial_edge<'a, 'b>(
        &self,
        skeleton_interface: &'b SkeletonInterface3D<'a>,
    ) -> Option<IterPartialEdge<'a, 'b>> {
        if let Some(ind_pedge_last) = self.opt_ind_pedge_last {
            Some(skeleton_interface.get_partial_edge_uncheck(ind_pedge_last))
        } else {
            None
        }
    }

    pub fn ind_partial_edges(&self) -> Vec<usize> {
        self.components
            .iter()
            .filter_map(|pp| match pp {
                &PathPart::PartialEdge(ind_pedge) => Some(ind_pedge),
                &PathPart::PartialNode(_) => None,
            })
            .collect()
    }

    pub fn ind_alveolae(&self, skeleton_interface: &SkeletonInterface3D) -> Vec<usize> {
        self.components
            .iter()
            .filter_map(|pp| match pp {
                &PathPart::PartialEdge(ind_pedge) => Some(
                    skeleton_interface
                        .get_partial_edge_uncheck(ind_pedge)
                        .partial_alveola()
                        .alveola()
                        .ind(),
                ),
                &PathPart::PartialNode(_) => None,
            })
            .collect()
    }

    pub fn closable_path(&self, skeleton_interface: &SkeletonInterface3D) -> Result<bool> {
        let mut has_deg1 = false;
        for ind in 0..self.components.len() {
            let ind_next = (ind + 1) % self.components.len();
            match (self.components[ind], self.components[ind_next]) {
                (PathPart::PartialEdge(ind_pedge), PathPart::PartialNode(_)) => {
                    let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge);
                    if pedge.partial_alveola().alveola().label().is_some() {
                        return Ok(false);
                    }
                    let deg = pedge
                        .partial_edge_next()
                        .ok_or(anyhow::Error::msg("Next edge does not exist"))?
                        .edge()
                        .degree();
                    if deg != 1 && deg != 2 {
                        return Ok(false);
                    }
                    if deg == 1 {
                        has_deg1 = true;
                    }
                }
                (_, _) => (),
            };
        }
        return Ok(has_deg1);
    }

    pub fn nodes(&self, skeleton_interface: &SkeletonInterface3D) -> Vec<usize> {
        let mut nodes: Vec<usize> = self
            .components
            .iter()
            .filter_map(|&cmp| match cmp {
                PathPart::PartialNode(ind_pnode) => Some(
                    skeleton_interface
                        .get_partial_node_uncheck(ind_pnode)
                        .node()
                        .ind(),
                ),
                PathPart::PartialEdge(_) => None,
            })
            .collect();
        nodes.sort();
        nodes
    }

    pub fn basis_spheres(
        &self,
        skeleton_interface: &SkeletonInterface3D,
    ) -> Result<(MatrixXx3<f32>, MatrixXx1<f32>)> {
        let ind_nodes = self.nodes(skeleton_interface);
        let mut center_mat = MatrixXx3::<f32>::zeros(ind_nodes.len());
        let mut radius_mat = MatrixXx1::<f32>::zeros(ind_nodes.len());

        for i in 0..ind_nodes.len() {
            let ind_node = ind_nodes[i];
            let tet_vert: Vec<Vector3<f32>> = skeleton_interface
                .get_node_uncheck(ind_node)
                .delaunay_tetrahedron()
                .iter()
                .map(|&ind| {
                    skeleton_interface
                        .get_mesh()
                        .get_vertex(ind)
                        .unwrap()
                        .vertex()
                })
                .collect();
            let (center, radius) = geometry_operations::center_and_radius(
                [tet_vert[0], tet_vert[1], tet_vert[2], tet_vert[3]],
                None,
            )
            .ok_or(anyhow::Error::msg("Could not find sphere center"))?;
            center_mat[(i, 0)] = center[0];
            center_mat[(i, 1)] = center[1];
            center_mat[(i, 2)] = center[2];
            radius_mat[i] = radius;
        }

        Ok((center_mat, radius_mat))
    }

    pub fn follow_singular_path(
        &mut self,
        skeleton_interface: &mut SkeletonInterface3D,
    ) -> Result<()> {
        loop {
            if let Some(pedge) = self.last_partial_edge(skeleton_interface) {
                if pedge.is_singular() {
                    self.append_last(skeleton_interface)?;
                } else {
                    self.rotate_last(skeleton_interface)?;
                }
            } else {
                break;
            }
        }
        Ok(())
    }
}
