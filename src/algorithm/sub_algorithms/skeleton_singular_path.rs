use anyhow::Result;
use nalgebra::*;

use super::SkeletonInterface3D;

#[derive(Copy, Clone)]
pub enum PathPart {
    PartialNode(usize),
    PartialEdge(usize),
}

pub enum State {
    Computing,
    Closed,
}

pub struct SkeletonSingularPath {
    components: Vec<PathPart>,
    opt_ind_pedge_last: Option<usize>,
}

impl SkeletonSingularPath {
    pub fn create(ind_pedge: usize) -> SkeletonSingularPath {
        SkeletonSingularPath {
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

    pub fn components(&self) -> &Vec<PathPart> {
        &self.components
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

            self.check_loop(&skeleton_interface)
        } else {
            Ok(State::Closed)
        }
    }

    fn check_loop(&mut self, skeleton_interface: &SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_last) = self.opt_ind_pedge_last {
            let part_first = self.components.first().unwrap();
            let looped = match part_first {
                &PathPart::PartialNode(ind_pnode) => {
                    ind_pnode
                        == skeleton_interface
                            .get_partial_edge_uncheck(ind_pedge_last)
                            .partial_node_first()
                            .unwrap()
                            .ind()
                }
                &PathPart::PartialEdge(ind_pedge) => ind_pedge == ind_pedge_last,
            };

            if looped {
                let part_last = self.components.first().unwrap();
                if let (&PathPart::PartialNode(ind_pnode1), &PathPart::PartialNode(ind_pnode2)) =
                    (part_first, part_last)
                {
                    if ind_pnode1 == ind_pnode2 {
                        self.components.pop();
                    }
                }
                self.opt_ind_pedge_last = None;
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
            self.check_loop(skeleton_interface)
        } else {
            Ok(State::Closed)
        }
    }

    pub fn ind_last_partial_edge(&self) -> Option<usize> {
        self.opt_ind_pedge_last
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

    pub fn basis_spheres_matrices(
        &self,
        skeleton_interface: &SkeletonInterface3D,
    ) -> Result<(MatrixXx3<f32>, MatrixXx1<f32>)> {
        let ind_nodes = self.nodes(skeleton_interface);
        let mut center_mat = MatrixXx3::<f32>::zeros(ind_nodes.len());
        let mut radius_mat = MatrixXx1::<f32>::zeros(ind_nodes.len());

        for i in 0..ind_nodes.len() {
            let ind_node = ind_nodes[i];
            let (center, radius) = skeleton_interface
                .get_node_uncheck(ind_node)
                .center_and_radius()?;

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
            if let Some(ind_pedge_last) = self.ind_last_partial_edge() {
                let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge_last);
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

    pub fn halfedges_path(&self, skeleton_interface: &SkeletonInterface3D) -> Result<Vec<usize>> {
        let vec_bnd_path_vert = self.mesh_path(skeleton_interface);
        let mut mesh_path_hedge = Vec::new();
        for ind1 in 0..vec_bnd_path_vert.len() {
            let ind2 = (ind1 + 1) % vec_bnd_path_vert.len();
            let ind_vertex1 = vec_bnd_path_vert[ind1];
            let ind_vertex2 = vec_bnd_path_vert[ind2];

            let hedge = skeleton_interface
                .get_mesh()
                .is_edge_in(ind_vertex1, ind_vertex2)
                .ok_or(anyhow::Error::msg(format!(
                    "Halfedge ({}, {}) is not on the boundary",
                    ind_vertex1, ind_vertex2
                )))?;
            mesh_path_hedge.push(hedge.ind());
        }
        Ok(mesh_path_hedge)
    }

    pub fn alveolae_path(&self, skeleton_interface: &SkeletonInterface3D) -> Result<Vec<usize>> {
        let mesh_path = self.mesh_path(skeleton_interface);
        let mut palve_path = Vec::new();
        for ind1 in 0..mesh_path.len() {
            let ind2 = (ind1 + 1) % mesh_path.len();
            let ind_vertex1 = mesh_path[ind1];
            let ind_vertex2 = mesh_path[ind2];
            let seg = if ind_vertex1 < ind_vertex2 {
                [ind_vertex1, ind_vertex2]
            } else {
                [ind_vertex2, ind_vertex1]
            };

            let &ind_alveola = skeleton_interface
                .del_seg
                .get(&seg)
                .ok_or(anyhow::Error::msg("Alveola not in skeleton"))?;
            let alve = skeleton_interface.get_alveola_uncheck(ind_alveola);
            let ind_palve = if alve.partial_alveolae()[0].corner() == ind_vertex1 {
                alve.partial_alveolae()[0].ind()
            } else {
                alve.partial_alveolae()[1].ind()
            };
            palve_path.push(ind_palve);
        }

        Ok(palve_path)
    }

    pub fn print(&self, skeleton_interface: &SkeletonInterface3D) -> () {
        for &part in self.components.iter() {
            match part {
                PathPart::PartialEdge(ind_pedge) => {
                    let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge);
                    println!(
                        "({} -> {}): {}",
                        pedge.partial_node_first().unwrap().node().ind(),
                        pedge.partial_node_last().unwrap().node().ind(),
                        pedge.corner()
                    );
                }
                PathPart::PartialNode(ind_pnode) => {
                    let pnode = skeleton_interface.get_partial_node_uncheck(ind_pnode);
                    println!("({}): {}", pnode.node().ind(), pnode.corner());
                }
            };
        }
    }
}
