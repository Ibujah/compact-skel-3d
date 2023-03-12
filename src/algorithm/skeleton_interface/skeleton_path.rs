use crate::algorithm::skeleton_interface::SkeletonInterface3D;
use anyhow::Result;

use super::skeleton_interface::IterPartialEdge;

#[derive(Copy, Clone)]
pub enum PathPart {
    PartialNode(usize),
    PartialEdge(usize),
}

pub enum State {
    Computing,
    Closed,
}

pub struct SkeletonPath<'a, 'b> {
    skeleton_interface: &'b mut SkeletonInterface3D<'a>,
    components: Vec<PathPart>,
    opt_ind_pedge_last: Option<usize>,
}

impl<'a, 'b> SkeletonPath<'a, 'b> {
    pub fn new(
        skeleton_interface: &'b mut SkeletonInterface3D<'a>,
        ind_pedge: usize,
    ) -> SkeletonPath<'a, 'b> {
        SkeletonPath {
            skeleton_interface,
            components: Vec::new(),
            opt_ind_pedge_last: Some(ind_pedge),
        }
    }

    pub fn skeleton_interface(&self) -> &SkeletonInterface3D {
        self.skeleton_interface
    }

    pub fn mesh_path(&self) -> Vec<usize> {
        self.components
            .iter()
            .map(|part| match part {
                &PathPart::PartialNode(pnode) => self
                    .skeleton_interface
                    .get_partial_node_uncheck(pnode)
                    .corner(),
                &PathPart::PartialEdge(pedge) => self
                    .skeleton_interface
                    .get_partial_edge_uncheck(pedge)
                    .corner(),
            })
            .collect()
    }

    pub fn append_last(&mut self) -> Result<State> {
        if let Some(ind_pedge_last) = self.opt_ind_pedge_last {
            let edge = self
                .skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_last)
                .edge();
            if !edge.is_computed() {
                let ind_edge = edge.ind();
                self.skeleton_interface.propagate_edge(ind_edge)?;
            }
            self.components.push(PathPart::PartialEdge(ind_pedge_last));

            let ind_pedge_next = self
                .skeleton_interface
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

    pub fn rotate_last(&mut self) -> Result<State> {
        if let Some(ind_pedge_last) = self.opt_ind_pedge_last {
            let pedge = self
                .skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_last);
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
            let ind_pedge_new = pedge_next.ind();
            self.components.push(PathPart::PartialNode(ind_pnode));
            self.opt_ind_pedge_last = Some(ind_pedge_new);
            self.check_loop()
        } else {
            Ok(State::Closed)
        }
    }

    pub fn ind_last_partial_edge(&self) -> Option<usize> {
        self.opt_ind_pedge_last
    }

    pub fn last_partial_edge(&self) -> Option<IterPartialEdge> {
        if let Some(ind_pedge_last) = self.opt_ind_pedge_last {
            Some(
                self.skeleton_interface
                    .get_partial_edge_uncheck(ind_pedge_last),
            )
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

    pub fn ind_alveolae(&self) -> Vec<usize> {
        self.components
            .iter()
            .filter_map(|pp| match pp {
                &PathPart::PartialEdge(ind_pedge) => Some(
                    self.skeleton_interface
                        .get_partial_edge_uncheck(ind_pedge)
                        .partial_alveola()
                        .alveola()
                        .ind(),
                ),
                &PathPart::PartialNode(_) => None,
            })
            .collect()
    }

    pub fn close_path(&mut self) -> Result<()> {
        for ind in 0..self.components.len() {
            let ind_next = (ind + 1) % self.components.len();
            match (self.components[ind], self.components[ind_next]) {
                (PathPart::PartialEdge(ind_pedge), PathPart::PartialNode(_)) => {
                    let pedge = self.skeleton_interface.get_partial_edge_uncheck(ind_pedge);
                    let [ind_vertex1, ind_vertex2, ind_vertex3] = pedge
                        .partial_edge_next()
                        .ok_or(anyhow::Error::msg("Next edge does not exist"))?
                        .edge()
                        .delaunay_triangle();
                    self.skeleton_interface
                        .close_face(ind_vertex1, ind_vertex2, ind_vertex3)?
                }
                (_, _) => (),
            };
        }
        Ok(())
    }
}
