use crate::algorithm::skeleton_interface::SkeletonInterface3D;
use anyhow::Result;

#[derive(Copy, Clone)]
pub enum PathPart {
    PartialNode(usize),
    PartialEdge(usize),
}

pub enum State {
    Computing,
    Closed,
}

pub struct SkeletonPath<'a, 'b, 'c> {
    skeleton_interface: &'c mut SkeletonInterface3D<'a, 'b>,
    components: Vec<PathPart>,
    opt_ind_pedge_last: Option<usize>,
}

impl<'a, 'b, 'c> SkeletonPath<'a, 'b, 'c> {
    pub fn new(
        skeleton_interface: &'c mut SkeletonInterface3D<'a, 'b>,
        ind_pedge: usize,
    ) -> SkeletonPath<'a, 'b, 'c> {
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
            .map(|part| part.corner(self.skeleton_interface))
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

    pub fn partial_edges(&self) -> Vec<usize> {
        self.components
            .iter()
            .filter_map(|pp| match pp {
                &PathPart::PartialEdge(ind_pedge) => Some(ind_pedge),
                &PathPart::PartialNode(_) => None,
            })
            .collect()
    }

    pub fn alveolae(&self) -> Vec<usize> {
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
}

impl PathPart {
    pub fn corner(&self, skeleton_interface: &SkeletonInterface3D) -> usize {
        match self {
            &PathPart::PartialNode(pnode) => {
                skeleton_interface.get_partial_node_uncheck(pnode).corner()
            }
            &PathPart::PartialEdge(pedge) => {
                skeleton_interface.get_partial_edge_uncheck(pedge).corner()
            }
        }
    }
}
