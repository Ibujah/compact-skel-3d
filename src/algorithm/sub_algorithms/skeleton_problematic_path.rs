use anyhow::Result;

use super::SkeletonInterface3D;

pub enum State {
    Computing,
    Closed,
}

pub struct SkeletonProblematicPath {
    components: Vec<usize>,
    opt_ind_pedge_first: Option<usize>,
    opt_ind_pedge_last: Option<usize>,
}

impl SkeletonProblematicPath {
    pub fn create(
        ind_pedge: usize,
        skeleton_interface: &mut SkeletonInterface3D,
    ) -> Result<SkeletonProblematicPath> {
        let edge = skeleton_interface
            .get_partial_edge_uncheck(ind_pedge)
            .edge();
        if !edge.is_computed() {
            let ind_edge = edge.ind();
            skeleton_interface.propagate_edge(ind_edge)?;
        }
        let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge);
        let ind_prev = pedge.partial_edge_prev().unwrap().ind();
        let ind_next = pedge.partial_edge_next().unwrap().ind();
        Ok(SkeletonProblematicPath {
            components: vec![ind_pedge],
            opt_ind_pedge_first: Some(ind_prev),
            opt_ind_pedge_last: Some(ind_next),
        })
    }

    pub fn components(&self) -> &Vec<usize> {
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
            self.components.push(ind_pedge_last);

            let ind_pedge_next = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_last)
                .partial_edge_next()
                .ok_or(anyhow::Error::msg("No next partial edge"))?
                .ind();

            self.opt_ind_pedge_last = Some(ind_pedge_next);
            self.check_end(skeleton_interface)
        } else {
            Ok(State::Closed)
        }
    }

    fn check_end(&mut self, skeleton_interface: &SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_last) = self.opt_ind_pedge_last {
            let edge_last = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_last)
                .edge();
            let ind_pedg = *self.components.last().unwrap();
            let edge = skeleton_interface.get_partial_edge_uncheck(ind_pedg).edge();

            if edge_last.ind() == edge.ind() {
                self.opt_ind_pedge_last = None;
                Ok(State::Closed)
            } else {
                Ok(State::Computing)
            }
        } else {
            Ok(State::Closed)
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
            let ind_pedge_new = pedge_next.ind();
            self.opt_ind_pedge_last = Some(ind_pedge_new);
            self.check_end(skeleton_interface)
        } else {
            Ok(State::Closed)
        }
    }

    pub fn ind_last_partial_edge(&self) -> Option<usize> {
        self.opt_ind_pedge_last
    }

    pub fn append_first(&mut self, skeleton_interface: &mut SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_first) = self.opt_ind_pedge_first {
            let edge = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_first)
                .edge();
            if !edge.is_computed() {
                let ind_edge = edge.ind();
                skeleton_interface.propagate_edge(ind_edge)?;
            }
            self.components.push(ind_pedge_first);

            let ind_pedge_prev = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_first)
                .partial_edge_prev()
                .ok_or(anyhow::Error::msg("No prev partial edge"))?
                .ind();

            self.opt_ind_pedge_first = Some(ind_pedge_prev);
            self.check_end(skeleton_interface)
        } else {
            Ok(State::Closed)
        }
    }

    fn check_beg(&mut self, skeleton_interface: &SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_first) = self.opt_ind_pedge_first {
            let edge_first = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_first)
                .edge();
            let ind_pedg = *self.components.first().unwrap();
            let edge = skeleton_interface.get_partial_edge_uncheck(ind_pedg).edge();

            if edge_first.ind() == edge.ind() {
                self.opt_ind_pedge_first = None;
                Ok(State::Closed)
            } else {
                Ok(State::Computing)
            }
        } else {
            Ok(State::Closed)
        }
    }

    pub fn rotate_first(&mut self, skeleton_interface: &SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_first) = self.opt_ind_pedge_first {
            let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge_first);
            let pedge_neigh = pedge.partial_edge_neighbor();

            let alve = pedge_neigh.partial_alveola().alveola();
            let pedge_prev = if alve.is_full() {
                pedge_neigh
                    .partial_edge_prev()
                    .ok_or(anyhow::Error::msg("No next partial edge"))?
            } else {
                pedge_neigh.partial_edge_opposite()
            };
            let ind_pedge_new = pedge_prev.ind();
            self.opt_ind_pedge_first = Some(ind_pedge_new);
            self.check_beg(skeleton_interface)
        } else {
            Ok(State::Closed)
        }
    }

    pub fn ind_first_partial_edge(&self) -> Option<usize> {
        self.opt_ind_pedge_first
    }

    pub fn follow_problematic_path(
        &mut self,
        skeleton_interface: &mut SkeletonInterface3D,
    ) -> Result<()> {
        loop {
            if let Some(ind_pedge_last) = self.ind_last_partial_edge() {
                let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge_last);
                if pedge.edge().is_non_manifold() {
                    self.append_last(skeleton_interface)?;
                } else {
                    self.rotate_last(skeleton_interface)?;
                }
            } else {
                break;
            }
        }
        loop {
            if let Some(ind_pedge_first) = self.ind_first_partial_edge() {
                let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge_first);
                if pedge.edge().is_non_manifold() {
                    self.append_first(skeleton_interface)?;
                } else {
                    self.rotate_first(skeleton_interface)?;
                }
            } else {
                break;
            }
        }
        Ok(())
    }

    pub fn removable_from_first(
        self: &SkeletonProblematicPath,
        skeleton_interface: &SkeletonInterface3D,
    ) -> Result<bool> {
        let ind_pedge_last = self.components.last().unwrap();
        let ind_pedge_first = self.components.first().unwrap();

        let pedge_last = skeleton_interface.get_partial_edge(*ind_pedge_last)?;
        let pedge_first = skeleton_interface.get_partial_edge(*ind_pedge_first)?;

        if pedge_last.partial_edge_next().unwrap().edge().degree() == 1
            && pedge_first.partial_edge_prev().unwrap().edge().degree() == 2
        {
            Ok(true)
        } else {
            Ok(false)
        }
    }

    pub fn removable_from_last(
        self: &SkeletonProblematicPath,
        skeleton_interface: &SkeletonInterface3D,
    ) -> Result<bool> {
        let ind_pedge_last = self.components.last().unwrap();
        let ind_pedge_first = self.components.first().unwrap();

        let pedge_last = skeleton_interface.get_partial_edge(*ind_pedge_last)?;
        let pedge_first = skeleton_interface.get_partial_edge(*ind_pedge_first)?;

        if pedge_last.partial_edge_next().unwrap().edge().degree() == 2
            && pedge_first.partial_edge_prev().unwrap().edge().degree() == 1
        {
            Ok(true)
        } else {
            Ok(false)
        }
    }

    // pub fn print(&self, skeleton_interface: &SkeletonInterface3D) -> () {
    //     for &part in self.components.iter() {
    //         match part {
    //             PathPart::PartialEdge(ind_pedge) => {
    //                 let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge);
    //                 println!(
    //                     "({} -> {}): {}",
    //                     pedge.partial_node_first().unwrap().node().ind(),
    //                     pedge.partial_node_last().unwrap().node().ind(),
    //                     pedge.corner()
    //                 );
    //             }
    //             PathPart::PartialNode(ind_pnode) => {
    //                 let pnode = skeleton_interface.get_partial_node_uncheck(ind_pnode);
    //                 println!("({}): {}", pnode.node().ind(), pnode.corner());
    //             }
    //         };
    //     }
    // }
}

pub fn last_to_boundary(
    skel_prob: &SkeletonProblematicPath,
    skeleton_interface: &mut SkeletonInterface3D,
) -> Result<Vec<usize>> {
    let pedge_last =
        skeleton_interface.get_partial_edge(*skel_prob.components().last().unwrap())?;
    let pedge_path = pedge_last.partial_edge_next().unwrap();
    let mut ind_pedge_path = pedge_path.ind();
    let direction = (pedge_path
        .partial_node_last()
        .unwrap()
        .node()
        .center_and_radius()?
        .0
        - pedge_path
            .partial_node_first()
            .unwrap()
            .node()
            .center_and_radius()?
            .0)
        .normalize();

    let mut vec_path = vec![ind_pedge_path];

    loop {
        let pedge_path = skeleton_interface.get_partial_edge(ind_pedge_path)?;
        let pnode_path_next = pedge_path.partial_node_last().unwrap();

        let mut reached_bnd = false;
        for pedge in pnode_path_next.partial_edge_next() {
            if pedge.edge().degree() == 1 {
                reached_bnd = true;
                break;
            }
        }
        if reached_bnd {
            break;
        }

        let (opt_ind_pedge_next, _) = pnode_path_next
            .partial_edge_next()
            .iter()
            .map(|pedg| {
                let dir = pedg
                    .partial_node_last()
                    .unwrap()
                    .node()
                    .center_and_radius()
                    .unwrap()
                    .0
                    - pedg
                        .partial_node_first()
                        .unwrap()
                        .node()
                        .center_and_radius()
                        .unwrap()
                        .0;
                (pedg.ind(), direction.dot(&dir))
            })
            .fold((None, 0.0), |(imax, cmax), (icur, ccur)| {
                if imax.is_some() {
                    if ccur > cmax {
                        (Some(icur), ccur)
                    } else {
                        (imax, cmax)
                    }
                } else {
                    (Some(icur), ccur)
                }
            });

        ind_pedge_path = opt_ind_pedge_next.unwrap();
        vec_path.push(ind_pedge_path);
    }

    Ok(vec_path)
}

pub fn first_to_boundary(
    skel_prob: &SkeletonProblematicPath,
    skeleton_interface: &mut SkeletonInterface3D,
) -> Result<Vec<usize>> {
    let pedge_first =
        skeleton_interface.get_partial_edge(*skel_prob.components().first().unwrap())?;
    let pedge_path = pedge_first.partial_edge_prev().unwrap();
    let mut ind_pedge_path = pedge_path.ind();
    let direction = (pedge_path
        .partial_node_first()
        .unwrap()
        .node()
        .center_and_radius()?
        .0
        - pedge_path
            .partial_node_last()
            .unwrap()
            .node()
            .center_and_radius()?
            .0)
        .normalize();

    let mut vec_path = vec![ind_pedge_path];

    loop {
        let pedge_path = skeleton_interface.get_partial_edge(ind_pedge_path)?;
        let pnode_path_prev = pedge_path.partial_node_first().unwrap();

        let mut reached_bnd = false;
        for pedge in pnode_path_prev.partial_edge_prev() {
            if pedge.edge().degree() == 1 {
                reached_bnd = true;
                break;
            }
        }
        if reached_bnd {
            break;
        }

        let (opt_ind_pedge_prev, _) = pnode_path_prev
            .partial_edge_prev()
            .iter()
            .map(|pedg| {
                let dir = pedg
                    .partial_node_first()
                    .unwrap()
                    .node()
                    .center_and_radius()
                    .unwrap()
                    .0
                    - pedg
                        .partial_node_last()
                        .unwrap()
                        .node()
                        .center_and_radius()
                        .unwrap()
                        .0;
                (pedg.ind(), direction.dot(&dir))
            })
            .fold((None, 0.0), |(imax, cmax), (icur, ccur)| {
                if imax.is_some() {
                    if ccur > cmax {
                        (Some(icur), ccur)
                    } else {
                        (imax, cmax)
                    }
                } else {
                    (Some(icur), ccur)
                }
            });

        ind_pedge_path = opt_ind_pedge_prev.unwrap();
        vec_path.push(ind_pedge_path);
    }

    Ok(vec_path)
}
