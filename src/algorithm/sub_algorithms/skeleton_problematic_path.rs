use std::collections::HashMap;

use anyhow::Result;
use nalgebra::base::*;

use super::SkeletonInterface3D;

pub enum State {
    Computing,
    Closed,
}

pub struct SkeletonProblematicPath {
    components_non_manifold: Vec<usize>,
    components_boundary: Vec<usize>,
    opt_ind_pedge_first: Option<usize>,
    opt_ind_pedge_last: Option<usize>,
    opt_ind_pedge_before_first: Option<usize>,
    opt_ind_pedge_after_last: Option<usize>,
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
            components_non_manifold: vec![ind_pedge],
            components_boundary: Vec::new(),
            opt_ind_pedge_first: Some(ind_prev),
            opt_ind_pedge_last: Some(ind_next),
            opt_ind_pedge_before_first: None,
            opt_ind_pedge_after_last: None,
        })
    }

    pub fn components(&self) -> &Vec<usize> {
        &self.components_non_manifold
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
            self.components_non_manifold.push(ind_pedge_last);

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

    pub fn append_first(&mut self, skeleton_interface: &mut SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_first) = self.opt_ind_pedge_first {
            let edge = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_first)
                .edge();
            if !edge.is_computed() {
                let ind_edge = edge.ind();
                skeleton_interface.propagate_edge(ind_edge)?;
            }
            self.components_non_manifold.insert(0, ind_pedge_first);

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

    fn check_end(&mut self, skeleton_interface: &SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_last) = self.opt_ind_pedge_last {
            let edge_last = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_last)
                .edge();
            let ind_pedg_beg = *self.components_non_manifold.first().unwrap();
            let edge_beg = skeleton_interface
                .get_partial_edge_uncheck(ind_pedg_beg)
                .edge();
            let ind_pedg_end = *self.components_non_manifold.last().unwrap();
            let edge_end = skeleton_interface
                .get_partial_edge_uncheck(ind_pedg_end)
                .edge();

            if edge_last.ind() == edge_end.ind() {
                self.opt_ind_pedge_last = None;
                Ok(State::Closed)
            } else if edge_last.ind() == edge_beg.ind() {
                self.opt_ind_pedge_first = None;
                self.opt_ind_pedge_last = None;
                Ok(State::Closed)
            } else {
                Ok(State::Computing)
            }
        } else {
            Ok(State::Closed)
        }
    }

    fn check_beg(&mut self, skeleton_interface: &SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_first) = self.opt_ind_pedge_first {
            let edge_first = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_first)
                .edge();
            let ind_pedg_beg = *self.components_non_manifold.first().unwrap();
            let edge_beg = skeleton_interface
                .get_partial_edge_uncheck(ind_pedg_beg)
                .edge();
            let ind_pedg_end = *self.components_non_manifold.last().unwrap();
            let edge_end = skeleton_interface
                .get_partial_edge_uncheck(ind_pedg_end)
                .edge();

            if edge_first.ind() == edge_beg.ind() {
                self.opt_ind_pedge_first = None;
                Ok(State::Closed)
            } else if edge_first.ind() == edge_end.ind() {
                self.opt_ind_pedge_first = None;
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

    pub fn ind_last_partial_edge(&self) -> Option<usize> {
        self.opt_ind_pedge_last
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

    pub fn append_after_last(
        &mut self,
        skeleton_interface: &mut SkeletonInterface3D,
    ) -> Result<State> {
        if let Some(ind_pedge_after_last) = self.opt_ind_pedge_after_last {
            let edge = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_after_last)
                .edge();
            if !edge.is_computed() {
                let ind_edge = edge.ind();
                skeleton_interface.propagate_edge(ind_edge)?;
            }
            self.components_boundary.push(ind_pedge_after_last);

            let ind_pedge_next = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_after_last)
                .partial_edge_next()
                .ok_or(anyhow::Error::msg("No next partial edge"))?
                .ind();

            self.opt_ind_pedge_after_last = Some(ind_pedge_next);
            self.check_end(skeleton_interface)
        } else {
            Ok(State::Closed)
        }
    }

    pub fn append_before_first(
        &mut self,
        skeleton_interface: &mut SkeletonInterface3D,
    ) -> Result<State> {
        if let Some(ind_pedge_before_first) = self.opt_ind_pedge_before_first {
            let edge = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_before_first)
                .edge();
            if !edge.is_computed() {
                let ind_edge = edge.ind();
                skeleton_interface.propagate_edge(ind_edge)?;
            }
            self.components_boundary.insert(0, ind_pedge_before_first);

            let ind_pedge_prev = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_before_first)
                .partial_edge_prev()
                .ok_or(anyhow::Error::msg("No prev partial edge"))?
                .ind();

            self.opt_ind_pedge_before_first = Some(ind_pedge_prev);
            self.check_end(skeleton_interface)
        } else {
            Ok(State::Closed)
        }
    }

    fn check_after_end(&mut self, skeleton_interface: &SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_after_last) = self.opt_ind_pedge_after_last {
            if self.components_boundary.len() == 0 {
                return Ok(State::Computing);
            }
            let edge_last = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_after_last)
                .edge();
            let ind_pedg_beg = *self.components_boundary.first().unwrap();
            let edge_beg = skeleton_interface
                .get_partial_edge_uncheck(ind_pedg_beg)
                .edge();
            let ind_pedg_end = *self.components_boundary.last().unwrap();
            let edge_end = skeleton_interface
                .get_partial_edge_uncheck(ind_pedg_end)
                .edge();

            if edge_last.ind() == edge_end.ind() {
                self.opt_ind_pedge_last = None;
                Ok(State::Closed)
            } else if edge_last.ind() == edge_beg.ind() {
                self.opt_ind_pedge_first = None;
                self.opt_ind_pedge_last = None;
                Ok(State::Closed)
            } else {
                Ok(State::Computing)
            }
        } else {
            Ok(State::Closed)
        }
    }

    fn check_before_beg(&mut self, skeleton_interface: &SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_before_first) = self.opt_ind_pedge_before_first {
            let edge_first = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_before_first)
                .edge();
            let ind_pedg_beg = *self.components_non_manifold.first().unwrap();
            let edge_beg = skeleton_interface
                .get_partial_edge_uncheck(ind_pedg_beg)
                .edge();
            let ind_pedg_end = *self.components_non_manifold.last().unwrap();
            let edge_end = skeleton_interface
                .get_partial_edge_uncheck(ind_pedg_end)
                .edge();

            if edge_first.ind() == edge_beg.ind() {
                self.opt_ind_pedge_first = None;
                Ok(State::Closed)
            } else if edge_first.ind() == edge_end.ind() {
                self.opt_ind_pedge_first = None;
                self.opt_ind_pedge_last = None;
                Ok(State::Closed)
            } else {
                Ok(State::Computing)
            }
        } else {
            Ok(State::Closed)
        }
    }

    pub fn rotate_after_last(&mut self, skeleton_interface: &SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_after_last) = self.opt_ind_pedge_after_last {
            let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge_after_last);
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
            self.opt_ind_pedge_after_last = Some(ind_pedge_new);
            self.check_after_end(skeleton_interface)
        } else {
            Ok(State::Closed)
        }
    }

    pub fn rotate_before_first(
        &mut self,
        skeleton_interface: &SkeletonInterface3D,
    ) -> Result<State> {
        if let Some(ind_pedge_before_first) = self.opt_ind_pedge_before_first {
            let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge_before_first);
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
            self.opt_ind_pedge_before_first = Some(ind_pedge_new);
            self.check_before_beg(skeleton_interface)
        } else {
            Ok(State::Closed)
        }
    }

    pub fn follow_boundary_path_from_first(
        self: &mut SkeletonProblematicPath,
        skeleton_interface: &mut SkeletonInterface3D,
    ) -> Result<bool> {
        let ind_pedge_first = self.components_non_manifold.first().unwrap();
        let pedge_first = skeleton_interface
            .get_partial_edge(*ind_pedge_first)?
            .partial_edge_prev()
            .unwrap();

        self.opt_ind_pedge_before_first = Some(pedge_first.ind());

        loop {
            if let Some(ind_before_first) = self.opt_ind_pedge_before_first {
                let pedge_first = skeleton_interface.get_partial_edge(ind_before_first)?;
                if pedge_first.edge().is_boundary() {
                    self.append_before_first(skeleton_interface)?;
                } else if pedge_first.edge().is_singular() {
                    self.opt_ind_pedge_before_first = None;
                    break;
                } else {
                    self.rotate_before_first(skeleton_interface)?;
                }
            } else {
                break;
            }
        }
        Ok(self.components_boundary.len() != 0)
    }

    pub fn follow_boundary_path_from_last(
        self: &mut SkeletonProblematicPath,
        skeleton_interface: &mut SkeletonInterface3D,
    ) -> Result<bool> {
        let ind_pedge_last = self.components_non_manifold.last().unwrap();
        let pedge_last = skeleton_interface
            .get_partial_edge(*ind_pedge_last)?
            .partial_edge_next()
            .unwrap();

        self.opt_ind_pedge_after_last = Some(pedge_last.ind());

        loop {
            if let Some(ind_after_last) = self.opt_ind_pedge_after_last {
                let pedge_last = skeleton_interface.get_partial_edge(ind_after_last)?;
                if pedge_last.edge().is_boundary() {
                    self.append_after_last(skeleton_interface)?;
                } else if pedge_last.edge().is_singular() {
                    self.opt_ind_pedge_after_last = None;
                    break;
                } else {
                    self.rotate_after_last(skeleton_interface)?;
                }
            } else {
                break;
            }
        }
        Ok(self.components_boundary.len() != 0)
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

fn next_pedges_to_eval(
    ind_node: usize,
    skeleton_interface: &SkeletonInterface3D,
    label: usize,
) -> Vec<usize> {
    let mut vec_res = Vec::new();
    for edg in skeleton_interface.get_node_uncheck(ind_node).edges().iter() {
        if !edg.is_full() {
            continue;
        }
        for pedg in edg.partial_edges().iter() {
            if let Some(pnode) = pedg.partial_node_first() {
                if pnode.node().ind() != ind_node {
                    continue;
                }
            } else {
                continue;
            }
            if pedg.partial_alveola().alveola().label() == Some(label) {
                vec_res.push(pedg.ind());
            }
        }
    }
    vec_res
}

fn dist_min(ctr: &Vector3<f32>, vec_centers: &Vec<Vector3<f32>>) -> f32 {
    vec_centers
        .iter()
        .map(|ctr_cur| (ctr - ctr_cur).norm())
        .fold(None, |opt_dmin, dcur| {
            if let Some(dmin) = opt_dmin {
                if dcur < dmin {
                    Some(dcur)
                } else {
                    opt_dmin
                }
            } else {
                Some(dcur)
            }
        })
        .unwrap()
}

pub(super) fn last_to_boundary(
    skel_prob: &SkeletonProblematicPath,
    skeleton_interface: &mut SkeletonInterface3D,
    label: usize,
) -> Result<Vec<usize>> {
    let vec_centers: Vec<Vector3<f32>> = skel_prob
        .components_boundary
        .iter()
        .map(|&ind_pedge| {
            let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge);
            let node = pedge.partial_node_first().unwrap().node();
            node.center_and_radius().unwrap().0
        })
        .collect();

    let mut map_nodes_dist = HashMap::new();
    let mut map_nodes_prev = HashMap::new();
    let mut map_nodes_next_to_eval = HashMap::new();
    let mut map_nodes_ctr = HashMap::new();
    let mut map_nodes_dist_to_bnd = HashMap::new();

    let pedge_last =
        skeleton_interface.get_partial_edge(*skel_prob.components().last().unwrap())?;

    let node_init = pedge_last.partial_node_last().unwrap().node();

    map_nodes_dist.insert(node_init.ind(), 0.0);
    map_nodes_prev.insert(node_init.ind(), None);
    map_nodes_next_to_eval.insert(
        node_init.ind(),
        next_pedges_to_eval(node_init.ind(), skeleton_interface, label),
    );
    let ctr_init = node_init.center_and_radius().unwrap().0;
    map_nodes_ctr.insert(node_init.ind(), ctr_init);
    map_nodes_dist_to_bnd.insert(node_init.ind(), dist_min(&ctr_init, &vec_centers));

    let mut opt_last_node = None;

    loop {
        let opt_min = map_nodes_dist_to_bnd
            .iter()
            .fold(None, |opt_min, (icur, dcur)| {
                if let Some((_, dmin)) = opt_min {
                    if dcur < dmin {
                        Some((icur, dcur))
                    } else {
                        opt_min
                    }
                } else {
                    Some((icur, dcur))
                }
            });
        if let Some((&ind_node_cur, _)) = opt_min {
            let node = skeleton_interface.get_node_uncheck(ind_node_cur);
            if node.ind() != node_init.ind() {
                if node.edges().iter().any(|edg| edg.is_boundary()) {
                    opt_last_node = Some(ind_node_cur);
                    break;
                }
            }
            map_nodes_dist_to_bnd.remove(&ind_node_cur);
            let to_eval = map_nodes_next_to_eval.remove(&ind_node_cur).unwrap();
            let &ctr_cur = map_nodes_ctr.get(&ind_node_cur).unwrap();
            let &dist_cur = map_nodes_dist.get(&ind_node_cur).unwrap();
            for &ind_pedge in to_eval.iter() {
                let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge);
                let node_aft = pedge.partial_node_last().unwrap().node();
                let ind_node_aft = node_aft.ind();
                let ctr_aft = node_aft.center_and_radius().unwrap().0;
                let dist_aft = dist_cur + (ctr_aft - ctr_cur).norm();
                let should_add = if let Some(&dist) = map_nodes_dist.get(&ind_node_aft) {
                    dist_aft < dist
                } else {
                    true
                };
                if should_add {
                    map_nodes_dist.insert(ind_node_aft, dist_aft);
                    map_nodes_prev.insert(ind_node_aft, Some(ind_node_cur));
                    map_nodes_next_to_eval.insert(
                        ind_node_aft,
                        next_pedges_to_eval(ind_node_aft, skeleton_interface, label),
                    );
                    map_nodes_ctr.insert(ind_node_aft, ctr_aft);
                    map_nodes_dist_to_bnd.insert(ind_node_aft, dist_min(&ctr_aft, &vec_centers));
                }
            }
        } else {
            break;
        }
    }

    let mut vec_edges = Vec::new();

    loop {
        if let Some(last_node) = opt_last_node {
            let opt_ind_prev = map_nodes_prev
                .get(&last_node)
                .ok_or(anyhow::Error::msg("No previous node"))?;
            if let Some(ind_prev) = opt_ind_prev {
                let node = skeleton_interface.get_node_uncheck(last_node);
                for edge in node.edges() {
                    let vec_nod = edge.nodes();
                    if vec_nod.len() == 2 {
                        if vec_nod[0].ind() == *ind_prev || vec_nod[1].ind() == *ind_prev {
                            vec_edges.push(edge.ind());
                            break;
                        }
                    }
                }
                opt_last_node = Some(*ind_prev);
            } else {
                break;
            }
        } else {
            break;
        }
    }

    Ok(vec_edges)
}

pub(super) fn first_to_boundary(
    skel_prob: &SkeletonProblematicPath,
    skeleton_interface: &mut SkeletonInterface3D,
    label: usize,
) -> Result<Vec<usize>> {
    let vec_centers: Vec<Vector3<f32>> = skel_prob
        .components_boundary
        .iter()
        .map(|&ind_pedge| {
            let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge);
            let node = pedge.partial_node_last().unwrap().node();
            node.center_and_radius().unwrap().0
        })
        .collect();

    let mut map_nodes_dist = HashMap::new();
    let mut map_nodes_prev = HashMap::new();
    let mut map_nodes_next_to_eval = HashMap::new();
    let mut map_nodes_ctr = HashMap::new();
    let mut map_nodes_dist_to_bnd = HashMap::new();

    let pedge_first =
        skeleton_interface.get_partial_edge(*skel_prob.components().first().unwrap())?;

    let node_init = pedge_first.partial_node_first().unwrap().node();

    map_nodes_dist.insert(node_init.ind(), 0.0);
    map_nodes_prev.insert(node_init.ind(), None);
    map_nodes_next_to_eval.insert(
        node_init.ind(),
        next_pedges_to_eval(node_init.ind(), skeleton_interface, label),
    );
    let ctr_init = node_init.center_and_radius().unwrap().0;
    map_nodes_ctr.insert(node_init.ind(), ctr_init);
    map_nodes_dist_to_bnd.insert(node_init.ind(), dist_min(&ctr_init, &vec_centers));

    let mut opt_first_node = None;

    loop {
        let opt_min = map_nodes_dist_to_bnd
            .iter()
            .fold(None, |opt_min, (icur, dcur)| {
                if let Some((_, dmin)) = opt_min {
                    if dcur < dmin {
                        Some((icur, dcur))
                    } else {
                        opt_min
                    }
                } else {
                    Some((icur, dcur))
                }
            });
        if let Some((&ind_node_cur, _)) = opt_min {
            let node = skeleton_interface.get_node_uncheck(ind_node_cur);
            if node.ind() != node_init.ind() {
                if node.edges().iter().any(|edg| edg.is_boundary()) {
                    opt_first_node = Some(ind_node_cur);
                    break;
                }
            }
            map_nodes_dist_to_bnd.remove(&ind_node_cur);
            let to_eval = map_nodes_next_to_eval.remove(&ind_node_cur).unwrap();
            let &ctr_cur = map_nodes_ctr.get(&ind_node_cur).unwrap();
            let &dist_cur = map_nodes_dist.get(&ind_node_cur).unwrap();
            for &ind_pedge in to_eval.iter() {
                let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge);
                let node_aft = pedge.partial_node_last().unwrap().node();
                let ind_node_aft = node_aft.ind();
                let ctr_aft = node_aft.center_and_radius().unwrap().0;
                let dist_aft = dist_cur + (ctr_aft - ctr_cur).norm();
                let should_add = if let Some(&dist) = map_nodes_dist.get(&ind_node_aft) {
                    dist_aft < dist
                } else {
                    true
                };
                if should_add {
                    map_nodes_dist.insert(ind_node_aft, dist_aft);
                    map_nodes_prev.insert(ind_node_aft, Some(ind_node_cur));
                    map_nodes_next_to_eval.insert(
                        ind_node_aft,
                        next_pedges_to_eval(ind_node_aft, skeleton_interface, label),
                    );
                    map_nodes_ctr.insert(ind_node_aft, ctr_aft);
                    map_nodes_dist_to_bnd.insert(ind_node_aft, dist_min(&ctr_aft, &vec_centers));
                }
            }
        } else {
            break;
        }
    }

    let mut vec_edges = Vec::new();

    loop {
        if let Some(last_node) = opt_first_node {
            let opt_ind_prev = map_nodes_prev
                .get(&last_node)
                .ok_or(anyhow::Error::msg("No previous node"))?;
            if let Some(ind_prev) = opt_ind_prev {
                let node = skeleton_interface.get_node_uncheck(last_node);
                for edge in node.edges() {
                    let vec_nod = edge.nodes();
                    if vec_nod.len() == 2 {
                        if vec_nod[0].ind() == *ind_prev || vec_nod[1].ind() == *ind_prev {
                            vec_edges.push(edge.ind());
                            break;
                        }
                    }
                }
                opt_first_node = Some(*ind_prev);
            } else {
                break;
            }
        } else {
            break;
        }
    }

    Ok(vec_edges)
}
