use std::collections::HashSet;

use anyhow::Result;
use nalgebra::*;

use super::SkeletonInterface3D;
use super::SkeletonSingularPath;

pub enum State {
    Computing,
    Closed,
}

pub struct SkeletonBoundaryPath {
    components: Vec<usize>,
    opt_ind_pedge_first: Option<usize>,
    opt_ind_pedge_last: Option<usize>,
}

impl SkeletonBoundaryPath {
    pub fn create(
        ind_pedge: usize,
        skeleton_interface: &mut SkeletonInterface3D,
    ) -> Result<SkeletonBoundaryPath> {
        let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge);
        if !pedge.is_boundary() {
            return Err(anyhow::Error::msg("Partial edge should be boundary"));
        }
        let edge = pedge.edge();
        if !edge.is_computed() {
            let ind_edge = edge.ind();
            skeleton_interface.propagate_edge(ind_edge)?;
        }
        let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge);

        let ind_pedge_next = pedge.partial_edge_next().unwrap().ind();
        let ind_pedge_prev = pedge.partial_edge_prev().unwrap().ind();

        Ok(SkeletonBoundaryPath {
            components: vec![ind_pedge],
            opt_ind_pedge_first: Some(ind_pedge_prev),
            opt_ind_pedge_last: Some(ind_pedge_next),
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
                .unwrap()
                .ind();

            self.opt_ind_pedge_last = Some(ind_pedge_next);

            self.check_reach_end(&skeleton_interface)
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
                pedge_neigh.partial_edge_next().unwrap()
            } else {
                pedge_neigh.partial_edge_opposite()
            };
            let ind_pedge_new = pedge_next.ind();
            self.opt_ind_pedge_last = Some(ind_pedge_new);
            self.check_reach_end(skeleton_interface)
        } else {
            Ok(State::Closed)
        }
    }

    fn check_reach_end(&mut self, skeleton_interface: &SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_last) = self.opt_ind_pedge_last {
            let ind_pedge_opp = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_last)
                .partial_edge_opposite()
                .ind();

            if let Some(&ind_pedge) = self.components.last() {
                if ind_pedge == ind_pedge_opp {
                    self.opt_ind_pedge_last = None;
                }
            }
        }

        if self.opt_ind_pedge_last.is_none() {
            Ok(State::Closed)
        } else {
            Ok(State::Computing)
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
            self.components.insert(0, ind_pedge_first);

            let ind_pedge_prev = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_first)
                .partial_edge_prev()
                .unwrap()
                .ind();

            self.opt_ind_pedge_first = Some(ind_pedge_prev);

            self.check_reach_beg(&skeleton_interface)
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
                pedge_neigh.partial_edge_prev().unwrap()
            } else {
                pedge_neigh.partial_edge_opposite()
            };
            let ind_pedge_new = pedge_prev.ind();
            self.opt_ind_pedge_first = Some(ind_pedge_new);
            self.check_reach_beg(skeleton_interface)
        } else {
            Ok(State::Closed)
        }
    }

    fn check_reach_beg(&mut self, skeleton_interface: &SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_first) = self.opt_ind_pedge_first {
            let ind_pedge_opp = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_first)
                .partial_edge_opposite()
                .ind();

            if let Some(&ind_pedge) = self.components.first() {
                if ind_pedge_opp == ind_pedge {
                    self.opt_ind_pedge_first = None;
                }
            }
        }

        if self.opt_ind_pedge_first.is_none() {
            Ok(State::Closed)
        } else {
            Ok(State::Computing)
        }
    }

    pub fn ind_last_partial_edge(&self) -> Option<usize> {
        self.opt_ind_pedge_last
    }

    pub fn ind_first_partial_edge(&self) -> Option<usize> {
        self.opt_ind_pedge_first
    }

    pub fn follow_boundary_path(
        &mut self,
        skeleton_interface: &mut SkeletonInterface3D,
    ) -> Result<()> {
        loop {
            if let Some(ind_pedge_last) = self.ind_last_partial_edge() {
                let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge_last);
                if pedge.is_boundary() {
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
                if pedge.is_boundary() {
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

    pub fn most_salient(&self, skeleton_interface: &SkeletonInterface3D) -> Result<usize> {
        todo!();
    }

    pub fn create_path_excluding(
        &self,
        skeleton_interface: &mut SkeletonInterface3D,
        ind_in_path: usize,
    ) -> Result<SkeletonSingularPath> {
        if ind_in_path >= self.components.len() - 1 {
            return Err(anyhow::Error::msg("Index in path out of bounds"));
        }
        let mut ind_first = ind_in_path;
        let mut ind_last = ind_in_path;
        loop {
            if ind_first == 0 {
                break;
            }
            let ind_alve_cur = skeleton_interface
                .get_partial_edge_uncheck(self.components[ind_first])
                .partial_alveola()
                .ind();
            let ind_alve_bef = skeleton_interface
                .get_partial_edge_uncheck(self.components[ind_first - 1])
                .partial_alveola()
                .ind();

            if ind_alve_bef == ind_alve_cur {
                ind_first = ind_first - 1;
            } else {
                break;
            }
        }
        loop {
            if ind_last >= self.components.len() - 1 {
                break;
            }
            let ind_alve_cur = skeleton_interface
                .get_partial_edge_uncheck(self.components[ind_last])
                .partial_alveola()
                .ind();
            let ind_alve_aft = skeleton_interface
                .get_partial_edge_uncheck(self.components[ind_last + 1])
                .partial_alveola()
                .ind();

            if ind_alve_aft == ind_alve_cur {
                ind_last = ind_last + 1;
            } else {
                break;
            }
        }
        let ind_alve1 = skeleton_interface
            .get_partial_edge_uncheck(self.components[ind_first])
            .partial_alveola()
            .alveola()
            .ind();
        let ind_alve2 = skeleton_interface
            .get_partial_edge_uncheck(self.components[ind_last])
            .partial_alveola()
            .alveola()
            .ind();

        let mut edges_path = HashSet::new();
        if ind_alve1 == ind_alve2 {
            let edges1: HashSet<usize> = skeleton_interface
                .get_alveola_uncheck(ind_alve1)
                .edges()
                .iter()
                .filter_map(|&edge| {
                    if edge.degree() != 1 {
                        Some(edge.ind())
                    } else {
                        None
                    }
                })
                .collect();
            for ind_edge in edges1 {
                edges_path.insert(ind_edge);
            }
        } else {
            let edges1: HashSet<usize> = skeleton_interface
                .get_alveola_uncheck(ind_alve1)
                .edges()
                .iter()
                .filter_map(|&edge| {
                    if edge.degree() != 1 {
                        Some(edge.ind())
                    } else {
                        None
                    }
                })
                .collect();
            let edges2: HashSet<usize> = skeleton_interface
                .get_alveola_uncheck(ind_alve2)
                .edges()
                .iter()
                .filter_map(|&edge| {
                    if edge.degree() != 1 {
                        Some(edge.ind())
                    } else {
                        None
                    }
                })
                .collect();
            for &ind_edge in edges1.iter() {
                if !edges2.contains(&ind_edge) {
                    edges_path.insert(ind_edge);
                }
            }
            for &ind_edge in edges2.iter() {
                if !edges1.contains(&ind_edge) {
                    edges_path.insert(ind_edge);
                }
            }
        }

        let ind_pedge_first = skeleton_interface
            .get_partial_edge_uncheck(self.components[ind_last])
            .partial_edge_next()
            .unwrap()
            .ind();

        let mut sing_path = SkeletonSingularPath::create(ind_pedge_first, skeleton_interface);
        sing_path.append_last(skeleton_interface)?;

        loop {
            if let Some(ind_pedge_last) = sing_path.ind_last_partial_edge() {
                let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge_last);
                if edges_path.contains(&pedge.edge().ind()) {
                    sing_path.append_last(skeleton_interface)?;
                } else {
                    sing_path.rotate_last(skeleton_interface)?;
                }
            } else {
                break;
            }
        }

        Ok(sing_path)
    }
}
