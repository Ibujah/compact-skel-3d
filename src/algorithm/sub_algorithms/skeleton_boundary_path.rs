use anyhow::Result;
use nalgebra::*;

use super::path_part::PathPart;
use super::SkeletonInterface3D;

pub enum State {
    Computing,
    Closed,
}

pub struct SkeletonBoundaryPath {
    components: Vec<PathPart>,
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
            components: vec![PathPart::PartialEdge(ind_pedge)],
            opt_ind_pedge_first: Some(ind_pedge_prev),
            opt_ind_pedge_last: Some(ind_pedge_next),
        })
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

            self.check_reach_end(&skeleton_interface)
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

            if self.components.iter().any(|&path_part| {
                if let PathPart::PartialEdge(ind_pedge) = path_part {
                    ind_pedge_opp == ind_pedge
                } else {
                    false
                }
            }) {
                self.opt_ind_pedge_last = None;
            }
        }

        if self.opt_ind_pedge_last.is_none() {
            Ok(State::Closed)
        } else {
            Ok(State::Computing)
        }
    }

    fn check_reach_beg(&mut self, skeleton_interface: &SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_first) = self.opt_ind_pedge_first {
            let ind_pedge_opp = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_first)
                .partial_edge_opposite()
                .ind();

            if self.components.iter().any(|&path_part| {
                if let PathPart::PartialEdge(ind_pedge) = path_part {
                    ind_pedge_opp == ind_pedge
                } else {
                    false
                }
            }) {
                self.opt_ind_pedge_first = None;
            }
        }

        if self.opt_ind_pedge_first.is_none() {
            Ok(State::Closed)
        } else {
            Ok(State::Computing)
        }
    }
}
