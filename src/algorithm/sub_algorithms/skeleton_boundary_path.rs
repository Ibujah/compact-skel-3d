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
}
