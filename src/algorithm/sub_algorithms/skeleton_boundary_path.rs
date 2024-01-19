use std::collections::HashSet;

use anyhow::Result;

use super::SkeletonInterface3D;
use super::SkeletonSingularPath;

pub fn next_boundary_edge(
    ind_pedge: usize,
    skeleton_interface: &SkeletonInterface3D,
) -> Option<usize> {
    let pedge_cur = skeleton_interface.get_partial_edge_uncheck(ind_pedge);
    if pedge_cur.edge().degree() != 1 {
        return None;
    }
    if pedge_cur.partial_alveola().alveola().label().is_none() {
        return None;
    }
    let pnode_last = pedge_cur.partial_node_last().unwrap();
    let pedge_next: Vec<usize> = pnode_last
        .partial_edge_next()
        .iter()
        .filter_map(|pedge| {
            if pedge.partial_alveola().alveola().label().is_some()
                && pedge.edge().degree() == 1
                && pedge.edge().ind() != pedge_cur.edge().ind()
            {
                Some(pedge.ind())
            } else {
                None
            }
        })
        .collect();

    if pedge_next.len() != 1 {
        return None;
    }

    Some(pedge_next[0])
}

pub fn compute_saliency(
    ind_pedge: usize,
    skeleton_interface: &SkeletonInterface3D,
) -> Result<Option<f64>> {
    if let Some(ind_aft) = next_boundary_edge(ind_pedge, skeleton_interface) {
        let pedge_cur = skeleton_interface.get_partial_edge_uncheck(ind_pedge);
        let pedge_aft = skeleton_interface.get_partial_edge_uncheck(ind_aft);

        let pedge_nex = pedge_cur.partial_edge_next().unwrap();
        let (vert1, _) = pedge_cur
            .partial_node_first()
            .unwrap()
            .node()
            .center_and_radius()?;
        let (vert2, _) = pedge_cur
            .partial_node_last()
            .unwrap()
            .node()
            .center_and_radius()?;
        let (vert3, _) = pedge_aft
            .partial_node_last()
            .unwrap()
            .node()
            .center_and_radius()?;
        let (vert4, _) = pedge_nex
            .partial_node_last()
            .unwrap()
            .node()
            .center_and_radius()?;

        let vec_cur = (vert2 - vert1).normalize();
        let vec_nex = vert4 - vert2;
        let nor = vec_cur.cross(&vec_nex).normalize();
        let vec_aft = (vert3 - vert2).normalize();
        let vec_cross = vec_cur.cross(&vec_aft);
        if vec_cross.dot(&nor) > 0.0 {
            let sin_ang = vec_cross.norm();
            if vec_cur.dot(&vec_aft) < 0.0 {
                return Ok(Some(sin_ang));
            }
        }
    }
    Ok(None)
}

pub fn excluded_alveolae(
    ind_pedge: usize,
    skeleton_interface: &SkeletonInterface3D,
) -> HashSet<usize> {
    let pnode = skeleton_interface
        .get_partial_edge_uncheck(ind_pedge)
        .partial_node_last()
        .unwrap();
    pnode
        .partial_edge_next()
        .iter()
        .filter_map(|pedge| {
            let alve = pedge.partial_alveola().alveola();
            if alve.is_full() {
                Some(alve.ind())
            } else {
                None
            }
        })
        .collect()
}

pub fn singular_path_to_exclude_alveolae(
    set_alve_to_exclude: &HashSet<usize>,
    skeleton_interface: &mut SkeletonInterface3D,
) -> Result<Option<SkeletonSingularPath>> {
    let mut opt_ind_first_pedge = None;
    for &ind_alve in set_alve_to_exclude {
        let palve = skeleton_interface
            .get_alveola_uncheck(ind_alve)
            .partial_alveolae()[0];
        for pedge in palve.partial_edges().iter() {
            if pedge.edge().degree() == 1 {
                continue;
            }
            let neigh_alve: HashSet<usize> = pedge
                .edge()
                .alveolae()
                .iter()
                .filter_map(|alve| {
                    if set_alve_to_exclude.contains(&alve.ind()) {
                        Some(alve.ind())
                    } else {
                        None
                    }
                })
                .collect();
            if neigh_alve.len() == 1 {
                opt_ind_first_pedge = Some(pedge.ind());
                break;
            }
        }
        if opt_ind_first_pedge.is_some() {
            break;
        }
    }

    if let Some(ind_first_pedge) = opt_ind_first_pedge {
        let mut sing_path = SkeletonSingularPath::create(ind_first_pedge);
        sing_path.append_last(skeleton_interface)?;

        loop {
            if let Some(ind_pedge_last) = sing_path.ind_last_partial_edge() {
                let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge_last);
                let should_append = if pedge.edge().degree() == 1 {
                    false
                } else {
                    let neigh_alve: HashSet<usize> = pedge
                        .edge()
                        .alveolae()
                        .iter()
                        .filter_map(|alve| {
                            if set_alve_to_exclude.contains(&alve.ind()) {
                                Some(alve.ind())
                            } else {
                                None
                            }
                        })
                        .collect();
                    if neigh_alve.len() == 1 {
                        true
                    } else {
                        false
                    }
                };

                if should_append {
                    sing_path.append_last(skeleton_interface)?;
                } else {
                    sing_path.rotate_last(skeleton_interface)?;
                }
            } else {
                break;
            }
        }

        Ok(Some(sing_path))
    } else {
        Ok(None)
    }
}
