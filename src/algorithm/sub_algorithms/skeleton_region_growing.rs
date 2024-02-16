use anyhow::Result;
use std::collections::HashMap;

use super::SkeletonInterface3D;

pub fn region_grow(
    skeleton_interface: &SkeletonInterface3D,
    passed_alveolae: &mut HashMap<usize, usize>,
    near_alveolae: &mut Vec<(usize, usize, f64)>,
) -> Result<()> {
    while let Some((ind_alveola, ind_region)) = next_to_add(passed_alveolae, near_alveolae) {
        passed_alveolae.insert(ind_alveola, ind_region);
        let to_add_near = neighbors_to_add(skeleton_interface, passed_alveolae, ind_alveola)?;
        for &ind_to_add in to_add_near.iter() {
            let score = score_alveola(skeleton_interface, passed_alveolae, ind_to_add, ind_region)?;
            near_alveolae.push((ind_to_add, ind_region, score));
        }
    }
    Ok(())
}

pub fn init_near_alveolae(
    skeleton_interface: &SkeletonInterface3D,
    passed_alveolae: &HashMap<usize, usize>,
) -> Result<Vec<(usize, usize, f64)>> {
    let mut near_alveolae = Vec::new();
    for (&ind_alveola, &ind_region) in passed_alveolae.iter() {
        let to_add_near = neighbors_to_add(skeleton_interface, passed_alveolae, ind_alveola)?;
        for &ind_to_add in to_add_near.iter() {
            let score = score_alveola(skeleton_interface, passed_alveolae, ind_to_add, ind_region)?;
            near_alveolae.push((ind_to_add, ind_region, score));
        }
    }
    Ok(near_alveolae)
}

pub fn next_to_add(
    passed_alveolae: &HashMap<usize, usize>,
    near_alveolae: &mut Vec<(usize, usize, f64)>,
) -> Option<(usize, usize)> {
    // Loop until a suitable alveolus is found
    loop {
        // Get the minimum score from the `near_alveolae` vector
        if let Some((ind_min, _)) = near_alveolae
            .iter()
            .map(|(_, _, score)| score)
            .enumerate()
            .fold(None, |curr_min, (ind, score)| {
                // If a current minimum is found, check if the current score is greater than the current minimum
                if let Some((_, score_curr)) = curr_min {
                    if score_curr > score {
                        Some((ind, score))
                    } else {
                        curr_min
                    }
                // If no current minimum is found, set the index and score as the current minimum
                } else {
                    Some((ind, score))
                }
            })
        {
            // Get the alveola and region indices from the `near_alveolae` vector at the minimum index
            let (ind_alveola, ind_region, _) = near_alveolae.remove(ind_min);
            // Check if the alveola is already passed
            if passed_alveolae.contains_key(&ind_alveola) {
                continue;
            } else {
                break Some((ind_alveola, ind_region));
            }
        } else {
            // If no suitable alveola is found, return None
            break None;
        }
    }
}

pub fn neighbors_to_add(
    skeleton_interface: &SkeletonInterface3D,
    passed_alveolae: &HashMap<usize, usize>,
    ind_alveola: usize,
) -> Result<Vec<usize>> {
    let alveola = skeleton_interface.get_alveola(ind_alveola)?;
    let mut to_add = Vec::new();
    for edge in alveola.edges().iter() {
        if !edge.is_regular() {
            continue;
        }
        for alveola in edge.alveolae().iter() {
            if !alveola.is_full() {
                continue;
            }
            if passed_alveolae.contains_key(&alveola.ind()) {
                continue;
            }
            to_add.push(alveola.ind());
        }
    }
    Ok(to_add)
}

pub fn score_alveola(
    skeleton_interface: &SkeletonInterface3D,
    passed_alveolae: &HashMap<usize, usize>,
    ind_alveola: usize,
    ind_region: usize,
) -> Result<f64> {
    let alveola = skeleton_interface.get_alveola(ind_alveola)?;

    let mut score = 0.0;

    for pedge in alveola.partial_alveolae()[0].partial_edges() {
        if !pedge.edge().is_regular() {
            continue;
        }

        // compute edge length
        let v1 = pedge
            .partial_node_first()
            .unwrap()
            .node()
            .center_and_radius()?
            .0;
        let v2 = pedge
            .partial_node_last()
            .unwrap()
            .node()
            .center_and_radius()?
            .0;

        let mut length = (v1 - v2).norm();

        let mut pedge_opp = pedge.partial_edge_neighbor();
        while !pedge_opp.partial_alveola().alveola().is_full() {
            pedge_opp = pedge_opp.partial_edge_opposite().partial_edge_neighbor();
        }

        let alveola_neigh = pedge_opp.partial_alveola().alveola();

        if let Some(&ind_reg) = passed_alveolae.get(&alveola_neigh.ind()) {
            if ind_region == ind_reg {
                length = -length;
            }
        }

        score = score + length;
    }
    Ok(score)
}
