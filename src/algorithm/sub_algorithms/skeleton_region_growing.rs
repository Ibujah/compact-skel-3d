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
        let (mut to_add_near, mut ind_region_near) =
            neighbors_to_add(skeleton_interface, passed_alveolae, ind_alveola)?;
        // while let Some((ind_reg, _)) = ind_region_near
        //     .into_iter()
        //     .filter(|&ind_reg_near| {
        //         can_merge_region(
        //             skeleton_interface,
        //             passed_alveolae,
        //             ind_region,
        //             ind_reg_near,
        //         )
        //     })
        //     .map(|ind_reg_near| {
        //         (
        //             ind_reg_near,
        //             score_fusion(
        //                 skeleton_interface,
        //                 passed_alveolae,
        //                 ind_alveola,
        //                 ind_reg_near,
        //             )
        //             .unwrap(),
        //         )
        //     })
        //     .fold(None, |curr_min, (ind_reg_near, score)| {
        //         if let Some((_, score_min)) = curr_min {
        //             if score_min > score {
        //                 Some((ind_reg_near, score))
        //             } else {
        //                 curr_min
        //             }
        //         } else {
        //             Some((ind_reg_near, score))
        //         }
        //     })
        // {
        //     for (_, ind_region_i) in passed_alveolae.iter_mut() {
        //         if *ind_region_i == ind_reg {
        //             *ind_region_i = ind_region;
        //         }
        //     }
        //     (to_add_near, ind_region_near) =
        //         neighbors_to_add(skeleton_interface, passed_alveolae, ind_alveola)?;
        // }
        for &ind_to_add in to_add_near.iter() {
            let score = score_alveola(skeleton_interface, passed_alveolae, ind_to_add, ind_region)?;
            near_alveolae.push((ind_to_add, ind_region, score));
        }
    }
    Ok(())
}

pub fn region_merge(
    skeleton_interface: &SkeletonInterface3D,
    passed_alveolae: &mut HashMap<usize, usize>,
) -> Result<()> {
    let mut neighboring_score: HashMap<(usize, usize), (f64, usize)> = HashMap::new();

    init_neighboring_score(
        skeleton_interface,
        passed_alveolae,
        &mut neighboring_score,
        None,
    )?;

    // get mimimum score
    while let Some((&(ind_region1, ind_region2), _)) =
        neighboring_score
            .iter()
            .fold(None, |curr_min, (ind, &(score_sum_tst, nb_tst))| {
                // If a current minimum is found, check if the current score is greater than the current minimum
                let score_tst = score_sum_tst; // / nb_tst as f64;
                if let Some((_, score_curr)) = curr_min {
                    if score_curr < score_tst {
                        Some((ind, score_tst))
                    } else {
                        curr_min
                    }
                // If no current minimum is found, set the index and score as the current minimum
                } else {
                    Some((ind, score_tst))
                }
            })
    {
        // let (score_sum, nb) =
        neighboring_score
            .remove(&(ind_region1, ind_region2))
            .unwrap();
        // let score = score_sum / nb as f64;
        // if score < 0.8 {
        //     break;
        // }
        if !can_merge_region(
            skeleton_interface,
            passed_alveolae,
            ind_region1,
            ind_region2,
        ) {
            continue;
        }
        // merge regions
        passed_alveolae.iter_mut().for_each(|(_, ind_reg)| {
            if *ind_reg == ind_region2 {
                *ind_reg = ind_region1
            }
        });

        let (no_region2, with_region2): (
            HashMap<(usize, usize), (f64, usize)>,
            HashMap<(usize, usize), (f64, usize)>,
        ) = neighboring_score
            .into_iter()
            .partition(|((ind_r1, ind_r2), _)| *ind_r1 != ind_region2 && *ind_r2 != ind_region2);

        neighboring_score = no_region2;

        let mut with_region2: HashMap<usize, (f64, usize)> = with_region2
            .iter()
            .map(|(&(ind_r1, ind_r2), &(sc, nb))| {
                if ind_r1 == ind_region2 {
                    (ind_r2, (sc, nb))
                } else {
                    (ind_r1, (sc, nb))
                }
            })
            .collect();

        // update region1 existing neighborhood
        for (&(ind_r1, ind_r2), (sc, nb)) in neighboring_score.iter_mut() {
            if ind_r1 == ind_region1 || ind_r2 == ind_region1 {
                let ind_reg_near = if ind_r1 == ind_region1 {
                    ind_r2
                } else {
                    ind_r1
                };
                if let Some((sc_up, nb_up)) = with_region2.remove(&ind_reg_near) {
                    *sc += sc_up;
                    *nb += nb_up;
                }
            }
        }

        // include new region1 neighborhood
        for (&ind_reg_near, &(sc, nb)) in with_region2.iter() {
            if ind_reg_near < ind_region1 {
                neighboring_score.insert((ind_reg_near, ind_region1), (sc, nb));
            } else {
                neighboring_score.insert((ind_region1, ind_reg_near), (sc, nb));
            }
        }
    }
    Ok(())
}

pub fn init_neighboring_score(
    skeleton_interface: &SkeletonInterface3D,
    passed_alveolae: &HashMap<usize, usize>,
    neighboring_score: &mut HashMap<(usize, usize), (f64, usize)>,
    only_region: Option<usize>,
) -> Result<()> {
    for (&ind_alveola, &ind_region) in passed_alveolae.iter() {
        let alveola = skeleton_interface.get_alveola_uncheck(ind_alveola);
        let palveola = alveola.partial_alveolae()[0];
        for pedge in palveola.partial_edges() {
            if !pedge.edge().is_regular() {
                continue;
            }

            let mut pedge_opp = pedge.partial_edge_neighbor().partial_edge_opposite();
            while !pedge_opp.partial_alveola().alveola().is_full() {
                pedge_opp = pedge_opp.partial_edge_neighbor().partial_edge_opposite();
            }
            let alveola_near = pedge_opp.partial_alveola().alveola();
            let ind_alveola_near = alveola_near.ind();
            let ind_region_near = passed_alveolae[&ind_alveola_near];
            if let Some(ind_only) = only_region {
                if ind_region != ind_only && ind_region_near != ind_only {
                    continue;
                }
            }
            if ind_region > ind_region_near {
                continue;
            }

            // let seg = palveola.alveola().delaunay_segment();
            // let v1 = skeleton_interface.get_mesh().vertices()[&seg[0]];
            // let v2 = skeleton_interface.get_mesh().vertices()[&seg[1]];
            // let normal = (v2 - v1).normalize();

            // let seg_near = palveola.alveola().delaunay_segment();
            // let v1_near = skeleton_interface.get_mesh().vertices()[&seg_near[0]];
            // let v2_near = skeleton_interface.get_mesh().vertices()[&seg_near[1]];
            // let normal_near = (v2_near - v1_near).normalize();

            // let cos_ang = normal.dot(&normal_near).abs();

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

            let length = (v1 - v2).norm();

            neighboring_score
                .entry((ind_region, ind_region_near))
                .and_modify(|(v, nb)| {
                    *v += length;
                    *nb += 1
                })
                .or_insert((length, 1));
        }
    }
    Ok(())
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
) -> Result<(Vec<usize>, Vec<usize>)> {
    let alveola = skeleton_interface.get_alveola(ind_alveola)?;
    let &ind_region_curr = passed_alveolae.get(&ind_alveola).unwrap();
    let mut to_add = Vec::new();
    let mut in_region = Vec::new();
    for edge in alveola.edges().iter() {
        if !edge.is_regular() {
            continue;
        }
        for alveola in edge.alveolae().iter() {
            if !alveola.is_full() {
                continue;
            }
            if let Some(&ind_region) = passed_alveolae.get(&alveola.ind()) {
                if ind_region != ind_region_curr {
                    in_region.push(ind_region);
                }
            } else {
                to_add.push(alveola.ind());
            }
        }
    }
    Ok((to_add, in_region))
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

        score += length;
    }
    Ok(score)
}

pub fn score_fusion(
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

        let length = (v1 - v2).norm();

        let mut pedge_opp = pedge.partial_edge_neighbor();
        while !pedge_opp.partial_alveola().alveola().is_full() {
            pedge_opp = pedge_opp.partial_edge_opposite().partial_edge_neighbor();
        }

        let alveola_neigh = pedge_opp.partial_alveola().alveola();

        if let Some(&ind_reg) = passed_alveolae.get(&alveola_neigh.ind()) {
            if ind_region == ind_reg {
                score -= length;
            }
        }
    }
    Ok(score)
}

pub fn alveola_regular_perimeter(
    skeleton_interface: &SkeletonInterface3D,
    ind_alveola: usize,
) -> Result<f64> {
    let alveola = skeleton_interface.get_alveola(ind_alveola)?;

    let mut regular_perimeter = 0.0;

    for pedge in alveola.partial_alveolae()[0].partial_edges() {
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

        let length = (v1 - v2).norm();

        regular_perimeter += length;
    }
    Ok(regular_perimeter)
}

pub fn can_merge_region(
    skeleton_interface: &SkeletonInterface3D,
    passed_alveolae: &HashMap<usize, usize>,
    ind_region1: usize,
    ind_region2: usize,
) -> bool {
    !passed_alveolae
        .iter()
        .filter_map(|(&ind_alveola, &ind_region)| {
            // alveolae from first region
            if ind_region == ind_region1 {
                Some(ind_alveola)
            } else {
                None
            }
        })
        .flat_map(|ind_alveola| {
            // edges
            skeleton_interface.get_alveola_uncheck(ind_alveola).edges()
        })
        .any(|edge| {
            if edge.is_regular() {
                false
            } else {
                let mut cpt = 0;
                for alve in edge.alveolae() {
                    if let Some(&ind_reg) = passed_alveolae.get(&alve.ind()) {
                        if ind_reg == ind_region1 || ind_reg == ind_region2 {
                            cpt += 1;
                        }
                    }
                }
                cpt == 3
            }
        })
}
