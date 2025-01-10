use anyhow::Result;
use std::collections::HashMap;
use std::collections::HashSet;
use std::vec;

use crate::{mesh3d::ManifoldMesh3D, skeleton3d::Skeleton3D};

/// Returns singular edges on skeleton
pub fn seeds(skeleton: &Skeleton3D) -> Vec<usize> {
    let mut vec_seeds = Vec::new();
    for (_, vec_alve) in skeleton.get_edges_alveolae().iter() {
        if vec_alve.len() > 2 {
            for &i in vec_alve.iter() {
                vec_seeds.push(i);
            }
        }
    }

    vec_seeds
}

/// Computes permeter of an alveola
pub fn alveola_regular_perimeter(skeleton: &Skeleton3D, ind_alveola: usize) -> Result<f64> {
    let alveola_edges = skeleton.get_alveolae_edges().get(&ind_alveola).unwrap();

    let mut regular_perimeter = 0.0;

    for ind_edge in alveola_edges.iter() {
        // compute edge length
        let [ind1, ind2] = skeleton.get_edges_on_alv().get(ind_edge).unwrap();

        let v1 = skeleton.get_nodes().get(ind1).unwrap().center;
        let v2 = skeleton.get_nodes().get(ind2).unwrap().center;

        let length = (v1 - v2).norm();

        regular_perimeter += length;
    }
    Ok(regular_perimeter)
}

/// Finds next alveola to add for region growing
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

/// Finds neighbors alveolae to add to given alveola
pub fn neighbors_to_add(
    skeleton: &Skeleton3D,
    passed_alveolae: &HashMap<usize, usize>,
    ind_alveola: usize,
) -> Result<(Vec<usize>, Vec<usize>)> {
    let alveola_edges = skeleton.get_alveolae_edges().get(&ind_alveola).unwrap();
    let &ind_region_curr = passed_alveolae.get(&ind_alveola).unwrap();
    let mut to_add = Vec::new();
    let mut in_region = Vec::new();
    for ind_edge in alveola_edges.iter() {
        let edge_alveolae = skeleton.get_edges_alveolae().get(ind_edge).unwrap();
        if edge_alveolae.len() != 2 {
            continue;
        }
        for &alveola in edge_alveolae.iter() {
            if let Some(&ind_region) = passed_alveolae.get(&alveola) {
                if ind_region != ind_region_curr {
                    in_region.push(ind_region);
                }
            } else {
                to_add.push(alveola);
            }
        }
    }
    Ok((to_add, in_region))
}

/// Computes score of ading a given alveola in a given region:
/// difference between bonudary length after and before
pub fn score_alveola(
    skeleton: &Skeleton3D,
    passed_alveolae: &HashMap<usize, usize>,
    ind_alveola: usize,
    ind_region: usize,
) -> Result<f64> {
    let alveola_edges = skeleton.get_alveolae_edges().get(&ind_alveola).unwrap();

    let mut score = 0.0;

    for ind_edge in alveola_edges.iter() {
        let vec_edg_alve = skeleton.get_edges_alveolae().get(ind_edge).unwrap();
        if vec_edg_alve.len() != 2 {
            continue;
        }

        // compute edge length
        let [ind1, ind2] = skeleton.get_edges_on_alv().get(ind_edge).unwrap();

        let v1 = skeleton.get_nodes().get(ind1).unwrap().center;
        let v2 = skeleton.get_nodes().get(ind2).unwrap().center;

        let mut length = (v1 - v2).norm();

        let ind_alve_neigh = if vec_edg_alve[0] == ind_alveola {
            vec_edg_alve[1]
        } else {
            vec_edg_alve[0]
        };

        if let Some(&ind_reg) = passed_alveolae.get(&ind_alve_neigh) {
            if ind_region == ind_reg {
                length = -length;
            }
        }

        score += length;
    }
    Ok(score)
}

/// Skeletal sheet region growing function
pub fn region_grow_skel(
    skeleton: &Skeleton3D,
    passed_alveolae: &mut HashMap<usize, usize>,
    near_alveolae: &mut Vec<(usize, usize, f64)>,
) -> Result<()> {
    while let Some((ind_alveola, ind_region)) = next_to_add(passed_alveolae, near_alveolae) {
        passed_alveolae.insert(ind_alveola, ind_region);
        let (to_add_near, _) = neighbors_to_add(skeleton, passed_alveolae, ind_alveola)?;
        for &ind_to_add in to_add_near.iter() {
            let score = score_alveola(skeleton, passed_alveolae, ind_to_add, ind_region)?;
            near_alveolae.push((ind_to_add, ind_region, score));
        }
    }
    Ok(())
}
/// Skeletal sheet region growing function
pub fn search_forgot_alveolae(
    skeleton: &Skeleton3D,
    passed_alveolae: &mut HashMap<usize, usize>,
    near_alveolae: &mut Vec<(usize, usize, f64)>,
    nb_regions: &mut usize,
) -> () {
    for (ind_alv, _) in skeleton.get_alveolae().iter() {
        if !passed_alveolae.contains_key(ind_alv) {
            let perimeter = alveola_regular_perimeter(skeleton, *ind_alv).unwrap();
            near_alveolae.push((*ind_alv, *nb_regions, -perimeter));
            *nb_regions += 1;
            break;
        }
    }
}

/// Inits region growing score
pub fn init_neighboring_score(
    skeleton: &Skeleton3D,
    passed_alveolae: &HashMap<usize, usize>,
    neighboring_score: &mut HashMap<(usize, usize), (f64, usize)>,
    only_region: Option<usize>,
) -> Result<()> {
    for (&ind_alveola, &ind_region) in passed_alveolae.iter() {
        let alveola_edges = skeleton.get_alveolae_edges().get(&ind_alveola).unwrap();
        for ind_edge in alveola_edges.iter() {
            let vec_edg_alve = skeleton.get_edges_alveolae().get(ind_edge).unwrap();
            if vec_edg_alve.len() != 2 {
                continue;
            }

            let ind_alveola_near = if vec_edg_alve[0] == ind_alveola {
                vec_edg_alve[1]
            } else {
                vec_edg_alve[0]
            };

            let ind_region_near = passed_alveolae[&ind_alveola_near];
            if let Some(ind_only) = only_region {
                if ind_region != ind_only && ind_region_near != ind_only {
                    continue;
                }
            }
            if ind_region > ind_region_near {
                continue;
            }

            let [ind1, ind2] = skeleton.get_edges_on_alv().get(ind_edge).unwrap();

            let v1 = skeleton.get_nodes().get(ind1).unwrap().center;
            let v2 = skeleton.get_nodes().get(ind2).unwrap().center;

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

/// Check if two regions can be merged, i.e. if it does not create problematic edges
pub fn can_merge_region(
    skeleton: &Skeleton3D,
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
            skeleton.get_alveolae_edges().get(&ind_alveola).unwrap()
        })
        .any(|ind_edge| {
            let vec_edg_alve = skeleton.get_edges_alveolae().get(ind_edge).unwrap();
            if vec_edg_alve.len() == 2 {
                false
            } else {
                let mut cpt = 0;
                for ind_alve in vec_edg_alve.iter() {
                    if let Some(&ind_reg) = passed_alveolae.get(&ind_alve) {
                        if ind_reg == ind_region1 || ind_reg == ind_region2 {
                            cpt += 1;
                        }
                    }
                }
                cpt == 3
            }
        })
}

/// Skeletal sheet region merging function
pub fn region_merge(
    skeleton: &Skeleton3D,
    passed_alveolae: &mut HashMap<usize, usize>,
) -> Result<()> {
    let mut neighboring_score: HashMap<(usize, usize), (f64, usize)> = HashMap::new();

    init_neighboring_score(skeleton, passed_alveolae, &mut neighboring_score, None)?;

    // get mimimum score
    while let Some((&(ind_region1, ind_region2), _)) =
        neighboring_score
            .iter()
            .fold(None, |curr_min, (ind, &(score_sum_tst, _))| {
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
        if !can_merge_region(skeleton, passed_alveolae, ind_region1, ind_region2) {
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

/// associates a skeleton label for each mesh vertex
pub fn get_label_per_vertex(
    skeleton: &Skeleton3D,
    mesh: &ManifoldMesh3D,
) -> Result<HashMap<usize, HashSet<usize>>> {
    let mut labels_per_node: HashMap<usize, HashSet<usize>> = HashMap::new();
    let mut labels_per_vert: HashMap<usize, HashSet<usize>> = HashMap::new();
    let mut node_per_vert: Vec<usize> = vec![0; mesh.get_nb_vertices()];

    for ind_vertex in 0..mesh.get_nb_vertices() {
        let coord = mesh.get_vertex(ind_vertex)?.vertex();
        let mut cur_ind_opt = None;
        let mut cur_dist_opt = None;
        for (ind_nod, sph) in skeleton.get_nodes() {
            let dist = (coord - sph.center).norm() - sph.radius;
            if let Some(cur_dist) = cur_dist_opt {
                if dist < cur_dist {
                    cur_ind_opt = Some(ind_nod);
                    cur_dist_opt = Some(dist);
                }
            } else {
                cur_ind_opt = Some(ind_nod);
                cur_dist_opt = Some(dist);
            }
        }
        if let Some(&cur_ind) = cur_ind_opt {
            node_per_vert[ind_vertex] = cur_ind;
        }
    }

    for (ind_alv, nods) in skeleton.get_alveolae().iter() {
        if let Some(label_opt) = skeleton.get_labels().get(ind_alv) {
            if let Some(label) = label_opt {
                for &nod in nods.iter() {
                    labels_per_node
                        .entry(nod)
                        .or_insert(HashSet::new())
                        .insert(*label);
                }
            }
        }
    }

    for ind_vert in 0..node_per_vert.len() {
        let nod = node_per_vert[ind_vert];
        // get label
        if let Some(lab) = labels_per_node.get(&nod) {
            labels_per_vert.insert(ind_vert, lab.clone());
        } else {
            labels_per_vert.insert(ind_vert, HashSet::new());
        }
    }

    Ok(labels_per_vert)
}

/// Estimates regions on given skeleton
pub fn compute_regions(
    skeleton: &mut Skeleton3D,
    mesh_opt: &mut Option<ManifoldMesh3D>,
) -> Result<usize> {
    let mut seeds = seeds(skeleton);
    seeds.sort();
    seeds.dedup();

    let mut passed_alveolae = HashMap::new();
    let mut near_alveolae = Vec::new();
    seeds
        .iter()
        .enumerate()
        .for_each(|(ind_region, &ind_alveola)| {
            let perimeter = alveola_regular_perimeter(skeleton, ind_alveola).unwrap();
            near_alveolae.push((ind_alveola, ind_region, -perimeter));
        });

    let mut nb_regions = seeds.len();
    println!("region grow");
    while !near_alveolae.is_empty() {
        region_grow_skel(skeleton, &mut passed_alveolae, &mut near_alveolae)?;

        search_forgot_alveolae(
            skeleton,
            &mut passed_alveolae,
            &mut near_alveolae,
            &mut nb_regions,
        );
    }

    println!("region merge");
    region_merge(skeleton, &mut passed_alveolae)?;

    for (&ind_alveola, &ind_region) in passed_alveolae.iter() {
        skeleton.set_label(ind_alveola, ind_region + 1);
    }

    if let Some(mesh) = mesh_opt {
        println!("Computing labels");

        let label_per_vertex = get_label_per_vertex(&skeleton, &mesh)?;
        let mut assignment: Vec<(usize, usize)> = Vec::new();
        for ind_face in 0..mesh.get_nb_faces() {
            let vert_inds = mesh.get_face(ind_face)?.vertices_inds();
            let mut nb_vote_per_lab = HashMap::new();
            for ind_v in vert_inds.iter() {
                if let Some(list_lab) = label_per_vertex.get(ind_v) {
                    for &lab in list_lab.iter() {
                        nb_vote_per_lab
                            .entry(lab)
                            .and_modify(|c| *c += 1)
                            .or_insert(1);
                    }
                }
            }

            let (opt_lab, _) =
                nb_vote_per_lab
                    .iter()
                    .fold((None, 0), |(lab, nb), (&lab_cur, &nb_cur)| {
                        if nb_cur > nb {
                            (Some(lab_cur), nb_cur)
                        } else {
                            (lab, nb)
                        }
                    });

            if let Some(lab) = opt_lab {
                assignment.push((ind_face, lab));
            }
        }
        for (ind_face, lab) in assignment.iter() {
            mesh.set_face_in_group(*ind_face, *lab);
        }
    }
    Ok(seeds.len())
}
