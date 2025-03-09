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
    label_alveolae: &Vec<Option<usize>>,
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
            if label_alveolae[ind_alveola].is_some() {
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
    passed_alveolae: &Vec<Option<usize>>,
    ind_alveola: usize,
) -> Result<(Vec<usize>, Vec<usize>)> {
    let alveola_edges = skeleton.get_alveolae_edges().get(&ind_alveola).unwrap();
    let ind_region_curr = passed_alveolae[ind_alveola].unwrap();
    let mut to_add = Vec::new();
    let mut in_region = Vec::new();
    for ind_edge in alveola_edges.iter() {
        let edge_alveolae = skeleton.get_edges_alveolae().get(ind_edge).unwrap();
        if edge_alveolae.len() != 2 {
            continue;
        }
        for &alveola in edge_alveolae.iter() {
            if let Some(ind_region) = passed_alveolae[alveola] {
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
    passed_alveolae: &Vec<Option<usize>>,
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

        if let Some(ind_reg) = passed_alveolae[ind_alve_neigh] {
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
    label_alveolae: &mut Vec<Option<usize>>,
    near_alveolae: &mut Vec<(usize, usize, f64)>,
) -> Result<()> {
    while let Some((ind_alveola, ind_region)) = next_to_add(label_alveolae, near_alveolae) {
        label_alveolae[ind_alveola] = Some(ind_region);
        let (to_add_near, _) = neighbors_to_add(skeleton, label_alveolae, ind_alveola)?;
        for &ind_to_add in to_add_near.iter() {
            let score = score_alveola(skeleton, label_alveolae, ind_to_add, ind_region)?;
            near_alveolae.push((ind_to_add, ind_region, score));
        }
    }
    Ok(())
}
/// Skeletal sheet region growing function
pub fn search_forgot_alveolae(
    skeleton: &Skeleton3D,
    passed_alveolae: &mut Vec<Option<usize>>,
    near_alveolae: &mut Vec<(usize, usize, f64)>,
    nb_regions: &mut usize,
) -> () {
    for (&ind_alv, _) in skeleton.get_alveolae().iter() {
        if passed_alveolae[ind_alv].is_none() {
            let perimeter = alveola_regular_perimeter(skeleton, ind_alv).unwrap();
            near_alveolae.push((ind_alv, *nb_regions, -perimeter));
            *nb_regions += 1;
            break;
        }
    }
}

/// Inits region growing score
pub fn init_neighboring_score(
    skeleton: &Skeleton3D,
    passed_alveolae: &Vec<Option<usize>>,
    neighboring_score: &mut HashMap<(usize, usize), f64>,
) -> Result<(Vec<usize>, Vec<usize>, Vec<f64>)> {
    let mut i_score = Vec::new();
    let mut j_score = Vec::new();
    let mut score = Vec::new();

    for ind_alveola in 0..passed_alveolae.len() {
        let ind_region = passed_alveolae[ind_alveola].unwrap();
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

            let ind_region_near = passed_alveolae[ind_alveola_near].unwrap();
            if ind_region > ind_region_near {
                continue;
            }

            let [ind1, ind2] = skeleton.get_edges_on_alv().get(ind_edge).unwrap();

            let v1 = skeleton.get_nodes().get(ind1).unwrap().center;
            let v2 = skeleton.get_nodes().get(ind2).unwrap().center;

            let length = (v1 - v2).norm();

            neighboring_score
                .entry((ind_region, ind_region_near))
                .and_modify(|v| {
                    *v += length;
                })
                .or_insert(length);
        }
    }
    for (&(i, j), &sc) in neighboring_score.iter() {
        i_score.push(i);
        j_score.push(j);
        score.push(sc);
    }

    Ok((i_score, j_score, score))
}

/// Check if two regions can be merged, i.e. if it does not create problematic edges
pub fn can_merge_region(
    skeleton: &Skeleton3D,
    passed_alveolae: &Vec<Option<usize>>,
    ind_region1: usize,
    ind_region2: usize,
) -> bool {
    !passed_alveolae
        .iter()
        .enumerate()
        .filter_map(|(ind_alveola, &ind_region)| {
            // alveolae from first region
            if ind_region == Some(ind_region1) {
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
                for &ind_alve in vec_edg_alve.iter() {
                    if let Some(ind_reg) = passed_alveolae[ind_alve] {
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
pub fn region_merge(skeleton: &Skeleton3D, passed_alveolae: &mut Vec<Option<usize>>) -> Result<()> {
    // neighboring score is a sparse matrix, representing border distance between regions
    let mut neighboring_score: HashMap<(usize, usize), f64> = HashMap::new();

    println!("init_neighboring_score");
    let (mut i_score, mut j_score, mut score) =
        init_neighboring_score(skeleton, passed_alveolae, &mut neighboring_score)?;

    println!("loop");
    while let Some((ind_cpl_min, _)) =
        score
            .iter()
            .enumerate()
            .fold(None, |curr_min, (ind_tst, score_tst)| {
                // If a current minimum is found, check if the current score is greater than the current minimum
                if let Some((_, score_curr)) = curr_min {
                    if score_curr < score_tst {
                        Some((ind_tst, score_tst))
                    } else {
                        curr_min
                    }
                // If no current minimum is found, set the index and score as the current minimum
                } else {
                    Some((ind_tst, score_tst))
                }
            })
    {
        let ind_region1 = i_score[ind_cpl_min];
        let ind_region2 = j_score[ind_cpl_min];

        // remove entry
        i_score.remove(ind_cpl_min);
        j_score.remove(ind_cpl_min);
        score.remove(ind_cpl_min);

        // check if regions can be merged
        if !can_merge_region(skeleton, passed_alveolae, ind_region1, ind_region2) {
            continue;
        }

        // merge regions
        passed_alveolae.iter_mut().for_each(|ind_reg| {
            if *ind_reg == Some(ind_region2) {
                *ind_reg = Some(ind_region1)
            }
        });

        // update neighboring score: region 2 is removed
        // first get all neighboring scores with region 2, remove them for i_score, j_score and score

        let mut to_update = Vec::new();

        for i in 0..score.len() {
            if i_score[i] == ind_region2 {
                to_update.push((j_score[i], score[i]));
                score[i] = 0.;
            }
            if j_score[i] == ind_region2 {
                to_update.push((i_score[i], score[i]));
                score[i] = 0.;
            }
        }

        // update already existing neighborhood
        for i in 0..score.len() {
            if i_score[i] == ind_region1 {
                let reg_nei = j_score[i];
                // search reg_nei in to_update
                if let Some((ind, _)) = to_update
                    .iter()
                    .enumerate()
                    .find(|(_, &(reg, _))| reg == reg_nei)
                {
                    let (_, sc) = to_update[ind];
                    score[i] += sc;
                    to_update.remove(ind);
                }
            }
        }
        // add new neighborhood
        for i in 0..to_update.len() {
            let (reg_nei, sc) = to_update[i];
            i_score.push(ind_region1);
            j_score.push(reg_nei);
            score.push(sc);
        }

        // remove entries from i_score, j_score and score
        i_score = i_score
            .iter()
            .enumerate()
            .filter_map(|(ind, &r)| if score[ind] != 0. { Some(r) } else { None })
            .collect();
        j_score = j_score
            .iter()
            .enumerate()
            .filter_map(|(ind, &r)| if score[ind] != 0. { Some(r) } else { None })
            .collect();
        score = score
            .iter()
            .filter_map(|&sc| if sc != 0. { Some(sc) } else { None })
            .collect();
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

    let mut label_alveolae = vec![None; skeleton.get_alveolae().len()];
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
        region_grow_skel(skeleton, &mut label_alveolae, &mut near_alveolae)?;

        search_forgot_alveolae(
            skeleton,
            &mut label_alveolae,
            &mut near_alveolae,
            &mut nb_regions,
        );
    }

    println!("region merge");
    region_merge(skeleton, &mut label_alveolae)?;

    for ind_alveola in 0..label_alveolae.len() {
        let ind_region = label_alveolae[ind_alveola].unwrap();
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
