use anyhow::Result;
use std::collections::HashMap;

use crate::algorithm::delaunay_alg;
use crate::algorithm::sub_algorithms::SkeletonSeparation;
use crate::mesh3d::GenericMesh3D;
use crate::mesh3d::ManifoldMesh3D;
use crate::skeleton3d::Skeleton3D;

use super::sub_algorithms::skeleton_operations;
use super::sub_algorithms::SkeletonInterface3D;

/// Computes the full skeletonization of a delaunay mesh
pub fn full_skeletonization(mesh: &mut ManifoldMesh3D) -> Result<Skeleton3D> {
    println!("Mesh to delaunay");
    let faces = delaunay_alg::to_delaunay(mesh, Some(std::f64::consts::PI * 20.0 / 180.0))?;
    println!("");

    println!("Init skeleton interface");
    let mut skeleton_interface = SkeletonInterface3D::init(mesh, faces);

    println!("Finding some first alveola");
    let ind_first_alveola = skeleton_operations::first_alveola_in(&mut skeleton_interface)?;
    let mut vec_alveola = Vec::new();
    vec_alveola.push(ind_first_alveola);

    println!("Propagating skeleton");
    loop {
        if let Some(ind_alveola) = vec_alveola.pop() {
            let alveola = skeleton_interface.get_alveola(ind_alveola)?;
            let alveola_in = alveola.is_full();
            if !alveola.is_computed() && alveola_in {
                skeleton_interface.compute_alveola(ind_alveola)?;
                let mut vec_neigh =
                    skeleton_operations::neighbor_alveolae(&mut skeleton_interface, ind_alveola)?;
                vec_alveola.append(&mut vec_neigh);
            }
            if alveola_in {
                skeleton_operations::include_alveola_in_skel(
                    &mut skeleton_interface,
                    ind_alveola,
                    None,
                )?;
            }
            print!("\r{} alveolae remaining     ", vec_alveola.len());
        } else {
            break;
        }
    }
    println!("");

    println!("Checking skeleton");
    skeleton_interface.check()?;

    Ok(skeleton_interface.get_skeleton().clone())
}

fn loop_skeletonization(
    skeleton_interface: &mut SkeletonInterface3D,
    opt_epsilon: Option<f64>,
) -> Result<()> {
    // println!("Finding some first alveola");
    // let mut ind_first_alveola = skeleton_operations::first_alveola_in(skeleton_interface)?;
    let mut cpt_loop = 0;
    let mut nb_sheets_prev = 0;
    let mut label;
    loop {
        cpt_loop = cpt_loop + 1;
        label = 1;
        let mut modif_done = false;

        skeleton_interface.reinit_skeleton();
        println!("Loop {}", cpt_loop);
        println!("Propagating first sheet");
        let ind_first_alveola = skeleton_operations::first_alveola_in(skeleton_interface)?;
        skeleton_operations::compute_sheet(skeleton_interface, ind_first_alveola, label)?;
        let current_sheet = skeleton_interface.get_sheet(label);
        let mut sheet_siz_max = current_sheet.len();
        for &ind_alveola in current_sheet.iter() {
            if skeleton_interface.get_alveola(ind_alveola)?.is_full() {
                skeleton_operations::include_alveola_in_skel(
                    skeleton_interface,
                    ind_alveola,
                    Some(label),
                )?;
            }
        }
        let mut vec_pedges =
            skeleton_operations::outer_partial_edges(&skeleton_interface, &current_sheet);
        vec_pedges.sort();
        vec_pedges.dedup();

        println!("Searching paths");
        loop {
            if let Some(ind_pedge) = vec_pedges.pop() {
                print!(
                    "\rSheet {},  {} pedges remaining                                   ",
                    label,
                    vec_pedges.len()
                );
                if skeleton_interface
                    .get_partial_edge(ind_pedge)?
                    .partial_alveola()
                    .alveola()
                    .label()
                    .is_some()
                {
                    continue;
                }
                if let Some(skeleton_separation) =
                    skeleton_operations::extract_skeleton_separation(skeleton_interface, ind_pedge)?
                {
                    let mut removed = false;
                    if let Some(epsilon) = opt_epsilon {
                        if skeleton_separation.closable_path()? {
                            if let Some(mesh_faces) = skeleton_operations::collect_mesh_faces_index(
                                &skeleton_separation,
                                epsilon,
                            )? {
                                if let Some(closing_faces) =
                                    skeleton_operations::collect_closing_faces(
                                        &skeleton_separation,
                                        &mesh_faces,
                                    )?
                                {
                                    if !mesh_faces.is_empty() && !closing_faces.is_empty() {
                                        if skeleton_operations::try_remove_and_add(
                                            skeleton_interface,
                                            &mesh_faces,
                                            &closing_faces,
                                        )? {
                                            removed = true;
                                            modif_done = true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if !removed {
                        label = label + 1;
                        let ind_alveola = skeleton_interface
                            .get_partial_edge(ind_pedge)?
                            .partial_alveola()
                            .alveola()
                            .ind();
                        skeleton_operations::compute_sheet(skeleton_interface, ind_alveola, label)?;
                        let current_sheet = skeleton_interface.get_sheet(label);

                        if current_sheet.len() > sheet_siz_max {
                            sheet_siz_max = current_sheet.len();
                            //ind_first_alveola = ind_alveola;
                        }

                        for &ind_alveola in current_sheet.iter() {
                            if skeleton_interface.get_alveola(ind_alveola)?.is_full() {
                                skeleton_operations::include_alveola_in_skel(
                                    skeleton_interface,
                                    ind_alveola,
                                    Some(label),
                                )?;
                            }
                        }
                        let mut vec_pedges_new = skeleton_operations::outer_partial_edges(
                            &skeleton_interface,
                            &current_sheet,
                        );
                        vec_pedges.append(&mut vec_pedges_new);
                        vec_pedges.sort();
                        vec_pedges.dedup();
                    }
                }
            } else {
                break;
            }
        }

        println!(
            "\r{} Sheets,  {} pedges remaining                                   ",
            label,
            vec_pedges.len()
        );

        if !modif_done || nb_sheets_prev == label {
            break;
        }
        nb_sheets_prev = label;

        println!("Boundary edges correction");
        let vec_pedges = skeleton_operations::boundary_partial_edges(skeleton_interface);
        let mut saliencies =
            skeleton_operations::estimate_saliencies(skeleton_interface, &vec_pedges)?;
        skeleton_operations::sort_saliencies(&mut saliencies);
        loop {
            print!(
                "\r{} boundary pedges remaining                                   ",
                saliencies.len()
            );
            if let Some((ind_pedge, _)) = saliencies.pop() {
                let pedge = skeleton_interface.get_partial_edge(ind_pedge)?;
                if pedge.edge().degree() != 1 {
                    continue;
                }
                if pedge.partial_alveola().alveola().label().is_none() {
                    continue;
                }
                if let Some((sing_path, vec_new_pedges, set_alve)) =
                    skeleton_operations::exclusion_singular_path(ind_pedge, skeleton_interface)?
                {
                    let skeleton_separation =
                        SkeletonSeparation::from_singular_path(skeleton_interface, sing_path);
                    if let Some(epsilon) = opt_epsilon {
                        if let Some(mesh_faces) = skeleton_operations::collect_mesh_faces_index(
                            &skeleton_separation,
                            epsilon,
                        )? {
                            if let Some(closing_faces) = skeleton_operations::collect_closing_faces(
                                &skeleton_separation,
                                &mesh_faces,
                            )? {
                                if !mesh_faces.is_empty() && !closing_faces.is_empty() {
                                    if skeleton_operations::try_remove_and_add(
                                        skeleton_interface,
                                        &mesh_faces,
                                        &closing_faces,
                                    )? {
                                        for &ind_alve in set_alve.iter() {
                                            if !skeleton_interface.get_alveola(ind_alve)?.is_full()
                                            {
                                                skeleton_interface
                                                    .set_alveola_label(ind_alve, None)?;
                                            }
                                        }
                                        let mut new_saliencies =
                                            skeleton_operations::estimate_saliencies(
                                                skeleton_interface,
                                                &vec_new_pedges,
                                            )?;
                                        saliencies.append(&mut new_saliencies);
                                        skeleton_operations::sort_saliencies(&mut saliencies);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                break;
            }
        }
        println!(
            "\r{} boundary pedges remaining                                   ",
            saliencies.len()
        );
    }
    println!("Problematic edges correction");

    let problematics = skeleton_operations::problematic_partial_edges(skeleton_interface);
    println!("{} problematic pedges", problematics.len());
    label = skeleton_operations::handle_all_problematic_pedge_by_region_growing(
        &problematics,
        skeleton_interface,
        label,
    )?;
    println!("{} Sheets", label,);
    let problematics = skeleton_operations::problematic_partial_edges(skeleton_interface);
    println!("{} problematic pedges", problematics.len());
    // loop {
    //     let nb_pb = problematics.len();
    //     loop {
    //         print!(
    //             "\r{} problematic pedges remaining                                   ",
    //             problematics.len()
    //         );
    //         if let Some(ind_pedge) = problematics.pop() {
    //             let pedge = skeleton_interface.get_partial_edge(ind_pedge)?;
    //             if !pedge.edge().is_non_manifold() {
    //                 continue;
    //             }
    //             if pedge.partial_alveola().alveola().label().is_none() {
    //                 continue;
    //             }
    //             label = skeleton_operations::handle_problematic_pedge(
    //                 ind_pedge,
    //                 skeleton_interface,
    //                 label,
    //             )?;
    //         } else {
    //             break;
    //         }
    //     }
    //     problematics = skeleton_operations::problematic_partial_edges(skeleton_interface);
    //     if nb_pb == problematics.len() {
    //         break;
    //     }
    // }
    // println!(
    //     "\r{} problematic pedges remaining                                   ",
    //     problematics.len()
    // );
    println!("Checking skeleton");
    skeleton_interface.check()?;
    Ok(())
}

/// Computes the sheet based skeletonization of a delaunay mesh
pub fn sheet_skeletonization(
    mesh: &mut ManifoldMesh3D,
    opt_epsilon: Option<f64>,
) -> Result<(Skeleton3D, ManifoldMesh3D, Vec<GenericMesh3D>, Vec<usize>)> {
    let mut mesh_cl = mesh.clone();

    println!("Mesh to delaunay");
    let faces = delaunay_alg::to_delaunay(&mut mesh_cl, Some(std::f64::consts::PI * 20.0 / 180.0))?;
    println!("");

    println!("Init skeleton interface");
    let mut skeleton_interface = SkeletonInterface3D::init(&mut mesh_cl, faces);
    skeleton_interface.check()?;

    if let Some(err) = loop_skeletonization(&mut skeleton_interface, opt_epsilon).err() {
        println!("{}", err);
    }
    let problematic_edges = skeleton_operations::problematic_edges(&skeleton_interface);

    println!("Computing labels");
    let label_per_vertex = skeleton_interface.get_label_per_vertex()?;
    let mut assignment: Vec<(usize, usize)> = Vec::new();
    for (&ind_face, _) in mesh.faces() {
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
        mesh.set_face_in_group(*ind_face, lab.clone());
    }

    Ok((
        skeleton_interface.get_skeleton().clone(),
        skeleton_interface.get_mesh().clone(),
        skeleton_interface.get_debug_meshes().clone(),
        problematic_edges,
    ))
}
