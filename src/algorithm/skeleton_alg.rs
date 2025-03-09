use anyhow::Result;
use std::collections::HashMap;
use std::time::Instant;

use crate::algorithm::delaunay_alg;
use crate::algorithm::sub_algorithms::skeleton_operations::include_edge_in_skel;
use crate::algorithm::sub_algorithms::skeleton_operations::remap_sheet_indices;
use crate::algorithm::sub_algorithms::SkeletonSeparation;
use crate::mesh3d::GenericMesh3D;
use crate::mesh3d::ManifoldMesh3D;
use crate::skeleton3d::Skeleton3D;

use super::sub_algorithms::skeleton_operations;
use super::sub_algorithms::SkeletonInterface3D;

/// Computes the full skeletonization of a delaunay mesh
pub fn full_skeletonization(mesh: &mut ManifoldMesh3D) -> Result<Skeleton3D> {
    println!("Mesh to delaunay");
    let (faces, tetras_in) =
        delaunay_alg::to_delaunay(mesh, Some(std::f64::consts::PI * 60.0 / 180.0))?;
    println!();

    println!("Init skeleton interface");
    let mut skeleton_interface = SkeletonInterface3D::init(mesh, faces, tetras_in);

    println!("Finding some first alveola");
    let (ind_first_alveola, _) = skeleton_operations::first_alveola_in(&mut skeleton_interface)?;
    let mut vec_alveola = Vec::new();
    vec_alveola.push(ind_first_alveola.unwrap());

    println!("Propagating skeleton");
    while let Some(ind_alveola) = vec_alveola.pop() {
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
    }
    println!();

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

    let epsilon = opt_epsilon.unwrap_or(0.);
    loop {
        cpt_loop += 1;
        label = 1;
        let mut modif_done = false;
        let mut sheet_siz_max = 0;

        let mut vec_pedges = Vec::new();
        let mut vec_lone_edges = Vec::new();

        skeleton_interface.reinit_skeleton();
        println!("Loop {}", cpt_loop);
        println!("Propagating first sheet");
        let (ind_first_alveola_opt, ind_first_edge_opt) =
            skeleton_operations::first_alveola_in(skeleton_interface)?;
        if let Some(ind_first_alveola) = ind_first_alveola_opt {
            println!("First alveola : {}", ind_first_alveola);
            skeleton_operations::compute_sheet(skeleton_interface, ind_first_alveola, label)?;
            let current_sheet = skeleton_interface.get_sheet(label);
            println!("{}", current_sheet.len());
            sheet_siz_max = current_sheet.len();
            for &ind_alveola in current_sheet.iter() {
                if skeleton_interface.get_alveola(ind_alveola)?.is_full() {
                    skeleton_operations::include_alveola_in_skel(
                        skeleton_interface,
                        ind_alveola,
                        Some(label),
                    )?;
                }
            }
            let mut vec_pedges_new =
                skeleton_operations::outer_partial_edges(skeleton_interface, &current_sheet);
            vec_pedges.append(&mut vec_pedges_new);
            let mut vec_lone_edges_new =
                skeleton_operations::lone_edges(skeleton_interface, &current_sheet);
            vec_lone_edges.append(&mut vec_lone_edges_new);
        } else if let Some(ind_first_edge) = ind_first_edge_opt {
            println!("First edge : {}", ind_first_edge);
            vec_lone_edges.push(ind_first_edge)
        }
        vec_pedges.sort();
        vec_pedges.dedup();

        vec_lone_edges.sort();
        vec_lone_edges.dedup();

        println!("Searching paths");
        loop {
            if let Some(ind_pedge) = vec_pedges.pop() {
                if skeleton_interface
                    .get_partial_edge(ind_pedge)?
                    .partial_alveola()
                    .alveola()
                    .label()
                    .is_some()
                {
                    continue;
                }
                if skeleton_interface
                    .get_partial_edge(ind_pedge)?
                    .edge()
                    .degree()
                    == 1
                {
                    label += 1;
                    let ind_alveola = skeleton_interface
                        .get_partial_edge(ind_pedge)?
                        .partial_alveola()
                        .alveola()
                        .ind();
                    skeleton_operations::compute_sheet(skeleton_interface, ind_alveola, label)?;
                    let current_sheet = skeleton_interface.get_sheet(label);

                    if current_sheet.len() > sheet_siz_max {
                        sheet_siz_max = current_sheet.len();
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
                        skeleton_interface,
                        &current_sheet,
                    );
                    let mut vec_lone_edges_new =
                        skeleton_operations::lone_edges(skeleton_interface, &current_sheet);
                    vec_pedges.append(&mut vec_pedges_new);
                    vec_pedges.sort();
                    vec_pedges.dedup();
                    vec_lone_edges.append(&mut vec_lone_edges_new);
                    vec_lone_edges.sort();
                    vec_lone_edges.dedup();
                    continue;
                }
                if let Some(skeleton_separation) =
                    skeleton_operations::extract_skeleton_separation(skeleton_interface, ind_pedge)?
                {
                    let mut removed = false;
                    if skeleton_separation.closable_path()? {
                        if let Some(mesh_faces) = skeleton_operations::collect_mesh_faces_index(
                            &skeleton_separation,
                            epsilon,
                        )? {
                            if let Some(closing_faces) = skeleton_operations::collect_closing_faces(
                                &skeleton_separation,
                                &mesh_faces,
                            )? {
                                if !mesh_faces.is_empty()
                                    && !closing_faces.is_empty()
                                    && skeleton_operations::try_remove_and_add(
                                        skeleton_interface,
                                        &mesh_faces,
                                        &closing_faces,
                                    )?
                                {
                                    removed = true;
                                    modif_done = true;
                                }
                            }
                        }
                    }
                    if !removed {
                        label += 1;
                        let ind_alveola = skeleton_interface
                            .get_partial_edge(ind_pedge)?
                            .partial_alveola()
                            .alveola()
                            .ind();
                        skeleton_operations::compute_sheet(skeleton_interface, ind_alveola, label)?;
                        let current_sheet = skeleton_interface.get_sheet(label);

                        if current_sheet.len() > sheet_siz_max {
                            sheet_siz_max = current_sheet.len();
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
                            skeleton_interface,
                            &current_sheet,
                        );
                        let mut vec_lone_edges_new =
                            skeleton_operations::lone_edges(skeleton_interface, &current_sheet);
                        vec_pedges.append(&mut vec_pedges_new);
                        vec_pedges.sort();
                        vec_pedges.dedup();
                        vec_lone_edges.append(&mut vec_lone_edges_new);
                        vec_lone_edges.sort();
                        vec_lone_edges.dedup();
                    }
                }
            } else if let Some(ind_edge) = vec_lone_edges.pop() {
                let edge = skeleton_interface.get_edge(ind_edge)?;
                if !edge.is_full() {
                    continue;
                }
                if edge.degree() != 0 {
                    continue;
                }
                if edge.label().is_some() {
                    continue;
                }
                skeleton_interface.propagate_edge(ind_edge)?;
                label += 1;
                include_edge_in_skel(skeleton_interface, ind_edge)?;
                skeleton_interface.set_edge_label(ind_edge, Some(label))?;
                for node in skeleton_interface.get_edge(ind_edge)?.nodes() {
                    for edge in node.edges() {
                        if !edge.is_full() {
                            continue;
                        }
                        if edge.label().is_some() {
                            continue;
                        }
                        if edge.degree() == 0 {
                            vec_lone_edges.push(edge.ind());
                        } else if edge.degree() >= 1 {
                            for pedge in edge.partial_edges() {
                                if pedge.partial_alveola().alveola().is_full() {
                                    vec_pedges.push(pedge.ind());
                                    break;
                                }
                            }
                        }
                    }
                }
                vec_pedges.sort();
                vec_pedges.dedup();
                vec_lone_edges.sort();
                vec_lone_edges.dedup();
            } else {
                break;
            }
        }

        println!(
            "\r{} Sheets,  {} + {} pedges remaining                                   ",
            label,
            vec_pedges.len(),
            vec_lone_edges.len(),
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
                if !pedge.edge().is_full() {
                    continue;
                }
                let palve = pedge.partial_alveola();
                if palve.alveola().label().is_none() {
                    continue;
                }
                if !palve.alveola().is_full() {
                    continue;
                }

                let mut contains_nod_junction = false;
                for pedg in palve
                    .partial_edges()
                    .iter()
                    .filter(|pe| pe.is_boundary() && pe.partial_edge_next().unwrap().is_boundary())
                {
                    if !pedg.partial_edge_next().unwrap().is_boundary() {
                        continue;
                    }

                    let ind_e1 = pedg.edge().ind();
                    let ind_e2 = pedg.partial_edge_next().unwrap().edge().ind();

                    let node_last = pedg.partial_node_last().unwrap().node();
                    for edg in node_last.edges() {
                        if edg.is_boundary() && edg.ind() != ind_e1 && edg.ind() != ind_e2 {
                            for alv in edg.alveolae() {
                                if alv.is_full() {
                                    contains_nod_junction = true;
                                }
                            }
                        }
                    }
                }
                if contains_nod_junction {
                    continue;
                }

                if let Some((sing_path, vec_new_pedges, set_alve)) =
                    skeleton_operations::exclusion_singular_path(ind_pedge, skeleton_interface)?
                {
                    let skeleton_separation =
                        SkeletonSeparation::from_singular_path(skeleton_interface, sing_path);
                    if let Some(mesh_faces) = skeleton_operations::collect_mesh_faces_index(
                        &skeleton_separation,
                        epsilon,
                    )? {
                        if let Some(closing_faces) = skeleton_operations::collect_closing_faces(
                            &skeleton_separation,
                            &mesh_faces,
                        )? {
                            if !mesh_faces.is_empty()
                                && !closing_faces.is_empty()
                                && skeleton_operations::try_remove_and_add(
                                    skeleton_interface,
                                    &mesh_faces,
                                    &closing_faces,
                                )?
                            {
                                for &ind_alve in set_alve.iter() {
                                    if !skeleton_interface.get_alveola(ind_alve)?.is_full() {
                                        skeleton_interface.set_alveola_label(ind_alve, None)?;
                                    }
                                }
                                let mut new_saliencies = skeleton_operations::estimate_saliencies(
                                    skeleton_interface,
                                    &vec_new_pedges,
                                )?;
                                saliencies.append(&mut new_saliencies);
                                skeleton_operations::sort_saliencies(&mut saliencies);
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
    skeleton_operations::relabel_all_skeleton(skeleton_interface)?;
    let nb_sheets = remap_sheet_indices(skeleton_interface);
    println!("{} Sheets", nb_sheets);
    let problematics = skeleton_operations::problematic_partial_edges(skeleton_interface);
    println!("{} problematic pedges", problematics.len());
    println!("Checking skeleton");
    skeleton_interface.check()?;
    // if !skeleton_interface.check_cocone() {
    //     println!("Invalid for cocone criterion");
    // }
    Ok(())
}

/// Computes the sheet based skeletonization of a delaunay mesh
pub fn sheet_skeletonization(
    mesh: &mut ManifoldMesh3D,
    opt_epsilon: Option<f64>,
) -> Result<(
    Skeleton3D,
    ManifoldMesh3D,
    Vec<GenericMesh3D>,
    Vec<usize>,
    u64,
    u64,
)> {
    let mut mesh_cl = mesh.clone();

    println!("Mesh to delaunay");
    let now = Instant::now();
    let (faces, tetras_in) =
        delaunay_alg::to_delaunay(&mut mesh_cl, Some(std::f64::consts::PI * 60.0 / 180.0))?;
    let duration = now.elapsed();
    let del_sec = duration.as_secs();
    println!();

    println!("Init skeleton interface");
    let now = Instant::now();
    let mut skeleton_interface = SkeletonInterface3D::init(&mut mesh_cl, faces, tetras_in);
    skeleton_interface.check()?;

    if let Some(err) = loop_skeletonization(&mut skeleton_interface, opt_epsilon).err() {
        println!("{}", err);
    }
    let problematic_edges = skeleton_operations::problematic_edges(&skeleton_interface);

    println!("Computing labels");
    let label_per_vertex = skeleton_interface.get_label_per_vertex()?;
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
    let duration = now.elapsed();
    let skel_sec = duration.as_secs();

    Ok((
        skeleton_interface.get_skeleton().clone(),
        skeleton_interface.get_mesh().clone(),
        skeleton_interface.get_debug_meshes().clone(),
        problematic_edges,
        del_sec,
        skel_sec,
    ))
}
