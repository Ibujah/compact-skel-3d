use anyhow::Result;

use crate::mesh3d::GenericMesh3D;
use crate::mesh3d::ManifoldMesh3D;
use crate::skeleton3d::Skeleton3D;

use super::sub_algorithms::skeleton_operations;
use super::sub_algorithms::SkeletonInterface3D;

/// Computes the full skeletonization of a delaunay mesh
pub fn full_skeletonization(mesh: &mut ManifoldMesh3D) -> Result<Skeleton3D> {
    println!("Init skeleton interface");
    let mut skeleton_interface = SkeletonInterface3D::init(mesh)?;

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
    opt_epsilon: Option<f32>,
) -> Result<()> {
    println!("Finding some first alveola");
    let mut ind_first_alveola = skeleton_operations::first_alveola_in(skeleton_interface)?;
    let mut cpt_loop = 0;
    let mut nb_sheets_prev = 0;
    loop {
        cpt_loop = cpt_loop + 1;
        let mut label = 1;
        let mut modif_done = false;

        skeleton_interface.reinit_skeleton();
        println!("Loop {}", cpt_loop);
        println!("Propagating first sheet");
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
                            ind_first_alveola = ind_alveola;
                        }

                        let mut vec_pedges_boundary = skeleton_operations::boundary_partial_edges(
                            &skeleton_interface,
                            &current_sheet,
                        );

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
    }
    println!("Checking skeleton");
    skeleton_interface.check()?;
    Ok(())
}

/// Computes the sheet based skeletonization of a delaunay mesh
pub fn sheet_skeletonization(
    mesh: &mut ManifoldMesh3D,
    opt_epsilon: Option<f32>,
) -> Result<(Skeleton3D, Vec<GenericMesh3D>)> {
    println!("Init skeleton interface");
    let mut skeleton_interface = SkeletonInterface3D::init(mesh)?;
    skeleton_interface.check()?;

    if let Some(err) = loop_skeletonization(&mut skeleton_interface, opt_epsilon).err() {
        println!("{}", err);
    }

    Ok((
        skeleton_interface.get_skeleton().clone(),
        skeleton_interface.get_debug_meshes().clone(),
    ))
}
