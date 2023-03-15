use crate::algorithm::skeleton_interface::{skeleton_operations, SkeletonInterface3D};
use crate::mesh3d::GenericMesh3D;
use crate::mesh3d::ManifoldMesh3D;
use crate::skeleton3d::Skeleton3D;
use anyhow::Result;

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

pub fn sheet_skeletonization(mesh: &mut ManifoldMesh3D) -> Result<(Skeleton3D, GenericMesh3D)> {
    println!("Init skeleton interface");
    let mut skeleton_interface = SkeletonInterface3D::init(mesh)?;
    skeleton_interface.check()?;

    println!("Finding some first alveola");
    let ind_first_alveola = skeleton_operations::first_alveola_in(&mut skeleton_interface)?;
    let mut vec_alveola = Vec::new();
    vec_alveola.push(ind_first_alveola);

    println!("Propagating sheet");
    let mut label = 1;
    loop {
        if let Some(ind_alveola) = vec_alveola.pop() {
            print!(
                "\rSheet {},  {} alveolae remaining     ",
                label,
                vec_alveola.len()
            );
            let alveola = skeleton_interface.get_alveola(ind_alveola)?;

            if alveola.is_full() && alveola.label().is_none() {
                skeleton_operations::compute_sheet(&mut skeleton_interface, ind_alveola, label)?;
                let current_sheet = skeleton_interface.get_sheet(label);
                let mut pedges_set =
                    skeleton_operations::outer_partial_edges(&skeleton_interface, &current_sheet);
                loop {
                    if let Some(skeleton_path) = skeleton_operations::extract_one_skeleton_path(
                        &mut skeleton_interface,
                        &pedges_set,
                    )? {
                        // if skeleton_path.closable_path()? {
                        //     if let Some(mesh_faces) = skeleton_path.collect_mesh_faces()? {
                        //         if let Some(closing_faces) =
                        //             skeleton_path.collect_closing_faces()?
                        //         {
                        //             skeleton_path.compute_debug_mesh()?;
                        //             todo!();
                        //         }
                        //     }
                        // }
                        for ind_pedge in skeleton_path.ind_partial_edges().iter() {
                            pedges_set.remove(&ind_pedge);
                            vec_alveola.push(
                                skeleton_interface
                                    .get_partial_edge(*ind_pedge)?
                                    .partial_alveola()
                                    .alveola()
                                    .ind(),
                            );
                        }
                        // pedges_set.retain(|&ind_pedge| {
                        //     skeleton_interface
                        //         .get_partial_edge(ind_pedge)
                        //         .unwrap()
                        //         .edge()
                        //         .is_singular()
                        // });
                    } else {
                        break;
                    }
                }

                for &ind_alveola in current_sheet.iter() {
                    if skeleton_interface.get_alveola(ind_alveola)?.is_full() {
                        skeleton_operations::include_alveola_in_skel(
                            &mut skeleton_interface,
                            ind_alveola,
                            Some(label),
                        )?;
                    }
                }
                label = label + 1;
            }
        } else {
            break;
        }
    }
    println!(
        "\r{} Sheets,  {} alveolae remaining     ",
        label - 1,
        vec_alveola.len()
    );

    println!("Checking skeleton");
    skeleton_interface.check()?;

    Ok((
        skeleton_interface.get_skeleton().clone(),
        skeleton_interface.get_debug_mesh().clone(),
    ))
}
