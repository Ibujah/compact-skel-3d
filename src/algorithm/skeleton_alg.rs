use crate::algorithm::skeleton_interface::{skeleton_operations, SkeletonInterface3D};
use crate::mesh3d::Mesh3D;
use crate::skeleton3d::Skeleton3D;
use anyhow::Result;

pub fn full_skeletonization(mesh: &mut Mesh3D, skeleton: &mut Skeleton3D) -> Result<()> {
    println!("Init skeleton interface");
    let mut skeleton_interface = SkeletonInterface3D::init(mesh, skeleton)?;

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
                skeleton_operations::include_alveola_in_skel(&mut skeleton_interface, ind_alveola)?;
            }
            print!("\r{} alveolae remaining     ", vec_alveola.len());
        } else {
            break;
        }
    }
    println!("");

    println!("Checking skeleton");
    skeleton_interface.check()?;

    Ok(())
}

pub fn sheet_skeletonization(mesh: &mut Mesh3D, skeleton: &mut Skeleton3D) -> Result<()> {
    println!("Init skeleton interface");
    let mut skeleton_interface = SkeletonInterface3D::init(mesh, skeleton)?;

    println!("Finding some first alveola");
    let ind_first_alveola = skeleton_operations::first_alveola_in(&mut skeleton_interface)?;
    let mut vec_alveola = Vec::new();
    vec_alveola.push(ind_first_alveola);

    println!("Propagating sheet");
    let label = 1;
    loop {
        if let Some(ind_alveola) = vec_alveola.pop() {
            let alveola = skeleton_interface.get_alveola(ind_alveola)?;
            let alveola_in = alveola.is_full();
            if alveola_in {
                skeleton_operations::compute_sheet(&mut skeleton_interface, ind_alveola, label)?;
                let current_sheet = skeleton_interface.get_sheet(label);
                skeleton_operations::extract_skeleton_paths(
                    &mut skeleton_interface,
                    &current_sheet,
                )?;
                for &ind_alve in current_sheet.iter() {
                    skeleton_operations::include_alveola_in_skel(
                        &mut skeleton_interface,
                        ind_alve,
                    )?;
                }

                break;
            }
        } else {
            break;
        }
    }

    println!("Checking skeleton");
    skeleton_interface.check()?;

    Ok(())
}
