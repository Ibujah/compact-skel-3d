use anyhow::Result;
use clap::Parser;
use env_logger;
use nalgebra::base::*;
use std::fs;
use std::time::Instant;

use std::fs::File;
use std::io::Write;

use compact_skel_3d::algorithm::skeleton_alg;
use compact_skel_3d::mesh3d::{self, ManifoldMesh3D};
use compact_skel_3d::skeleton3d;

fn generate_test_mesh() -> Result<ManifoldMesh3D> {
    let mut mesh = ManifoldMesh3D::new();
    let up_vert = mesh.add_vertex(&Vector3::new(0.0, 0.0, 0.5));
    let down_vert = mesh.add_vertex(&Vector3::new(0.0, 0.0, -0.5));
    let mut surr_vert = Vec::new();
    for i in 0..9 {
        let ang = 2.0 * std::f64::consts::PI * (i as f64) / 9.0;
        let cur_vert = mesh.add_vertex(&Vector3::new(ang.cos(), ang.sin(), 0.0));
        surr_vert.push(cur_vert);
    }

    for i in 0..9 {
        let j = (i + 1) % 9;
        mesh.add_face(up_vert, surr_vert[i], surr_vert[j])?;
        mesh.add_face(down_vert, surr_vert[j], surr_vert[i])?;
    }

    Ok(mesh)
}

#[derive(Parser)]
struct Cli {
    #[arg(long = "meshinfile")]
    mesh_in_path: Option<std::path::PathBuf>,
    #[arg(long = "epsilon")]
    epsilon: Option<f64>,
    #[arg(default_value = "./output/", long = "pathout")]
    out_path: std::path::PathBuf,
    #[arg(default_value = "mesh.ply", long = "meshoutfile")]
    mesh_out_name: std::path::PathBuf,
    #[arg(default_value = "skeleton.ply", long = "skeloutfile")]
    skel_out_name: std::path::PathBuf,
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Cli::parse();

    let mut mesh = if let Some(mesh_in_path) = args.mesh_in_path {
        let mesh_in_path_str = mesh_in_path.to_str().unwrap_or("");
        let extension = &mesh_in_path_str[mesh_in_path_str.len() - 3..];
        if extension == "obj" {
            mesh3d::io::load_obj_manifold(mesh_in_path_str)?
        } else if extension == "off" {
            mesh3d::io::load_off_manifold(mesh_in_path_str)?
        } else {
            return Err(anyhow::Error::msg("Extension not handled"));
        }
    } else {
        generate_test_mesh()?
    };

    let epsilon = args.epsilon;
    let out_path_str = args.out_path.to_str().unwrap();
    let mesh_out_name_str = args.mesh_out_name.to_str().unwrap();
    let skel_out_name_str = args.skel_out_name.to_str().unwrap();

    println!("Checking mesh");
    mesh.check_mesh()?;
    println!("");

    let epsilon = if let Some(val) = epsilon {
        let (bb_min, bb_max) = mesh
            .vertices()
            .iter()
            .fold(
                None,
                |bb_val: Option<(Vector3<f64>, Vector3<f64>)>, (_, &vert)| {
                    if let Some((mut bb_min, mut bb_max)) = bb_val {
                        if vert[0] < bb_min[0] {
                            bb_min[0] = vert[0];
                        }
                        if vert[1] < bb_min[1] {
                            bb_min[1] = vert[1];
                        }
                        if vert[2] < bb_min[2] {
                            bb_min[2] = vert[2];
                        }
                        if vert[0] > bb_max[0] {
                            bb_max[0] = vert[0];
                        }
                        if vert[1] > bb_max[1] {
                            bb_max[1] = vert[1];
                        }
                        if vert[2] > bb_max[2] {
                            bb_max[2] = vert[2];
                        }
                        Some((bb_min, bb_max))
                    } else {
                        Some((vert, vert))
                    }
                },
            )
            .ok_or(anyhow::Error::msg("No point in mesh"))?;

        let length = (bb_min - bb_max).norm();

        println!("Epsilon: {}% of diagonal = {}", val * 100.0, val * length);
        Some(val * length)
    } else {
        None
    };

    let nb_vert = mesh.get_nb_vertices();

    let now = Instant::now();
    println!("Sheet skeletonization");
    let (skeleton, _work_mesh, vec_debug_meshes, problematic_edges) =
        skeleton_alg::sheet_skeletonization(&mut mesh, epsilon)?;
    let duration = now.elapsed();
    let sec_all = duration.as_secs();
    let min = sec_all / 60;
    let sec = sec_all - min * 60;
    println!("Skeleton computed in {}m{}s", min, sec);
    println!("");

    println!("Saving skeleton and debug meshes");
    fs::create_dir_all(out_path_str)?;
    for i in 0..vec_debug_meshes.len() {
        mesh3d::io::save_obj_generic(
            &format!("{}debug{}.obj", out_path_str, i),
            &vec_debug_meshes[i],
        )?;
    }
    let vec_col = skeleton3d::io::save_ply(
        &format!("{}{}", out_path_str, skel_out_name_str),
        &skeleton,
        None,
    )?;
    mesh3d::io::save_ply_manifold(
        &format!("{}{}", out_path_str, mesh_out_name_str),
        &mesh,
        Some(vec_col),
    )?;
    skeleton3d::io::save_problematics_ply(
        &format!("{}problematics.ply", out_path_str),
        &skeleton,
        &problematic_edges,
    )?;

    let mut file_pb = File::create(&format!("{}problematics.txt", out_path_str))?;
    writeln!(file_pb, "{}", problematic_edges.len())?;

    let mut file_time = File::create(&format!("{}time.txt", out_path_str))?;
    writeln!(file_time, "{}", sec_all)?;

    let mut file_vert = File::create(&format!("{}nb_vert.txt", out_path_str))?;
    writeln!(file_vert, "{}", nb_vert)?;

    Ok(())
}
