use anyhow::Result;
use clap::Parser;
use env_logger;
use nalgebra::base::*;
use std::time::Instant;

use compact_skel_3d::algorithm::{delaunay_alg, skeleton_alg};
use compact_skel_3d::mesh3d::{self, ManifoldMesh3D};
use compact_skel_3d::skeleton3d;

fn generate_test_mesh() -> Result<ManifoldMesh3D> {
    let mut mesh = ManifoldMesh3D::new();
    let up_vert = mesh.add_vertex(&Vector3::new(0.0, 0.0, 0.5));
    let down_vert = mesh.add_vertex(&Vector3::new(0.0, 0.0, -0.5));
    let mut surr_vert = Vec::new();
    for i in 0..9 {
        let ang = 2.0 * std::f32::consts::PI * (i as f32) / 9.0;
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
    #[arg(long = "objinfile")]
    obj_in_path: Option<std::path::PathBuf>,
    #[arg(default_value = "./ressources/mesh.obj", long = "objoutfile")]
    obj_out_path: std::path::PathBuf,
    #[arg(default_value = "./ressources/skeleton.obj", long = "skeloutfile")]
    skel_out_path: std::path::PathBuf,
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Cli::parse();

    let obj_out_path_str = args.obj_out_path.to_str().unwrap_or("");
    let skel_out_path_str = args.skel_out_path.to_str().unwrap_or("");

    let mut mesh = if let Some(obj_in_path) = args.obj_in_path {
        let obj_in_path_str = obj_in_path.to_str().unwrap_or("");
        mesh3d::io::load_obj_manifold(obj_in_path_str)?
    } else {
        generate_test_mesh()?
    };

    println!("Checking mesh");
    mesh.check_mesh()?;
    println!("");

    println!("Mesh to delaunay");
    let now = Instant::now();
    delaunay_alg::to_delaunay(&mut mesh, Some(std::f32::consts::PI * 20.0 / 180.0))?;
    let duration = now.elapsed();
    let sec = duration.as_secs();
    let min = sec / 60;
    let sec = sec - min * 60;
    println!("Delaunay computed in {}m{}s", min, sec);
    println!("");

    let now = Instant::now();
    println!("Full skeletonization");
    let skeleton = skeleton_alg::full_skeletonization(&mut mesh)?;
    let duration = now.elapsed();
    let sec = duration.as_secs();
    let min = sec / 60;
    let sec = sec - min * 60;
    println!("Skeleton computed in {}m{}s", min, sec);
    println!("");

    println!("Saving skeleton and mesh");
    mesh3d::io::save_obj_manifold(obj_out_path_str, &mesh, None)?;
    skeleton3d::io::save_obj(skel_out_path_str, &skeleton, None)?;

    Ok(())
}
