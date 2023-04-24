use anyhow::Result;
use clap::Parser;
use std::time::Instant;

use compact_skel_3d::algorithm::delaunay_alg;
use compact_skel_3d::mesh3d::io;

#[derive(Parser)]
struct Cli {
    #[arg(default_value = "./ressources/hand.obj", long = "objinfile")]
    obj_in_path: std::path::PathBuf,
    #[arg(default_value = "./ressources/hand_del.obj", long = "objoutfile")]
    obj_out_path: std::path::PathBuf,
}

fn main() -> Result<()> {
    let args = Cli::parse();

    let obj_in_path_str = args.obj_in_path.to_str().unwrap_or("");
    let obj_out_path_str = args.obj_out_path.to_str().unwrap_or("");

    println!("Loading topo_mesh");
    let extension = &obj_in_path_str[obj_in_path_str.len() - 3..];
    println!("{}", extension);

    let mut mesh = if extension == "obj" {
        io::load_obj_manifold(obj_in_path_str)?
    } else if extension == "off" {
        io::load_off_manifold(obj_in_path_str)?
    } else {
        return Err(anyhow::Error::msg("Extension not handled"));
    };

    println!("Checking mesh");
    mesh.check_mesh()?;

    println!("Mesh to delaunay");
    let now = Instant::now();
    delaunay_alg::to_delaunay(&mut mesh, Some(std::f32::consts::PI * 20.0 / 180.0))?;
    let duration = now.elapsed();
    let sec = duration.as_secs();
    let min = sec / 60;
    let sec = sec - min * 60;
    println!("Delaunay mesh computed in {}m{}s", min, sec);

    println!("Checking mesh");
    mesh.check_mesh()?;

    println!("Save mesh");
    io::save_obj_manifold(obj_out_path_str, &mesh, None)?;

    Ok(())
}
