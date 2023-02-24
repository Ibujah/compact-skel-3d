use anyhow::Result;
use clap::Parser;
use std::time::Instant;

use compact_skel_3d::mesh3d::io;
use compact_skel_3d::algorithm::delaunay_alg;

#[derive(Parser)]
struct Cli {
    #[arg(default_value="./ressources/hand.obj", long = "objfile")]
    obj_path: std::path::PathBuf,
    #[arg(default_value="./ressources/hand_del.obj", long = "objoutfile")]
    obj_out_path: std::path::PathBuf,
}

fn main() -> Result<()> {
    let args = Cli::parse();
    
    let obj_path_str = args.obj_path.to_str().unwrap_or("");
    let obj_out_path_str = args.obj_out_path.to_str().unwrap_or("");
    
    println!("Loading topo_mesh");
    let mut mesh = io::load_obj(obj_path_str)?;
    
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
    io::save_obj(obj_out_path_str, &mesh)?;
    
    Ok(())
}
