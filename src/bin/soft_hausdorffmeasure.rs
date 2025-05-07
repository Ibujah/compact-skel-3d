use anyhow::Result;
use clap::Parser;

use std::fs::File;
use std::io::Write;

use compact_skel_3d::algorithm::hausdorff_distance;
use compact_skel_3d::mesh3d;
use compact_skel_3d::skeleton3d;

#[derive(Parser)]
struct Cli {
    #[arg(long = "skel_ply_in")]
    skel_ply_path: std::path::PathBuf,
    #[arg(long = "mesh_in")]
    mesh_in_path: std::path::PathBuf,
    #[arg(default_value = "./output/", long = "pathout")]
    out_path: std::path::PathBuf,
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Cli::parse();

    let skel_ply_path = args.skel_ply_path.to_str().unwrap();
    let mesh_in_path_str = args.mesh_in_path.to_str().unwrap();
    let out_path_str = args.out_path.to_str().unwrap();

    let skel = skeleton3d::io::import_from_ply(skel_ply_path)?;

    let extension = &mesh_in_path_str[mesh_in_path_str.len() - 3..];
    let mesh = if extension == "obj" {
        mesh3d::io::load_obj_manifold(mesh_in_path_str)?
    } else if extension == "off" {
        mesh3d::io::load_off_manifold(mesh_in_path_str)?
    } else {
        return Err(anyhow::Error::msg("Extension not handled"));
    };

    let dist = hausdorff_distance::hausdorff_distance(&mesh, &skel);

    let mut file_hdist = File::create(&format!("{}hdist.txt", out_path_str))?;
    writeln!(file_hdist, "{}", dist)?;

    Ok(())
}
