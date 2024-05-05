use anyhow::Result;
use clap::Parser;
use env_logger;

use std::fs::File;
use std::io::Write;

use compact_skel_3d::mesh3d;

#[derive(Parser)]
struct Cli {
    #[arg(long = "mesh")]
    mesh_path: Option<std::path::PathBuf>,
    #[arg(default_value = "./output/", long = "pathout")]
    out_path: std::path::PathBuf,
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Cli::parse();

    let out_path_str = args.out_path.to_str().unwrap();

    let mut file_manifold = File::create(&format!("{}manifold.txt", out_path_str))?;
    writeln!(file_manifold, "0")?;

    let mesh = if let Some(mesh_path) = args.mesh_path {
        let mesh_path_str = mesh_path.to_str().unwrap_or("");
        let extension = &mesh_path_str[mesh_path_str.len() - 3..];
        if extension == "obj" {
            mesh3d::io::load_obj_manifold(mesh_path_str)?
        } else if extension == "off" {
            mesh3d::io::load_off_manifold(mesh_path_str)?
        } else {
            return Err(anyhow::Error::msg("Extension not handled"));
        }
    } else {
        return Err(anyhow::Error::msg("No file provided"));
    };

    mesh.check_mesh()?;

    let mut file_manifold = File::create(&format!("{}manifold.txt", out_path_str))?;
    writeln!(file_manifold, "1")?;

    Ok(())
}
