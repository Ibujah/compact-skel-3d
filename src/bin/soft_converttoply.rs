use anyhow::{anyhow, Result};
use clap::Parser;
use env_logger;

use compact_skel_3d::algorithm::skeleton_region_grow::compute_regions;
use compact_skel_3d::skeleton3d;

#[derive(Parser)]
struct Cli {
    #[arg(long = "skel_ply_in")]
    skel_ply_in_path: Option<std::path::PathBuf>,
    #[arg(long = "rad_r_in")]
    rad_r_in_path: Option<std::path::PathBuf>,
    #[arg(long = "skel_moff_in")]
    skel_moff_in_path: Option<std::path::PathBuf>,
    #[arg(default_value = "./output/", long = "pathout")]
    out_path: std::path::PathBuf,
    #[arg(default_value = "skeleton.ply", long = "skeloutfile")]
    skel_out_name: std::path::PathBuf,
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Cli::parse();

    println!("Load skeleton");
    let mut skel = if let (Some(skel_ply_in_path), Some(rad_r_in_path)) =
        (args.skel_ply_in_path, args.rad_r_in_path)
    {
        skeleton3d::io::load_voxel_core(
            skel_ply_in_path.to_str().unwrap(),
            rad_r_in_path.to_str().unwrap(),
        )?
    } else if let Some(skel_moff_in_path) = args.skel_moff_in_path {
        skeleton3d::io::load_sat(skel_moff_in_path.to_str().unwrap())?
    } else {
        return Err(anyhow::Error::msg("no file given"));
    };

    println!("Compute regions");
    compute_regions(&mut skel)?;

    println!("Save skeleton");
    let skel_out_name_str = args.skel_out_name.to_str().unwrap();
    let out_path_str = args.out_path.to_str().unwrap();

    skeleton3d::io::save_ply(
        &format!("{}{}", out_path_str, skel_out_name_str),
        &skel,
        None,
    )?;

    Ok(())
}
