use anyhow::Result;
use clap::Parser;

use std::fs::File;
use std::io::Write;

use compact_skel_3d::algorithm::cc_count;
use compact_skel_3d::skeleton3d;

#[derive(Parser)]
struct Cli {
    #[arg(long = "skel_ply_in")]
    skel_ply_path: std::path::PathBuf,
    #[arg(default_value = "./output/", long = "pathout")]
    out_path: std::path::PathBuf,
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Cli::parse();

    let skel_ply_path = args.skel_ply_path.to_str().unwrap();
    let out_path_str = args.out_path.to_str().unwrap();

    let skel = skeleton3d::io::import_from_ply(skel_ply_path)?;

    let mut vec_lab: Vec<usize> = skel.get_labels().iter().map(|(_, &l)| l.unwrap()).collect();
    vec_lab.sort();
    vec_lab.dedup();
    let nb_sheets = vec_lab.len();
    let nb_vert = skel.get_nodes().len();
    let nb_alv = skel.get_alveolae().len();
    let nb_cc = cc_count::connected_components_count(&skel);

    let mut file_nbsheets = File::create(&format!("{}nbsheets.txt", out_path_str))?;
    writeln!(file_nbsheets, "{}", nb_sheets)?;

    let mut file_nbvert = File::create(&format!("{}nbvert.txt", out_path_str))?;
    writeln!(file_nbvert, "{}", nb_vert)?;

    let mut file_nbalv = File::create(&format!("{}nbalv.txt", out_path_str))?;
    writeln!(file_nbalv, "{}", nb_alv)?;

    let mut file_nbcc = File::create(&format!("{}nbcc.txt", out_path_str))?;
    writeln!(file_nbcc, "{}", nb_cc)?;

    Ok(())
}
