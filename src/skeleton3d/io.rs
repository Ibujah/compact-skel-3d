use anyhow::Result;
use rand::Rng;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;

use crate::skeleton3d::Skeleton3D;

fn write_alveola(
    file: &mut File,
    skel_ind_to_ind: &HashMap<&usize, i32>,
    alv: &Vec<usize>,
) -> Result<()> {
    for i in 1..(alv.len() >> 1) {
        writeln!(
            file,
            "f {}// {}// {}//",
            skel_ind_to_ind[&alv[alv.len() - i]],
            skel_ind_to_ind[&alv[i - 1]],
            skel_ind_to_ind[&alv[i]],
        )?;
        writeln!(
            file,
            "f {}// {}// {}//",
            skel_ind_to_ind[&alv[alv.len() - i - 1]],
            skel_ind_to_ind[&alv[alv.len() - i]],
            skel_ind_to_ind[&alv[i]],
        )?;
    }
    if alv.len() % 2 == 1 {
        let ind = alv.len() >> 1;
        writeln!(
            file,
            "f {}// {}// {}//",
            skel_ind_to_ind[&alv[ind - 1]],
            skel_ind_to_ind[&alv[ind]],
            skel_ind_to_ind[&alv[ind + 1]],
        )?;
    }
    Ok(())
}

/// Save skeleton as .obj file
pub fn save_obj(
    filename: &str,
    skeleton: &Skeleton3D,
    opt_material_file: Option<&str>,
) -> Result<()> {
    let mut file = File::create(filename)?;

    if let Some(material_file) = opt_material_file {
        writeln!(file, "mtllib {}", material_file)?;
    }

    let mut skel_ind_to_ind = HashMap::new();
    let mut ind = 1;
    for (skel_ind, sph) in skeleton.nodes.iter() {
        let vert = sph.center;
        writeln!(file, "v {} {} {}", vert[0], vert[1], vert[2])?;
        skel_ind_to_ind.insert(skel_ind, ind);
        ind = ind + 1;
    }

    let lab_max = skeleton.labels.iter().fold(0, |lm, (_, opt_lab)| {
        if let &Some(lab) = opt_lab {
            if lm > lab {
                lm
            } else {
                lab
            }
        } else {
            lm
        }
    }) + 1;

    let alv_ind_none: Vec<usize> = skeleton
        .labels
        .iter()
        .filter_map(
            |(&ind, &opt_lab)| {
                if opt_lab.is_some() {
                    None
                } else {
                    Some(ind)
                }
            },
        )
        .collect();
    if alv_ind_none.len() != 0 {
        writeln!(file, "g sheet_no_label")?;
        for ind in alv_ind_none.iter() {
            write_alveola(&mut file, &skel_ind_to_ind, &skeleton.alveolae[ind])?;
        }
    }

    for lab_curr in 0..lab_max {
        let alv_ind: Vec<usize> = skeleton
            .labels
            .iter()
            .filter_map(|(&ind, &opt_lab)| {
                if let Some(lab) = opt_lab {
                    if lab == lab_curr {
                        Some(ind)
                    } else {
                        None
                    }
                } else {
                    None
                }
            })
            .collect();
        if alv_ind.len() != 0 {
            writeln!(file, "g sheet{}", lab_curr)?;
            if opt_material_file.is_some() {
                writeln!(file, "usemtl mtl_sheet{}", lab_curr)?;
            }
            for ind in alv_ind.iter() {
                write_alveola(&mut file, &skel_ind_to_ind, &skeleton.alveolae[ind])?;
            }
        }
    }

    Ok(())
}

/// Save material file
pub fn save_mtl(filename: &str, skeleton: &Skeleton3D) -> Result<()> {
    let mut file = File::create(filename)?;

    let mut siz_sheet = HashMap::new();
    for (_, &opt_lab) in skeleton.labels.iter() {
        if let Some(lab) = opt_lab {
            siz_sheet
                .entry(lab)
                .and_modify(|e| *e = *e + 1)
                .or_insert(1);
        }
    }

    let mut rng = rand::thread_rng();
    for (&lab, _) in siz_sheet.iter() {
        let rand_r = rng.gen_range(0..11) as f32;
        let rand_g = rng.gen_range(0..11) as f32;
        let rand_b = rng.gen_range(0..11) as f32;
        let col_r = rand_r / 10.0;
        let col_g = rand_g / 10.0;
        let col_b = rand_b / 10.0;
        writeln!(file, "newmtl mtl_sheet{}", lab)?;
        writeln!(file, "Kd {} {} {}", col_r, col_g, col_b)?;
    }

    Ok(())
}
