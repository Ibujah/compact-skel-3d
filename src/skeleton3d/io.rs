use anyhow::Result;
use nalgebra::base::*;
use rand::Rng;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, Write};

use ply_rs::parser::Parser;
use ply_rs::ply::{DefaultElement, Property};

use crate::skeleton3d::Skeleton3D;

use super::skeleton3d::Sphere;

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
        let rand_r = rng.gen_range(0..11) as f64;
        let rand_g = rng.gen_range(0..11) as f64;
        let rand_b = rng.gen_range(0..11) as f64;
        let col_r = rand_r / 10.0;
        let col_g = rand_g / 10.0;
        let col_b = rand_b / 10.0;
        writeln!(file, "newmtl mtl_sheet{}", lab)?;
        writeln!(file, "Kd {} {} {}", col_r, col_g, col_b)?;
    }

    Ok(())
}

/// Save radii as .rad file
pub fn save_rad(filename: &str, skeleton: &Skeleton3D) -> Result<()> {
    let mut file = File::create(filename)?;

    for (_, sph) in skeleton.nodes.iter() {
        let rad = sph.radius;
        writeln!(file, "{}", rad)?;
    }

    Ok(())
}

/// Save skeleton as .ply file
pub fn save_ply(
    filename: &str,
    skeleton: &Skeleton3D,
    colors: Option<Vec<[u8; 3]>>,
) -> Result<Vec<[u8; 3]>> {
    let mut file = File::create(filename)?;

    writeln!(file, "ply")?;
    writeln!(file, "format ascii 1.0")?;

    writeln!(file, "element vertex {}", skeleton.nodes.len())?;
    writeln!(file, "property float x")?;
    writeln!(file, "property float y")?;
    writeln!(file, "property float z")?;
    writeln!(file, "property float radius")?;
    writeln!(file, "property uchar red")?;
    writeln!(file, "property uchar green")?;
    writeln!(file, "property uchar blue")?;
    writeln!(file, "property int vertex1")?;
    writeln!(file, "property int vertex2")?;
    writeln!(file, "property int vertex3")?;
    writeln!(file, "property int vertex4")?;

    writeln!(file, "element edge {}", skeleton.edges.len())?;
    writeln!(file, "property int vertex1")?;
    writeln!(file, "property int vertex2")?;

    writeln!(file, "element face {}", skeleton.alveolae.len())?;
    writeln!(file, "property list uchar int vertex_index")?;
    writeln!(file, "property int label")?;
    writeln!(file, "property uchar red")?;
    writeln!(file, "property uchar green")?;
    writeln!(file, "property uchar blue")?;

    writeln!(file, "end_header")?;

    let mut min_rad = -1.0;
    let mut max_rad = -1.0;
    for (_, sph) in skeleton.nodes.iter() {
        let rad = sph.radius;
        if min_rad < 0.0 || min_rad > rad {
            min_rad = rad;
        }
        if max_rad < 0.0 || max_rad < rad {
            max_rad = rad;
        }
    }

    let mut skel_ind_to_ind = HashMap::new();
    let mut ind = 0;
    for (skel_ind, sph) in skeleton.nodes.iter() {
        let vert = sph.center;
        let rad = sph.radius;
        let boundary_ind = skeleton.boundary_inds.get(skel_ind).unwrap();

        let p = (rad - min_rad) / (max_rad - min_rad);
        writeln!(
            file,
            "{} {} {} {} {} {} {} {} {} {} {}",
            vert[0],
            vert[1],
            vert[2],
            rad,
            (p * 255.0) as u8,
            0,
            ((1.0 - p) * 255.0) as u8,
            boundary_ind[0],
            boundary_ind[1],
            boundary_ind[2],
            boundary_ind[3],
        )?;
        skel_ind_to_ind.insert(skel_ind, ind);
        ind = ind + 1;
    }

    let vec_col = if let Some(col) = colors {
        col
    } else {
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
        let mut vec_col = Vec::new();
        let mut rng = rand::thread_rng();
        for _ in 0..(lab_max + 1) {
            let rand_r = rng.gen_range(0..11) as f64;
            let rand_g = rng.gen_range(0..11) as f64;
            let rand_b = rng.gen_range(0..11) as f64;
            let col_r = (255.0 * rand_r / 10.0) as u8;
            let col_g = (255.0 * rand_g / 10.0) as u8;
            let col_b = (255.0 * rand_b / 10.0) as u8;
            vec_col.push([col_r, col_g, col_b]);
        }
        vec_col
    };

    for (_, edge) in skeleton.edges.iter() {
        writeln!(
            file,
            "{} {}",
            skel_ind_to_ind[&edge[0]], skel_ind_to_ind[&edge[1]]
        )?;
    }

    for (alv_ind, alv_nods) in skeleton.alveolae.iter() {
        let label = skeleton.labels[alv_ind];
        write!(file, "{} ", alv_nods.len())?;
        for i in alv_nods {
            write!(file, "{} ", skel_ind_to_ind[i])?;
        }
        if let Some(lab) = label {
            writeln!(
                file,
                "{} {} {} {}",
                lab, vec_col[lab][0], vec_col[lab][1], vec_col[lab][2]
            )?;
        } else {
            let lab = vec_col.len() - 1;
            writeln!(
                file,
                "{} {} {} {}",
                lab, vec_col[lab][0], vec_col[lab][1], vec_col[lab][2]
            )?;
        }
    }

    Ok(vec_col)
}

/// Save problematic edges as .ply file
pub fn save_problematics_ply(
    filename: &str,
    skeleton: &Skeleton3D,
    problematic_edge: &Vec<usize>,
) -> Result<()> {
    let mut file = File::create(filename)?;

    writeln!(file, "ply")?;
    writeln!(file, "format ascii 1.0")?;

    writeln!(file, "element vertex {}", skeleton.nodes.len())?;
    writeln!(file, "property float x")?;
    writeln!(file, "property float y")?;
    writeln!(file, "property float z")?;
    writeln!(file, "property uchar red")?;
    writeln!(file, "property uchar green")?;
    writeln!(file, "property uchar blue")?;

    writeln!(file, "element edge {}", problematic_edge.len())?;
    writeln!(file, "property int vertex1")?;
    writeln!(file, "property int vertex2")?;

    writeln!(file, "end_header")?;

    let mut min_rad = -1.0;
    let mut max_rad = -1.0;
    for (_, sph) in skeleton.nodes.iter() {
        let rad = sph.radius;
        if min_rad < 0.0 || min_rad < rad {
            min_rad = rad;
        }
        if max_rad < 0.0 || max_rad > rad {
            max_rad = rad;
        }
    }

    let mut skel_ind_to_ind = HashMap::new();
    let mut ind = 0;
    for (skel_ind, sph) in skeleton.nodes.iter() {
        let vert = sph.center;

        writeln!(
            file,
            "{} {} {} {} {} {}",
            vert[0], vert[1], vert[2], 255, 0, 0,
        )?;
        skel_ind_to_ind.insert(skel_ind, ind);
        ind = ind + 1;
    }

    for ind_edge in problematic_edge.iter() {
        let edge = skeleton.edges[ind_edge];
        writeln!(
            file,
            "{} {}",
            skel_ind_to_ind[&edge[0]], skel_ind_to_ind[&edge[1]]
        )?;
    }

    Ok(())
}

pub fn load_voxel_core(filename_ply: &str, filename_rad: &str) -> Result<Skeleton3D> {
    let mut vec_vert = Vec::new();
    let mut vec_rad = Vec::new();
    let mut vec_edge = Vec::new();
    let mut vec_face = Vec::new();

    let file_ply = File::open(filename_ply)?;
    let lines_ply = io::BufReader::new(file_ply).lines();
    let mut opt_nb_vert = None;
    let mut opt_nb_edge = None;
    let mut opt_nb_face = None;
    let mut cur_vert = 0;
    let mut cur_edge = 0;
    let mut cur_face = 0;
    let mut header = true;
    for line_ in lines_ply {
        if let Ok(line) = line_ {
            if header {
                if line.len() > 15 && line[0..15].eq("element vertex ") {
                    opt_nb_vert = Some(line[15..].parse::<usize>()?);
                }
                if line.len() > 13 && line[0..13].eq("element edge ") {
                    opt_nb_edge = Some(line[13..].parse::<usize>()?);
                }
                if line.len() > 13 && line[0..13].eq("element face ") {
                    opt_nb_face = Some(line[13..].parse::<usize>()?);
                }
                if line.eq("end_header") {
                    header = false;
                }
            } else if cur_vert < opt_nb_vert.unwrap() {
                let mut line_split = line.split_whitespace();
                let mut vert: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
                for i in 0..3 {
                    let ind = line_split
                        .next()
                        .ok_or(anyhow::Error::msg("Expected value3"))?
                        .parse::<f64>()?;
                    vert[i] = ind;
                }
                vec_vert.push(vert);
                cur_vert += 1;
            } else if cur_edge < opt_nb_edge.unwrap() {
                let line_split: Vec<&str> = line.split_whitespace().collect();
                let mut edge: Vec<usize> = vec![0 as usize; 2];
                for i in 0..2 {
                    let ind = line_split[i].parse::<usize>()?;
                    edge[i] = ind;
                }
                vec_edge.push(edge);
                cur_edge += 1;
            } else if cur_face < opt_nb_face.unwrap() {
                let line_split: Vec<&str> = line.split_whitespace().collect();
                let mut face: Vec<usize> = vec![0 as usize; 3];
                for i in 0..3 {
                    let ind = line_split[line_split.len() + i - 3].parse::<usize>()?;
                    face[i] = ind;
                }
                vec_face.push(face);
                cur_face += 1;
            } else {
                break;
            }
        }
    }
    let file_rad = File::open(filename_rad)?;
    let lines_rad = io::BufReader::new(file_rad).lines();
    let mut header = true;
    for line_ in lines_rad {
        if let Ok(line) = line_ {
            if header {
                header = false;
            } else {
                let line_split: Vec<&str> = line.split_whitespace().collect();
                let rad = line_split[3].parse::<f64>()?;
                vec_rad.push(rad);
            }
        }
    }

    let mut skel = Skeleton3D::new();
    for i in 0..vec_vert.len() {
        let sphere = Sphere {
            center: vec_vert[i],
            radius: vec_rad[i],
        };
        skel.add_sphere(i, sphere)?;
    }
    for i in 0..vec_edge.len() {
        let edg = if vec_edge[i][0] < vec_edge[i][1] {
            [vec_edge[i][0], vec_edge[i][1]]
        } else {
            [vec_edge[i][1], vec_edge[i][0]]
        };
        skel.add_edge(i, edg);
    }
    for i in 0..vec_face.len() {
        skel.add_alveola(i, vec_face[i].clone());
    }

    Ok(skel)
}

pub fn load_sat(filename_moff: &str) -> Result<Skeleton3D> {
    let mut vec_vert = Vec::new();
    let mut vec_rad = Vec::new();
    let mut vec_face = Vec::new();

    let file = File::open(filename_moff)?;
    let lines = io::BufReader::new(file).lines();
    let mut opt_nb_vert = None;
    let mut opt_nb_face = None;
    let mut cur_vert = 0;
    let mut cur_face = 0;
    for line_ in lines {
        if let Ok(line) = line_ {
            if opt_nb_vert.is_none() {
                let mut line_split = line.split_whitespace();
                let moff = line_split
                    .next()
                    .ok_or(anyhow::Error::msg("Expected value1"))?;
                if moff != "MOFF" {
                    return Err(anyhow::Error::msg("Expected MOFF string"));
                }
                let nb_vert = line_split
                    .next()
                    .ok_or(anyhow::Error::msg("Expected value1"))?
                    .parse::<usize>()?;
                opt_nb_vert = Some(nb_vert);
                let nb_face = line_split
                    .next()
                    .ok_or(anyhow::Error::msg("Expected value2"))?
                    .parse::<usize>()?;
                opt_nb_face = Some(nb_face);
            } else {
                let nb_vert = opt_nb_vert.unwrap();
                let nb_face = opt_nb_face.unwrap();
                if cur_vert < nb_vert {
                    let mut line_split = line.split_whitespace();
                    let mut vert: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
                    for i in 0..3 {
                        let ind = line_split
                            .next()
                            .ok_or(anyhow::Error::msg("Expected value3"))?
                            .parse::<f64>()?;
                        vert[i] = ind;
                    }
                    let rad = line_split
                        .next()
                        .ok_or(anyhow::Error::msg("Expected value3"))?
                        .parse::<f64>()?;
                    vec_vert.push(vert);
                    vec_rad.push(rad);

                    cur_vert = cur_vert + 1;
                } else if cur_vert < nb_face {
                    let mut line_split = line.split_whitespace();
                    let mut face = Vec::new();
                    let nbv = line_split
                        .next()
                        .ok_or(anyhow::Error::msg("Expected value4"))?
                        .parse::<usize>()?;
                    for _ in 0..nbv {
                        let ind = line_split
                            .next()
                            .ok_or(anyhow::Error::msg("Expected value5"))?
                            .parse::<usize>()?;
                        face.push(ind);
                    }
                    vec_face.push(face);
                    cur_face = cur_face + 1;
                }
            }
        }
    }

    let mut skel = Skeleton3D::new();
    for i in 0..vec_vert.len() {
        let sphere = Sphere {
            center: vec_vert[i],
            radius: vec_rad[i],
        };
        skel.add_sphere(i, sphere)?;
    }
    for i in 0..vec_face.len() {
        skel.add_alveola(i, vec_face[i].clone());
    }

    Ok(skel)
}

pub fn import_from_ply(file_path: &str) -> Result<Skeleton3D> {
    let mut f = std::fs::File::open(file_path).unwrap();

    let p = Parser::<DefaultElement>::new();
    let ply = p.read_ply(&mut f)?;

    let mut skel = Skeleton3D::new();

    // load vertices
    if !ply.payload.contains_key("vertex") {
        return Err(anyhow::Error::msg("No vertex element in file"));
    }
    let mut ind_nod = 0;
    for v in ply.payload["vertex"].iter() {
        let mut x = None;
        let mut y = None;
        let mut z = None;
        let mut radius = None;
        let mut properties = HashMap::new();

        for (key, prop) in v.into_iter() {
            match (key.as_ref(), prop) {
                ("x", Property::Float(val)) => x = Some(val),
                ("y", Property::Float(val)) => y = Some(val),
                ("z", Property::Float(val)) => z = Some(val),
                ("radius", Property::Float(val)) => radius = Some(val),
                (k, p) => {
                    properties.insert(k.to_string(), p.clone());
                    ()
                }
            }
        }
        let x = *x.ok_or(anyhow::Error::msg("No x property in vertex"))?;
        let y = *y.ok_or(anyhow::Error::msg("No y property in vertex"))?;
        let z = *z.ok_or(anyhow::Error::msg("No z property in vertex"))?;
        let radius = *radius.ok_or(anyhow::Error::msg("No radius property in vertex"))?;
        let sphere = Sphere {
            center: Vector3::new(x.into(), y.into(), z.into()),
            radius: radius.into(),
        };
        skel.add_sphere(ind_nod, sphere)?;
        ind_nod = ind_nod + 1;
    }

    // load faces
    if !ply.payload.contains_key("face") {
        return Err(anyhow::Error::msg("No face element in file"));
    }
    let mut ind_alv = 0;
    for f in ply.payload["face"].iter() {
        let mut list_vertices = None;
        let mut properties = HashMap::new();

        for (key, prop) in f.into_iter() {
            match (key.as_ref(), prop) {
                ("vertex_index", Property::ListInt(val)) => {
                    list_vertices = Some(val.iter().map(|&v| usize::try_from(v).unwrap()).collect())
                }
                (k, p) => {
                    properties.insert(k.to_string(), p.clone());
                    ()
                }
            }
        }

        let list_vertices: Vec<usize> =
            list_vertices.ok_or(anyhow::Error::msg("No vertex_index property in face"))?;

        skel.add_alveola(ind_alv, list_vertices);
        ind_alv = ind_alv + 1;
    }

    Ok(skel)
}
