use anyhow::Result;
use nalgebra::base::*;
use rand::Rng;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, Write};

use crate::mesh3d::GenericMesh3D;
use crate::mesh3d::ManifoldMesh3D;

/// Loads obj file as manifold mesh
pub fn load_obj_manifold(filename: &str) -> Result<ManifoldMesh3D> {
    let mut mesh = ManifoldMesh3D::new();

    let file = File::open(filename)?;
    let lines = io::BufReader::new(file).lines();
    for line_ in lines {
        if let Ok(line) = line_ {
            if line.len() > 2 {
                if &line[..2] == "v " {
                    let mut line_split = line.split_whitespace();
                    let mut vert: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
                    line_split.next();
                    for i in 0..3 {
                        let cur = line_split
                            .next()
                            .ok_or(anyhow::Error::msg("Expected value"))?;
                        vert[i] = cur.parse::<f64>()?;
                    }

                    mesh.add_vertex(&vert);
                }
                if &line[..2] == "f " {
                    let mut line_split = line.split_whitespace();
                    let mut face: [usize; 3] = [0, 0, 0];
                    line_split.next();
                    for i in 0..3 {
                        let cur = line_split
                            .next()
                            .ok_or(anyhow::Error::msg("Expected value"))?;
                        let mut cur_split = cur.split('/');
                        let ind = cur_split
                            .next()
                            .ok_or(anyhow::Error::msg("Expected value"))?;
                        face[i] = ind.parse::<usize>()? - 1;
                    }

                    mesh.add_face(face[0], face[1], face[2])?;
                }
            }
        }
    }

    Ok(mesh)
}

/// Loads off file as manifold mesh
pub fn load_off_manifold(filename: &str) -> Result<ManifoldMesh3D> {
    let mut mesh = ManifoldMesh3D::new();

    let file = File::open(filename)?;
    let lines = io::BufReader::new(file).lines();
    let mut opt_nb_vert = None;
    let mut opt_nb_face = None;
    let mut cur_vert = 0;
    let mut cur_face = 0;
    for line_ in lines {
        if let Ok(line) = line_ {
            if opt_nb_vert.is_none() {
                if line == "OFF" {
                    continue;
                }
                let mut line_split = line.split_whitespace();
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

                    mesh.add_vertex(&vert);
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
                    if face.len() != 3 {
                        return Err(anyhow::Error::msg("Face with more than 3 vertices"));
                    }
                    mesh.add_face(face[0], face[1], face[2])?;
                    cur_face = cur_face + 1;
                }
            }
        }
    }

    Ok(mesh)
}

/// Save manifold mesh as obj file
pub fn save_obj_manifold(
    filename: &str,
    mesh: &ManifoldMesh3D,
    opt_material_file: Option<&str>,
) -> Result<()> {
    let mut file = File::create(filename)?;

    if let Some(material_file) = opt_material_file {
        writeln!(file, "mtllib {}", material_file)?;
    }

    let mut corresp: HashMap<usize, usize> = HashMap::new();
    let mut cpt = 1;

    for v in mesh.vertex_indices() {
        let vert = mesh.get_vertex(v)?.vertex();
        corresp.insert(v, cpt);
        cpt = cpt + 1;
        writeln!(file, "v {} {} {}", vert[0], vert[1], vert[2])?;
    }

    let mut groups: HashMap<String, Vec<usize>> = HashMap::new();
    let mut non_grouped = Vec::new();
    for (&ind_face, opt_lab) in mesh.groups.iter() {
        if let Some(lab) = opt_lab {
            groups
                .entry(format!("sheet{}", lab))
                .and_modify(|g| g.push(ind_face))
                .or_insert(vec![ind_face]);
        } else {
            non_grouped.push(ind_face);
        }
    }

    for &f in non_grouped.iter() {
        let face = mesh.get_face(f)?.vertices_inds();
        let ind0 = corresp.get(&face[0]).ok_or(anyhow::Error::msg(
            "save_obj(): vertex face does not exists",
        ))?;
        let ind1 = corresp.get(&face[1]).ok_or(anyhow::Error::msg(
            "save_obj(): vertex face does not exists",
        ))?;
        let ind2 = corresp.get(&face[2]).ok_or(anyhow::Error::msg(
            "save_obj(): vertex face does not exists",
        ))?;
        writeln!(file, "f {}// {}// {}//", ind0, ind1, ind2)?;
    }
    for (lab, group) in groups {
        writeln!(file, "g {}", lab)?;
        if opt_material_file.is_some() {
            writeln!(file, "usemtl mtl_{}", lab)?;
        }
        for &f in group.iter() {
            let face = mesh.get_face(f)?.vertices_inds();
            let ind0 = corresp.get(&face[0]).ok_or(anyhow::Error::msg(
                "save_obj(): vertex face does not exists",
            ))?;
            let ind1 = corresp.get(&face[1]).ok_or(anyhow::Error::msg(
                "save_obj(): vertex face does not exists",
            ))?;
            let ind2 = corresp.get(&face[2]).ok_or(anyhow::Error::msg(
                "save_obj(): vertex face does not exists",
            ))?;
            writeln!(file, "f {}// {}// {}//", ind0, ind1, ind2)?;
        }
    }

    Ok(())
}

/// Save non manifold mesh as obj file
pub fn save_obj_generic(filename: &str, mesh: &GenericMesh3D) -> Result<()> {
    let mut file = File::create(filename)?;

    for v in 0..mesh.get_nb_vertices() {
        let vert = mesh.get_vertex(v)?;
        writeln!(file, "v {} {} {}", vert[0], vert[1], vert[2])?;
    }

    for f in 0..mesh.get_nb_faces() {
        let face = mesh.get_face(f)?;
        writeln!(
            file,
            "f {}// {}// {}//",
            face[0] + 1,
            face[1] + 1,
            face[2] + 1
        )?;
    }

    Ok(())
}

/// Save manifold mesh as ply file
pub fn save_ply_manifold(
    filename: &str,
    mesh: &ManifoldMesh3D,
    colors: Option<Vec<[u8; 3]>>,
) -> Result<Vec<[u8; 3]>> {
    let mut file = File::create(filename)?;

    writeln!(file, "ply")?;
    writeln!(file, "format ascii 1.0")?;

    writeln!(file, "element vertex {}", mesh.vertices.len())?;
    writeln!(file, "property float x")?;
    writeln!(file, "property float y")?;
    writeln!(file, "property float z")?;

    writeln!(file, "element face {}", mesh.faces.len())?;
    writeln!(file, "property list uchar int vertex_index")?;
    writeln!(file, "property uchar label")?;
    writeln!(file, "property uchar red")?;
    writeln!(file, "property uchar green")?;
    writeln!(file, "property uchar blue")?;

    writeln!(file, "end_header")?;

    let mut corresp: HashMap<usize, usize> = HashMap::new();
    let mut cpt = 0;

    for v in mesh.vertex_indices() {
        let vert = mesh.get_vertex(v)?.vertex();
        corresp.insert(v, cpt);
        cpt = cpt + 1;
        writeln!(file, "{} {} {}", vert[0], vert[1], vert[2])?;
    }

    let vec_col = if let Some(col) = colors {
        col
    } else {
        let lab_max = mesh.groups.iter().fold(0, |lm, (_, opt_lab)| {
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

    for (&fac_ind, _) in mesh.faces.iter() {
        let face = mesh.get_face(fac_ind)?.vertices_inds();
        let label = mesh.groups[&fac_ind];
        write!(file, "{} ", face.len())?;
        for i in face {
            write!(file, "{} ", corresp[&i])?;
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
