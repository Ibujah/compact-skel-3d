use anyhow::Result;
use nalgebra::base::*;
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
                    let mut vert: Vector3<f32> = Vector3::new(0.0, 0.0, 0.0);
                    line_split.next();
                    for i in 0..3 {
                        let cur = line_split
                            .next()
                            .ok_or(anyhow::Error::msg("Expected value"))?;
                        vert[i] = cur.parse::<f32>()?;
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
                    let mut vert: Vector3<f32> = Vector3::new(0.0, 0.0, 0.0);
                    for i in 0..3 {
                        let ind = line_split
                            .next()
                            .ok_or(anyhow::Error::msg("Expected value3"))?
                            .parse::<f32>()?;
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
pub fn save_obj_manifold(filename: &str, mesh: &ManifoldMesh3D) -> Result<()> {
    let mut file = File::create(filename)?;

    let mut corresp: HashMap<usize, usize> = HashMap::new();
    let mut cpt = 1;

    for v in mesh.vertex_indices() {
        let vert = mesh.get_vertex(v)?.vertex();
        corresp.insert(v, cpt);
        cpt = cpt + 1;
        writeln!(file, "v {} {} {}", vert[0], vert[1], vert[2])?;
    }

    for (&f, _) in mesh.faces() {
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
