use anyhow::Result;
use std::collections::HashSet;

use crate::algorithm::delaunay_interface::DelaunayInterface;
use crate::mesh3d::{mesh3d, Mesh3D};

fn extract_physical_edges(
    mesh: &Mesh3D,
    ang_max: Option<f32>,
) -> Result<HashSet<mesh3d::HalfEdge>> {
    let ang_max = ang_max.unwrap_or(std::f32::consts::PI);
    let cos_min = ang_max.cos();

    let mut physical: HashSet<mesh3d::HalfEdge> = HashSet::new();
    // set physical edges
    for e in 0..mesh.get_nb_halfedges() {
        let he = mesh.get_halfedge(e)?;
        if he.halfedge()[0] > he.halfedge()[1] {
            continue;
        }

        if ang_max == std::f32::consts::PI {
            physical.insert(he.halfedge());
        }

        // compute angles between adjacent faces
        let face_a = he.face().ok_or(anyhow::Error::msg(
            "extract_physical_edges(): Halfedge should be linked to a face",
        ))?;
        let face_b = he
            .opposite_halfedge()
            .ok_or(anyhow::Error::msg(
                "extract_physical_edges(): Halfedge should have opposite halfedge",
            ))?
            .face()
            .ok_or(anyhow::Error::msg(
                "extract_physical_edges(): Opposite halfedge should be connected to face",
            ))?;

        // getting vertices
        let [vert_a_1, vert_a_2, vert_a_3] = face_a.vertices();
        let [vert_b_1, vert_b_2, vert_b_3] = face_b.vertices();
        let pt_a_1 = vert_a_1.vertex();
        let pt_a_2 = vert_a_2.vertex();
        let pt_a_3 = vert_a_3.vertex();
        let pt_b_1 = vert_b_1.vertex();
        let pt_b_2 = vert_b_2.vertex();
        let pt_b_3 = vert_b_3.vertex();

        // computing normals
        let vec_u_1 = pt_a_2 - pt_a_1;
        let vec_u_2 = pt_b_2 - pt_b_1;
        let vec_v_1 = pt_a_3 - pt_a_1;
        let vec_v_2 = pt_b_3 - pt_b_1;

        let nor_1 = vec_u_1.cross(&vec_v_1);
        let nor_2 = vec_u_2.cross(&vec_v_2);

        // cosinus between normals
        let cos_cur = nor_1.dot(&nor_2).abs();

        if cos_cur > cos_min {
            physical.insert(he.halfedge());
        }
    }
    Ok(physical)
}

fn compute_halfedge_split_vertex(
    deltet: &DelaunayInterface,
    vertex_inds: [usize; 2],
) -> Result<mesh3d::Vertex> {
    let (vert1, vert2) = if deltet.is_original_vertex(vertex_inds[0]) {
        (
            deltet.get_mesh().get_vertex(vertex_inds[0])?.vertex(),
            deltet.get_mesh().get_vertex(vertex_inds[1])?.vertex(),
        )
    } else {
        (
            deltet.get_mesh().get_vertex(vertex_inds[1])?.vertex(),
            deltet.get_mesh().get_vertex(vertex_inds[0])?.vertex(),
        )
    };

    let vert_mid = (vert1 + vert2) * 0.5;
    let dist = (0.5 * (vert1 - vert2)).norm();
    let k_f = dist.log2();
    let k_1 = k_f.floor();
    let k_2 = k_1 + 1.0;

    let vec = (vert2 - vert1).normalize();

    let vs1 = vert1 + vec * 2.0_f32.powf(k_1);
    let vs2 = vert1 + vec * 2.0_f32.powf(k_2);

    let d1 = (vert_mid - vs1).norm();
    let d2 = (vert_mid - vs2).norm();

    if d1 < d2 {
        Ok(vs1)
    } else {
        Ok(vs2)
    }
}

fn compute_face_split_vertex(face: mesh3d::IterFace) -> Result<mesh3d::Vertex> {
    let [vert1, vert2, vert3] = face.vertices();

    Ok((vert1.vertex() + vert2.vertex() + vert3.vertex()) / 3.0)
}

pub fn to_delaunay(mesh: &mut Mesh3D, ang_max: Option<f32>) -> Result<()> {
    let mut deltet = DelaunayInterface::from_mesh(mesh)?;

    let physical = extract_physical_edges(deltet.get_mesh(), ang_max)?;

    let mut nb_non_del_hedges = deltet.count_non_del_halfedges()?;
    let mut nb_non_del_faces = deltet.count_non_del_faces()?;
    println!("Vertices: {}", deltet.get_mesh().get_nb_vertices());
    println!(
        "Non delaunay edges: {}/{}",
        nb_non_del_hedges >> 1,
        deltet.get_mesh().get_nb_halfedges() >> 1
    );
    println!(
        "Non delaunay faces: {}/{}",
        nb_non_del_faces,
        deltet.get_mesh().get_nb_faces()
    );

    let mut num_split_edge = 0;
    let mut num_split_face = 0;
    let mut num_flip = 0;
    let mut shift_edge = 0;
    let mut shift_face = 0;
    let mut cpt_force_split = 0;

    print!(
        "\r{} flip(s), {} edge split(s), {} face split(s)",
        num_flip, num_split_edge, num_split_face
    );

    loop {
        if let Some(he) = deltet.get_non_del_halfedge(Some(shift_edge))? {
            shift_edge = (shift_edge + 1) % deltet.get_mesh().get_nb_halfedges();
            let mut he_inds = he.halfedge();
            he_inds.sort();
            let is_physical = physical.contains(&he_inds);
            let index_he = he.ind();
            let flipped = if !is_physical && cpt_force_split < nb_non_del_hedges {
                deltet.flip_halfedge(index_he)?
            } else {
                false
            };
            if flipped {
                num_flip = num_flip + 1;
                cpt_force_split = cpt_force_split + 1;
            } else {
                let vert_split = compute_halfedge_split_vertex(&deltet, he_inds)?;
                deltet.split_halfedge(&vert_split, index_he)?;
                num_split_edge = num_split_edge + 1;
                cpt_force_split = 0;
            }
        } else if let Some(face) = deltet.get_non_del_face(Some(shift_face))? {
            shift_face = shift_face + 1;
            let vert_split = compute_face_split_vertex(face)?;
            deltet.split_face(&vert_split, face.ind())?;
            num_split_face = num_split_face + 1;
        } else {
            break;
        }
        nb_non_del_hedges = deltet.count_non_del_halfedges()?;
        nb_non_del_faces = deltet.count_non_del_faces()?;
        print!("\r{} non del edges, {} non del faces, {} flip(s), {} edge split(s), {} face split(s)    ", 
               nb_non_del_hedges >> 1, nb_non_del_faces, num_flip, num_split_edge, num_split_face);
    }
    print!("\r{} flip(s), {} edge split(s), {} face split(s)                                                                          ", num_flip, num_split_edge, num_split_face);
    println!("");

    nb_non_del_hedges = deltet.count_non_del_halfedges()?;
    nb_non_del_faces = deltet.count_non_del_faces()?;
    println!("Vertices: {}", deltet.get_mesh().get_nb_vertices());
    println!(
        "Non delaunay edges: {}/{}",
        nb_non_del_hedges >> 1,
        deltet.get_mesh().get_nb_halfedges() >> 1
    );
    println!(
        "Non delaunay faces: {}/{}",
        nb_non_del_faces,
        deltet.get_mesh().get_nb_faces()
    );

    Ok(())
}
