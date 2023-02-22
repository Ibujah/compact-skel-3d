use anyhow::Result;
use std::collections::HashSet;

use crate::boundary3d::topo_mesh::{self, TopoMesh};
use crate::boundary3d::delaunay_struct::DelaunayStruct;

fn extract_physical_edges(topomesh: &TopoMesh, ang_max: Option<f32>) -> Result<HashSet<topo_mesh::HalfEdge>> {
    let ang_max = ang_max.unwrap_or(std::f32::consts::PI);
    let cos_min = ang_max.cos();

    let mut physical : HashSet<topo_mesh::HalfEdge> = HashSet::new();
    // set physical edges
    for e in 0..topomesh.get_nb_halfedges() {
        let he = topomesh.get_halfedge(e)?;
        if he[0] > he[1] {
            continue
        }
        
        if ang_max == std::f32::consts::PI {
            physical.insert(he);
        }
        
        // compute angles between adjacent faces
        let ind_fac1 = topomesh.get_associated_face(e)?;
        let ind_fac2 = topomesh.get_associated_face(topomesh.get_opposite_halfedge(e)?)?;
        
        // getting vertices
        let vert1 = topomesh.get_face_vertices(ind_fac1)?; 
        let pt_a_1 = topomesh.get_vertex(vert1[0])?;
        let pt_b_1 = topomesh.get_vertex(vert1[1])?;
        let pt_c_1 = topomesh.get_vertex(vert1[2])?;
        let vert2 = topomesh.get_face_vertices(ind_fac2)?; 
        let pt_a_2 = topomesh.get_vertex(vert2[0])?;
        let pt_b_2 = topomesh.get_vertex(vert2[1])?;
        let pt_c_2 = topomesh.get_vertex(vert2[2])?;
        
        // computing normals
        let vec_u_1 = pt_b_1 - pt_a_1;
        let vec_u_2 = pt_b_2 - pt_a_2;
        let vec_v_1 = pt_c_1 - pt_a_1;
        let vec_v_2 = pt_c_2 - pt_a_2;
        
        let nor_1 = vec_u_1.cross(&vec_v_1);
        let nor_2 = vec_u_2.cross(&vec_v_2);
        
        // cosinus between normals
        let cos_cur = nor_1.dot(&nor_2).abs();
        
        if cos_cur > cos_min {
            physical.insert(he);
        }
    }
    Ok(physical)
}

fn compute_halfedge_split_vertex(deltet: &DelaunayStruct, halfedge: [usize; 2]) -> Result<topo_mesh::Vertex> {
    let (vert1, vert2) = 
        if deltet.is_original_vertex(halfedge[0]) {
            (deltet.get_topomesh().get_vertex(halfedge[0])?,
            deltet.get_topomesh().get_vertex(halfedge[1])?)
        }
        else {
            (deltet.get_topomesh().get_vertex(halfedge[1])?,
            deltet.get_topomesh().get_vertex(halfedge[0])?)
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
    }
    else {
        Ok(vs2)
    }
}

fn compute_face_split_vertex(deltet: &DelaunayStruct, face: [usize; 3]) -> Result<topo_mesh::Vertex> {
    let vert1 = deltet.get_topomesh().get_vertex(face[0])?;
    let vert2 = deltet.get_topomesh().get_vertex(face[1])?;
    let vert3 = deltet.get_topomesh().get_vertex(face[2])?;

    Ok((vert1 + vert2 + vert3) / 3.0)
}

pub fn to_delaunay(topomesh: &mut TopoMesh, ang_max: Option<f32>) -> Result<()> {
    let mut deltet = DelaunayStruct::init(topomesh)?;

    let physical = extract_physical_edges(deltet.get_topomesh(), ang_max)?;
    
    let mut nb_non_del_hedges = deltet.count_non_del_halfedges()?;
    let mut nb_non_del_faces = deltet.count_non_del_faces()?;
    println!("Vertices: {}", deltet.get_topomesh().get_nb_vertices());
    println!("Non delaunay edges: {}/{}", nb_non_del_hedges >> 1, deltet.get_topomesh().get_nb_halfedges()>>1);
    println!("Non delaunay faces: {}/{}", nb_non_del_faces, deltet.get_topomesh().get_nb_faces());

    let mut num_split_edge = 0;
    let mut num_split_face = 0;
    let mut num_flip = 0;
    let mut shift_edge = 0;
    let mut shift_face = 0;
    let mut cpt_force_split = 0;

    print!("\r{} flip(s), {} edge split(s), {} face split(s)", num_flip, num_split_edge, num_split_face);

    loop {
        if let Some((ind_he, mut he)) = deltet.get_non_del_halfedge(Some(shift_edge))? {
            shift_edge = (shift_edge+1)%deltet.get_topomesh().get_nb_halfedges();
            he.sort();
            let is_physical = physical.contains(&he);
            let flipped = 
                if !is_physical && cpt_force_split < nb_non_del_hedges {
                    deltet.flip_halfedge(ind_he)?
                }
                else {
                    false
                };
            if flipped {
                num_flip = num_flip + 1;
                cpt_force_split = cpt_force_split + 1;
            }
            else {
                let vert_split = compute_halfedge_split_vertex(&deltet, he)?;
                deltet.split_halfedge(&vert_split, ind_he)?;
                num_split_edge = num_split_edge + 1;
                cpt_force_split = 0;
            }
        }
        else if let Some((ind_face, face)) = deltet.get_non_del_face(Some(shift_face))? {
            shift_face = shift_face+1;
            let vert_split = compute_face_split_vertex(&deltet, face)?;
            deltet.split_face(&vert_split, ind_face)?;
            num_split_face = num_split_face + 1;
        }
        else {
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
    println!("Vertices: {}", deltet.get_topomesh().get_nb_vertices());
    println!("Non delaunay edges: {}/{}", nb_non_del_hedges>>1, deltet.get_topomesh().get_nb_halfedges()>>1);
    println!("Non delaunay faces: {}/{}", nb_non_del_faces, deltet.get_topomesh().get_nb_faces());

    Ok(())
}
