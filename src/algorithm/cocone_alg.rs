use anyhow::Result;
use core::f64;
use nalgebra::Vector3;
use std::collections::HashMap;

use crate::mesh3d::{generic_mesh3d::Vertex, ManifoldMesh3D};

use super::sub_algorithms::DelaunayInterface;

/// Compute the cocone triangles of a mesh.
pub fn cocone_faces(mesh: &mut ManifoldMesh3D) -> Result<HashMap<[usize; 3], bool>> {
    // Compute the Delaunay graph from the mesh vertices.
    let deltet = DelaunayInterface::from_mesh(mesh)?;

    // Compute the cocone triangles of the mesh.
    let tri_cocone = deltet.compute_cocone()?;

    Ok(tri_cocone)
}

/// This function extracts edges from a set of triangles.
fn get_edges_from_triangles(triangles: &Vec<[usize; 3]>) -> HashMap<[usize; 2], Vec<[usize; 3]>> {
    // Initialize an empty hashmap to store the edges.
    let mut edges: HashMap<[usize; 2], Vec<[usize; 3]>> = HashMap::new();

    // Iterate over each triangle.
    for &triangle in triangles.iter() {
        // Extract the vertices of the current triangle.
        let mut edge1 = [triangle[0], triangle[1]];
        let mut edge2 = [triangle[1], triangle[2]];
        let mut edge3 = [triangle[2], triangle[0]];

        // Sort the vertices of each edge to ensure consistent ordering.
        edge1.sort();
        edge2.sort();
        edge3.sort();

        // Add the current triangle to the list of triangles associated with each edge.
        edges
            .entry(edge1)
            .and_modify(|v| v.push(triangle))
            .or_insert(vec![triangle]);
        edges
            .entry(edge2)
            .and_modify(|v| v.push(triangle))
            .or_insert(vec![triangle]);
        edges
            .entry(edge3)
            .and_modify(|v| v.push(triangle))
            .or_insert(vec![triangle]);
    }

    // Return the map of edges to their associated triangles.
    edges
}

fn is_sharp(edge: &[usize; 2], faces: &Vec<[usize; 3]>, verts: &Vec<Vector3<f64>>) -> bool {
    if faces.len() < 2 {
        return true;
    }
    // get edge direction
    let edge_dir = (verts[edge[1]] - verts[edge[0]]).normalize();
    let mid_edge = (verts[edge[1]] + verts[edge[0]]) / 2.0;

    // for all faces, get direction around edge
    let mut vec_dirs = Vec::new();
    for face in faces.iter() {
        let opp_vert = if face[0] != edge[0] && face[0] != edge[1] {
            verts[face[0]]
        } else if face[1] != edge[0] && face[1] != edge[1] {
            verts[face[1]]
        } else {
            verts[face[2]]
        };
        let mut dir = opp_vert - mid_edge;
        dir = (dir - dir.dot(&edge_dir) * edge_dir).normalize();
        vec_dirs.push(dir);
    }

    let base_vec0 = vec_dirs[0];
    let base_vec1 = edge_dir.cross(&base_vec0).normalize();

    let mut angles = Vec::new();
    for vec in vec_dirs.iter() {
        let x = base_vec0.dot(vec);
        let y = base_vec1.dot(vec);
        let angle = y.atan2(x);

        angles.push(angles);
    }

    angles.sort();

    // check if two consecutives angles are higher than 3pi/2
    for i in 0..angles.len() {
        let angle_cur = angles[i];
        let angle_next = if i == angles.len() - 1 {
            angles[0] + 2.0 * f64::consts::PI
        } else {
            angles[i + 1]
        };
        if angle_next - angle_cur > 3.0 * f64::consts::PI / 2.0 {
            return true;
        }
    }

    false
}

fn prune_sharp_edges(triangles: &Vec<[usize; 3]>, verts: &Vec<Vector3<f64>>) -> Vec<[usize; 3]> {
    let mut pruned_tri: Vec<[usize; 3]> = triangles.iter().map(|&f| f).collect();

    let mut nb_iter = 0;
    while nb_iter < 100 {
        // get edges from triangles
        let mut edges = get_edges_from_triangles(&pruned_tri);

        // for each edge, checks if it is sharp, and get the corresponding faces if so
        let mut rem_faces = Vec::new();
        let mut faces_removed = false;
        for (edge, faces) in edges.iter_mut() {
            if is_sharp(edge, faces, verts) {
                rem_faces.append(faces);
                faces_removed = true;
            }
        }
        if !faces_removed {
            return pruned_tri;
        }

        pruned_tri = pruned_tri
            .iter()
            .filter_map(|&f| {
                if rem_faces.contains(&f) {
                    None
                } else {
                    Some(f)
                }
            })
            .collect();
        nb_iter = nb_iter + 1;
    }
    pruned_tri
}

pub fn extract_manifold(triangles: &Vec<[usize; 3]>, verts: &Vec<Vector3<f64>>) {
    let pruned_tri = prune_sharp_edges(triangles, verts);
}
