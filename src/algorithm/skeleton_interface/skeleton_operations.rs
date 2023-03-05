use anyhow::Result;
use nalgebra::base::*;
use rand::Rng;
use std::collections::HashMap;

use crate::algorithm::skeleton_interface::SkeletonInterface3D;

pub fn first_node_in<'a, 'b, 'c>(
    skeleton_interface: &'c mut SkeletonInterface3D<'a, 'b>,
) -> Result<usize> {
    let mut rng = rand::thread_rng();
    let rand_fac = rng.gen_range(0..skeleton_interface.mesh.get_nb_faces());

    let face = skeleton_interface.get_mesh().get_face(rand_fac)?;

    let mut triangle = face.vertices_inds();
    triangle.sort();
    let hedges = face.halfedges();
    let vec1 = hedges[0].last_vertex().vertex() - hedges[0].first_vertex().vertex();
    let vec2 = hedges[1].last_vertex().vertex() - hedges[1].first_vertex().vertex();

    let normal = vec1.cross(&vec2).normalize();
    let pt_face = hedges[0].first_vertex().vertex();

    let vec_tets = skeleton_interface
        .faces
        .get(&triangle)
        .ok_or(anyhow::Error::msg("No tetrahedron surrounding triangle"))?;

    for i in 0..vec_tets.len() {
        let tet = vec_tets[i];
        let v0 = skeleton_interface.get_mesh().get_vertex(tet[0])?.vertex();
        let v1 = skeleton_interface.get_mesh().get_vertex(tet[1])?.vertex();
        let v2 = skeleton_interface.get_mesh().get_vertex(tet[2])?.vertex();
        let v3 = skeleton_interface.get_mesh().get_vertex(tet[3])?.vertex();

        let v_mean = (v0 + v1 + v2 + v3) * 0.25;

        let inside = normal.dot(&(v_mean - pt_face)) < 0.0;

        if inside {
            let node = skeleton_interface.add_node(&tet)?;
            return Ok(node.ind());
        }
    }

    Err(anyhow::Error::msg("No first node found"))
}

pub fn propagate_edge(skeleton_interface: &mut SkeletonInterface3D, ind_edge: usize) -> Result<()> {
    let del_tri = skeleton_interface.edge_tri[ind_edge];
    let del_tets = skeleton_interface.get_tetrahedra_from_triangle(del_tri)?;
    for del_tet in del_tets {
        skeleton_interface.add_node(&del_tet)?;
    }
    Ok(())
}

pub fn compute_alveola(
    skeleton_interface: &mut SkeletonInterface3D,
    ind_alveola: usize,
) -> Result<()> {
    let vec_pedg = skeleton_interface
        .get_alveola(ind_alveola)?
        .partial_alveolae()[0]
        .partial_edges();
    let pedg_first = vec_pedg[0].ind();
    let mut pedg_cur = pedg_first;
    loop {
        propagate_edge(skeleton_interface, skeleton_interface.pedge_edge[pedg_cur])?;
        pedg_cur = skeleton_interface
            .get_partial_edge(pedg_cur)?
            .partial_edge_next()?
            .ok_or(anyhow::Error::msg("Non complete partial edge"))?
            .ind();
        if pedg_cur == pedg_first {
            break;
        }
    }
    Ok(())
}

pub fn include_alveola_in_skel(
    skeleton_interface: &mut SkeletonInterface3D,
    ind_alveola: usize,
) -> Result<()> {
    let vec_pedg = skeleton_interface
        .get_alveola(ind_alveola)?
        .partial_alveolae()[0]
        .partial_edges();

    let mut bnd_pts = HashMap::new();
    let mut edges_map = HashMap::new();
    let mut lis_nods = Vec::new();
    let mut pedg_cur = vec_pedg[0];
    for _i in 0..vec_pedg.len() {
        let node = pedg_cur
            .partial_node_first()
            .ok_or(anyhow::Error::msg("Uncomputed node"))?
            .node();
        let boundary_points = node
            .delaunay_tetrahedron()
            .iter()
            .map(|&ind_vertex| Ok(skeleton_interface.mesh.get_vertex(ind_vertex)?.vertex()))
            .collect::<Result<Vec<Vector3<f32>>>>()?
            .try_into()
            .map_err(|_x: Vec<_>| anyhow::Error::msg("Could not convert vec to array"))
            .unwrap();
        bnd_pts.insert(node.ind(), boundary_points);
        let ind_edg_cur = pedg_cur.edge().ind();
        let nod_ext = [
            pedg_cur.edge().nodes()[0].ind(),
            pedg_cur.edge().nodes()[1].ind(),
        ];
        edges_map.insert(ind_edg_cur, nod_ext);
        lis_nods.push(node.ind());
        pedg_cur = pedg_cur
            .partial_edge_next()?
            .ok_or(anyhow::Error::msg("Non complete partial edge"))?;
    }
    for (ind_nod, boundary_points) in bnd_pts {
        skeleton_interface
            .skeleton
            .add_node(ind_nod, boundary_points)?;
    }
    for (ind_edge, ind_nodes) in edges_map {
        skeleton_interface.skeleton.add_edge(ind_edge, ind_nodes);
    }
    skeleton_interface
        .skeleton
        .add_alveola(ind_alveola, lis_nods);
    Ok(())
}

pub fn neighbor_alveolae(
    skeleton_interface: &mut SkeletonInterface3D,
    ind_alveola: usize,
) -> Result<Vec<usize>> {
    let mut vec_neigh = Vec::new();

    let alve = skeleton_interface.get_alveola(ind_alveola)?;

    for edge in alve.edges() {
        for alv in edge.alveolae() {
            if alv.ind() != ind_alveola {
                vec_neigh.push(alv.ind());
            }
        }
    }

    Ok(vec_neigh)
}
