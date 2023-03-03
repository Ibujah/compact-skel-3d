use anyhow::Result;
use nalgebra::base::*;
use std::collections::HashMap;

use crate::algorithm::skeleton_interface::SkeletonInterface3D;

pub fn propagate_edge(voro: &mut SkeletonInterface3D, ind_edge: usize) -> Result<()> {
    let del_tri = voro.edge_tri[ind_edge];
    let del_tets = voro.get_tetrahedra_from_triangle(del_tri)?;
    for del_tet in del_tets {
        voro.add_node(&del_tet)?;
    }
    Ok(())
}

pub fn compute_alveola(voro: &mut SkeletonInterface3D, ind_alveola: usize) -> Result<()> {
    let vec_pedg = voro.get_alveola(ind_alveola).partial_alveolae()[0].partial_edges();
    let pedg_first = vec_pedg[0].ind();
    let mut pedg_cur = pedg_first;
    loop {
        propagate_edge(voro, voro.pedge_edge[pedg_cur])?;
        pedg_cur = voro
            .get_partial_edge(pedg_cur)
            .partial_edge_next()?
            .ok_or(anyhow::Error::msg("Non complete partial edge"))?
            .ind();
        if pedg_cur == pedg_first {
            break;
        }
    }
    Ok(())
}

pub fn include_alveola_in_skel(voro: &mut SkeletonInterface3D, ind_alveola: usize) -> Result<()> {
    let vec_pedg = voro.get_alveola(ind_alveola).partial_alveolae()[0].partial_edges();

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
            .map(|&ind_vertex| Ok(voro.mesh.get_vertex(ind_vertex)?.vertex()))
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
        voro.skeleton.add_node(ind_nod, boundary_points)?;
    }
    for (ind_edge, ind_nodes) in edges_map {
        voro.skeleton.add_edge(ind_edge, ind_nodes);
    }
    voro.skeleton.add_alveola(ind_alveola, lis_nods);
    Ok(())
}
