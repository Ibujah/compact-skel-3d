use anyhow::Result;
use nalgebra::base::*;
use rand::Rng;
use std::collections::HashMap;

use crate::algorithm::skeleton_interface::{
    skeleton_path::PathPart, skeleton_separation::SkeletonSeparation, SkeletonInterface3D,
};
use crate::mesh3d::GenericMesh3D;

pub fn first_node_in(skeleton_interface: &mut SkeletonInterface3D) -> Result<usize> {
    let mut rng = rand::thread_rng();
    let rand_fac = rng.gen_range(0..skeleton_interface.mesh.get_nb_faces());
    println!("First face: {}", rand_fac);

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
            println!(
                "First node : ({}, {}, {}, {})",
                node.delaunay_tetrahedron()[0],
                node.delaunay_tetrahedron()[1],
                node.delaunay_tetrahedron()[2],
                node.delaunay_tetrahedron()[3]
            );
            return Ok(node.ind());
        }
    }

    Err(anyhow::Error::msg("No first node found"))
}

pub fn first_alveola_in(skeleton_interface: &mut SkeletonInterface3D) -> Result<usize> {
    let ind_first_node = first_node_in(skeleton_interface)?;

    let cur_node = skeleton_interface.get_node(ind_first_node)?;

    let edges = cur_node.edges();

    for edge in edges {
        let tri = edge.delaunay_triangle();

        if skeleton_interface
            .get_mesh()
            .is_face_in(tri[0], tri[1], tri[2])
            .is_none()
        {
            for alve in edge.alveolae() {
                let seg = alve.delaunay_segment();
                if skeleton_interface
                    .get_mesh()
                    .is_edge_in(seg[0], seg[1])
                    .is_none()
                {
                    return Ok(alve.ind());
                }
            }
        }
    }

    Err(anyhow::Error::msg("No first alveola found"))
}

pub fn include_alveola_in_skel(
    skeleton_interface: &mut SkeletonInterface3D,
    ind_alveola: usize,
    opt_label: Option<usize>,
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
            .partial_edge_next()
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
    if let Some(label) = opt_label {
        skeleton_interface.skeleton.set_label(ind_alveola, label);
    }
    Ok(())
}

pub fn neighbor_alveolae(
    skeleton_interface: &mut SkeletonInterface3D,
    ind_alveola: usize,
) -> Result<Vec<usize>> {
    let mut vec_neigh = Vec::new();

    let alve = skeleton_interface.get_alveola(ind_alveola)?;

    for edge in alve.edges() {
        if !edge.is_computed() {
            return Err(anyhow::Error::msg("Alveola not fully computed"));
        }
        for alv in edge.alveolae() {
            if alv.ind() != ind_alveola {
                vec_neigh.push(alv.ind());
            }
        }
    }

    Ok(vec_neigh)
}

pub fn compute_sheet(
    skeleton_interface: &mut SkeletonInterface3D,
    ind_alveola: usize,
    label: usize,
) -> Result<()> {
    skeleton_interface.get_alveola(ind_alveola)?;
    let mut to_compute = Vec::new();
    to_compute.push(ind_alveola);

    loop {
        if let Some(ind_alveola) = to_compute.pop() {
            let alve = skeleton_interface.get_alveola_uncheck(ind_alveola);
            if alve.label() != Some(label) {
                if !alve.is_computed() {
                    skeleton_interface.compute_alveola(ind_alveola)?;
                }

                for edge in skeleton_interface
                    .get_alveola_uncheck(ind_alveola)
                    .edges()
                    .iter()
                    .filter(|edge| edge.degree() == 2)
                {
                    edge.alveolae().iter().fold((), |_, alv| {
                        if alv.ind() != ind_alveola && alv.is_full() {
                            to_compute.push(alv.ind());
                        }
                    });
                }
                skeleton_interface.alve_label[ind_alveola] = Some(label);
            }
        } else {
            break;
        }
    }

    Ok(())
}

pub fn outer_partial_edges(
    skeleton_interface: &SkeletonInterface3D,
    current_sheet: &Vec<usize>,
) -> Vec<usize> {
    let mut vec_pedges = Vec::new();
    for &ind_alveola in current_sheet.iter() {
        for palve in skeleton_interface
            .get_alveola_uncheck(ind_alveola)
            .partial_alveolae()
        {
            for pedge in palve.partial_edges().iter() {
                let pedge_neigh = pedge.partial_edge_neighbor();
                if pedge_neigh.edge().is_singular() {
                    vec_pedges.push(pedge_neigh.ind());
                }
            }
        }
    }
    vec_pedges
}

pub fn extract_skeleton_path<'a, 'b>(
    skeleton_interface: &'b mut SkeletonInterface3D<'a>,
    ind_pedge: usize,
) -> Result<Option<SkeletonSeparation<'a, 'b>>> {
    let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge);
    if pedge.is_singular() {
        let mut skeleton_separation = SkeletonSeparation::new(skeleton_interface, ind_pedge);
        skeleton_separation.follow_separation()?;
        return Ok(Some(skeleton_separation));
    }

    Ok(None)
}

pub fn try_remove_and_add<'a, 'b>(
    skeleton_interface: &'b mut SkeletonInterface3D<'a>,
    vec_rem_faces: &Vec<usize>,
    vec_add_faces: &Vec<[usize; 3]>,
) -> Result<bool> {
    let mut vec_fac = Vec::new();
    for &ind_face in vec_rem_faces {
        vec_fac.push(skeleton_interface.mesh.get_face_vertices(ind_face).unwrap());
        skeleton_interface.mesh.remove_face(ind_face)?;
    }
    let mut vec_added = Vec::new();
    for &[ind_v1, ind_v2, ind_v3] in vec_add_faces {
        let res = skeleton_interface.mesh.add_face(ind_v1, ind_v2, ind_v3);
        match res {
            Ok(o) => vec_added.push(o),
            Err(_) => {
                for &f in vec_added.iter() {
                    skeleton_interface.mesh.remove_face(f)?;
                }
                for &[v1, v2, v3] in vec_fac.iter() {
                    skeleton_interface.mesh.add_face(v1, v2, v3)?;
                }
                return Ok(false);
            }
        }
    }
    Ok(true)
}

pub fn create_debug_meshes<'a, 'b>(
    skeleton_separation: &SkeletonSeparation<'a, 'b>,
    vec_rem_faces: &Vec<usize>,
    vec_add_faces: &Vec<[usize; 3]>,
) -> Result<Vec<GenericMesh3D>> {
    let mut debug_meshes = Vec::new();
    let mut debug_rem = GenericMesh3D::new();
    for &ind_face in vec_rem_faces {
        let [ind_v1, ind_v2, ind_v3] = skeleton_separation
            .skeleton_interface()
            .get_mesh()
            .get_face_vertices(ind_face)
            .unwrap();
        let pt1 = skeleton_separation
            .skeleton_interface()
            .get_mesh()
            .get_vertex(ind_v1)?
            .vertex();
        let pt2 = skeleton_separation
            .skeleton_interface()
            .get_mesh()
            .get_vertex(ind_v2)?
            .vertex();
        let pt3 = skeleton_separation
            .skeleton_interface()
            .get_mesh()
            .get_vertex(ind_v3)?
            .vertex();
        let i1 = debug_rem.add_vertex(&pt1);
        let i2 = debug_rem.add_vertex(&pt2);
        let i3 = debug_rem.add_vertex(&pt3);
        debug_rem.add_face(i1, i2, i3)?;
    }
    debug_meshes.push(debug_rem);
    let mut debug_add = GenericMesh3D::new();
    for &[ind_v1, ind_v2, ind_v3] in vec_add_faces {
        let pt1 = skeleton_separation
            .skeleton_interface()
            .get_mesh()
            .get_vertex(ind_v1)?
            .vertex();
        let pt2 = skeleton_separation
            .skeleton_interface()
            .get_mesh()
            .get_vertex(ind_v2)?
            .vertex();
        let pt3 = skeleton_separation
            .skeleton_interface()
            .get_mesh()
            .get_vertex(ind_v3)?
            .vertex();
        let i1 = debug_add.add_vertex(&pt1);
        let i2 = debug_add.add_vertex(&pt2);
        let i3 = debug_add.add_vertex(&pt3);
        debug_add.add_face(i1, i2, i3)?;
    }
    debug_meshes.push(debug_add);

    let mut debug_path_ext = GenericMesh3D::new();
    for ind1 in 0..skeleton_separation.external_path().components().len() {
        let ind2 = (ind1 + 1) % skeleton_separation.external_path().components().len();
        match (
            skeleton_separation.external_path().components()[ind1],
            skeleton_separation.external_path().components()[ind2],
        ) {
            (PathPart::PartialNode(ind_pnode1), PathPart::PartialNode(ind_pnode2)) => {
                let corner1 = skeleton_separation
                    .skeleton_interface()
                    .get_partial_node_uncheck(ind_pnode1)
                    .corner();
                let corner2 = skeleton_separation
                    .skeleton_interface()
                    .get_partial_node_uncheck(ind_pnode2)
                    .corner();
                if corner1 != corner2 {
                    let (center, _) = skeleton_separation
                        .skeleton_interface()
                        .get_partial_node_uncheck(ind_pnode1)
                        .node()
                        .center_and_radius()?;
                    let pt1 = skeleton_separation
                        .skeleton_interface()
                        .get_mesh()
                        .get_vertex(corner1)?
                        .vertex();
                    let pt2 = skeleton_separation
                        .skeleton_interface()
                        .get_mesh()
                        .get_vertex(corner2)?
                        .vertex();
                    let i1 = debug_path_ext.add_vertex(&pt1);
                    let i2 = debug_path_ext.add_vertex(&pt2);
                    let i3 = debug_path_ext.add_vertex(&center);
                    debug_path_ext.add_face(i1, i2, i3)?;
                }
            }
            (PathPart::PartialEdge(ind_pedge), _) => {
                let pedge = skeleton_separation
                    .skeleton_interface()
                    .get_partial_edge_uncheck(ind_pedge);
                let corner = skeleton_separation
                    .skeleton_interface()
                    .get_mesh()
                    .get_vertex(pedge.corner())?
                    .vertex();
                let (ctr1, _) = pedge
                    .partial_node_first()
                    .unwrap()
                    .node()
                    .center_and_radius()?;
                let (ctr2, _) = pedge
                    .partial_node_last()
                    .unwrap()
                    .node()
                    .center_and_radius()?;
                let i1 = debug_path_ext.add_vertex(&ctr1);
                let i2 = debug_path_ext.add_vertex(&ctr2);
                let i3 = debug_path_ext.add_vertex(&corner);
                debug_path_ext.add_face(i1, i2, i3)?;
            }
            (_, _) => (),
        }
    }
    debug_meshes.push(debug_path_ext);

    for path_int in skeleton_separation.internal_paths() {
        let mut debug_path_int = GenericMesh3D::new();
        for ind1 in 0..path_int.components().len() {
            let ind2 = (ind1 + 1) % path_int.components().len();
            match (path_int.components()[ind1], path_int.components()[ind2]) {
                (PathPart::PartialNode(ind_pnode1), PathPart::PartialNode(ind_pnode2)) => {
                    let corner1 = skeleton_separation
                        .skeleton_interface()
                        .get_partial_node_uncheck(ind_pnode1)
                        .corner();
                    let corner2 = skeleton_separation
                        .skeleton_interface()
                        .get_partial_node_uncheck(ind_pnode2)
                        .corner();
                    if corner1 != corner2 {
                        let (center, _) = skeleton_separation
                            .skeleton_interface()
                            .get_partial_node_uncheck(ind_pnode1)
                            .node()
                            .center_and_radius()?;
                        let pt1 = skeleton_separation
                            .skeleton_interface()
                            .get_mesh()
                            .get_vertex(corner1)?
                            .vertex();
                        let pt2 = skeleton_separation
                            .skeleton_interface()
                            .get_mesh()
                            .get_vertex(corner2)?
                            .vertex();
                        let i1 = debug_path_int.add_vertex(&pt1);
                        let i2 = debug_path_int.add_vertex(&pt2);
                        let i3 = debug_path_int.add_vertex(&center);
                        debug_path_int.add_face(i1, i2, i3)?;
                    }
                }
                (PathPart::PartialEdge(ind_pedge), _) => {
                    let pedge = skeleton_separation
                        .skeleton_interface()
                        .get_partial_edge_uncheck(ind_pedge);
                    let corner = skeleton_separation
                        .skeleton_interface()
                        .get_mesh()
                        .get_vertex(pedge.corner())?
                        .vertex();
                    let (ctr1, _) = pedge
                        .partial_node_first()
                        .unwrap()
                        .node()
                        .center_and_radius()?;
                    let (ctr2, _) = pedge
                        .partial_node_last()
                        .unwrap()
                        .node()
                        .center_and_radius()?;
                    let i1 = debug_path_int.add_vertex(&ctr1);
                    let i2 = debug_path_int.add_vertex(&ctr2);
                    let i3 = debug_path_int.add_vertex(&corner);
                    debug_path_int.add_face(i1, i2, i3)?;
                }
                (_, _) => (),
            }
        }
        debug_meshes.push(debug_path_int);
    }

    Ok(debug_meshes)
}
