use anyhow::Result;
use nalgebra::base::*;
use rand::Rng;
use std::collections::{HashMap, HashSet};

use crate::mesh3d::GenericMesh3D;

use super::path_part::PathPart;
use super::MovableDelaunayPath;
use super::SkeletonInterface3D;
use super::SkeletonSeparation;

/// Computes a random first node on skeleton
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

/// Computes a random first alveola on skeleton
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

/// Includes an alveola in the final skeleton
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

/// Returns neighbor alveola indices
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

/// Computes full sheet starting from an alveola
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

/// Returns neighbor partial edges to each singular edge on the sheet
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

/// Computes skeleton separation starting from a partial edge
pub fn extract_skeleton_separation<'a, 'b>(
    skeleton_interface: &'b mut SkeletonInterface3D<'a>,
    ind_pedge: usize,
) -> Result<Option<SkeletonSeparation<'a, 'b>>> {
    let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge);
    if pedge.is_singular() {
        let mut skeleton_separation = SkeletonSeparation::create(skeleton_interface, ind_pedge)?;
        skeleton_separation.follow_separation()?;
        return Ok(Some(skeleton_separation));
    }

    Ok(None)
}

/// Tries to remove a set of faces and add a set of face to the mesh
///
/// If operation fails, leaves mesh unchanged
pub fn try_remove_and_add<'a, 'b>(
    skeleton_interface: &'b mut SkeletonInterface3D<'a>,
    vec_rem_faces: &Vec<usize>,
    vec_add_faces: &Vec<[usize; 3]>,
) -> Result<bool> {
    let mut vec_fac = Vec::new();
    for &ind_face in vec_rem_faces {
        vec_fac.push(
            skeleton_interface
                .mesh
                .get_face(ind_face)
                .unwrap()
                .vertices_inds(),
        );
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

/// Collect list of faces on mesh portion described by separation
pub fn collect_mesh_faces_index(
    skeleton_separation: &SkeletonSeparation,
    epsilon: f32,
) -> Result<Option<Vec<usize>>> {
    fn last_hedge_deletion(
        mesh_paths_external: &mut Vec<Vec<usize>>,
        skeleton_separation: &SkeletonSeparation,
    ) -> Result<bool> {
        if let Some(mut mesh_path_external) = mesh_paths_external.pop() {
            if let Some(ind_hedge) = mesh_path_external.pop() {
                let hedge = skeleton_separation
                    .skeleton_interface()
                    .get_mesh()
                    .get_halfedge(ind_hedge)?;
                let ind_hedge_opp = hedge.opposite_halfedge().unwrap().ind();
                let opt_position = mesh_path_external
                    .iter()
                    .position(|&ind| ind == ind_hedge_opp);
                if let Some(position) = opt_position {
                    if position != 0 {
                        let mut path1 = vec![0 as usize; position];
                        path1.copy_from_slice(&mesh_path_external[..position]);
                        mesh_paths_external.push(path1);
                    }
                    if position != mesh_path_external.len() - 1 {
                        let mut path2 = vec![0 as usize; mesh_path_external.len() - 1 - position];
                        path2.copy_from_slice(&mesh_path_external[position + 1..]);
                        mesh_paths_external.push(path2);
                    }
                    return Ok(true);
                }
                mesh_path_external.push(ind_hedge);
            }
            mesh_paths_external.push(mesh_path_external);
        }
        Ok(false)
    }

    fn last_hedge_fusion(
        mesh_paths_external: &mut Vec<Vec<usize>>,
        mesh_paths_internal: &mut Vec<Vec<usize>>,
        skeleton_separation: &SkeletonSeparation,
    ) -> Result<bool> {
        if let Some(mut mesh_path_external) = mesh_paths_external.pop() {
            if let Some(ind_hedge) = mesh_path_external.pop() {
                let hedge = skeleton_separation
                    .skeleton_interface()
                    .get_mesh()
                    .get_halfedge(ind_hedge)?;
                let ind_hedge_opp = hedge.opposite_halfedge().unwrap().ind();
                let mut ind_pa_he = None;
                for ind_pa in 0..mesh_paths_internal.len() {
                    for ind_he in 0..mesh_paths_internal[ind_pa].len() {
                        if mesh_paths_internal[ind_pa][ind_he] == ind_hedge_opp {
                            ind_pa_he = Some((ind_pa, ind_he));
                        }
                    }
                }
                if let Some((ind_pa, ind_he)) = ind_pa_he {
                    let mesh_path = mesh_paths_external.remove(ind_pa);
                    for i in (ind_he + 1)..mesh_path.len() {
                        mesh_path_external.push(mesh_path[i]);
                    }
                    for i in 0..ind_he {
                        mesh_path_external.push(mesh_path[i]);
                    }
                    mesh_paths_external.push(mesh_path_external.clone());
                    return Ok(true);
                }
                mesh_path_external.push(ind_hedge);
            }
            mesh_paths_external.push(mesh_path_external);
        }
        Ok(false)
    }

    fn last_hedge_expansion(
        mesh_paths_external: &mut Vec<Vec<usize>>,
        skeleton_separation: &SkeletonSeparation,
        center_mat: &MatrixXx3<f32>,
        radius_mat: &MatrixXx1<f32>,
        epsilon: f32,
        faces: &mut Vec<usize>,
    ) -> Result<bool> {
        if let Some(mut mesh_path_external) = mesh_paths_external.pop() {
            if let Some(ind_hedge) = mesh_path_external.pop() {
                let hedge = skeleton_separation
                    .skeleton_interface()
                    .get_mesh()
                    .get_halfedge(ind_hedge)?;
                let vert_test = hedge
                    .next_halfedge()
                    .unwrap()
                    .last_vertex()
                    .vertex()
                    .transpose();

                if center_mat
                    .row_iter()
                    .zip(radius_mat.iter())
                    .find(|(row, &rad)| {
                        let diff = row - vert_test;
                        diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]
                            < rad * rad + 2.0 * rad * epsilon + epsilon * epsilon
                    })
                    .is_none()
                {
                    return Ok(false);
                }

                let face = hedge.face().unwrap();

                faces.push(face.ind());
                let hedge_rep1 = hedge.prev_halfedge().unwrap().opposite_halfedge().unwrap();
                let hedge_rep2 = hedge.next_halfedge().unwrap().opposite_halfedge().unwrap();
                mesh_path_external.push(hedge_rep1.ind());
                mesh_path_external.push(hedge_rep2.ind());
                mesh_paths_external.push(mesh_path_external.clone());
                return Ok(true);
            }
        }
        Err(anyhow::Error::msg("Paths should not be empty"))
    }

    let (center_mat, radius_mat) = skeleton_separation
        .external_path()
        .basis_spheres_matrices(&skeleton_separation.skeleton_interface())?;
    let mut mesh_paths_external = {
        let mesh_path_external = skeleton_separation
            .external_path()
            .halfedges_path(&skeleton_separation.skeleton_interface())?;
        vec![mesh_path_external]
    };
    let mut mesh_paths_internal = Vec::new();
    for skeleton_path_internal in skeleton_separation.internal_paths().iter() {
        let mesh_path_internal =
            skeleton_path_internal.halfedges_path(&skeleton_separation.skeleton_interface())?;
        mesh_paths_internal.push(mesh_path_internal);
    }

    let mut faces = Vec::new();
    loop {
        // println!("new iter");
        // for path in mesh_paths_hedge.iter() {
        //     for &ind_he in path.iter() {
        //         let hedge = skeleton_separation.skeleton_interface.get_mesh().get_halfedge(ind_he)?;
        //         print!(
        //             "({} -> {}), ",
        //             hedge.first_vertex().ind(),
        //             hedge.last_vertex().ind()
        //         );
        //     }
        //     println!("");
        // }
        // println!("");

        // emptying paths
        loop {
            if let Some(mesh_path_external) = mesh_paths_external.last() {
                if mesh_path_external.len() == 0 {
                    mesh_paths_external.pop();
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        if mesh_paths_external.is_empty() {
            break;
        }

        if last_hedge_deletion(&mut mesh_paths_external, &skeleton_separation)? {
            continue;
        }

        if last_hedge_fusion(
            &mut mesh_paths_external,
            &mut mesh_paths_internal,
            &skeleton_separation,
        )? {
            continue;
        }

        if !last_hedge_expansion(
            &mut mesh_paths_external,
            &skeleton_separation,
            &center_mat,
            &radius_mat,
            epsilon,
            &mut faces,
        )? {
            return Ok(None);
        }

        if faces.len()
            > skeleton_separation
                .skeleton_interface()
                .get_mesh()
                .get_nb_faces()
                * 2
        {
            return Ok(None);
        }
    }

    Ok(Some(faces))
}

/// Estimates Delaunay faces to add on mesh to close the given separation
pub fn collect_closing_faces(
    skeleton_separation: &SkeletonSeparation,
    removed_faces: &Vec<usize>,
) -> Result<Option<Vec<[usize; 3]>>> {
    let mut unfaced_hedges = HashSet::new();
    for &ind_fac in removed_faces {
        let fac = skeleton_separation
            .skeleton_interface()
            .get_mesh()
            .get_face(ind_fac)?;
        let [hedg0, hedg1, hedg2] = fac.halfedges();
        unfaced_hedges.insert(hedg0.halfedge());
        unfaced_hedges.insert(hedg1.halfedge());
        unfaced_hedges.insert(hedg2.halfedge());
    }

    let mut palve_paths_external = {
        let ind_palves = skeleton_separation
            .external_path()
            .alveolae_path(&skeleton_separation.skeleton_interface())?;
        vec![MovableDelaunayPath::create(
            skeleton_separation.skeleton_interface(),
            ind_palves,
            &unfaced_hedges,
        )?]
    };
    let mut palve_paths_internal = Vec::new();
    for internal_path in skeleton_separation.internal_paths().iter() {
        let ind_palves = internal_path.alveolae_path(&skeleton_separation.skeleton_interface())?;
        palve_paths_internal.push(MovableDelaunayPath::create(
            skeleton_separation.skeleton_interface(),
            ind_palves,
            &unfaced_hedges,
        )?);
    }

    let mut closing_faces = Vec::new();
    loop {
        // println!("new iter");
        // for path in palve_paths.iter() {
        //     for &ind_palve in path.iter() {
        //         let palve = skeleton_separation
        //             .skeleton_interface
        //             .get_partial_alveola_uncheck(ind_palve);
        //         let seg = palve.alveola().delaunay_segment();
        //         let seg = if seg[0] == palve.corner() {
        //             seg
        //         } else {
        //             [seg[1], seg[0]]
        //         };
        //         print!("({} -> {}), ", seg[0], seg[1]);
        //     }
        //     println!("");
        // }
        // println!("");

        if let Some(mut palve_path) = palve_paths_external.pop() {
            if palve_path.closed() {
                continue;
            }
            if closing_faces.len()
                > skeleton_separation
                    .skeleton_interface()
                    .get_mesh()
                    .get_nb_faces()
            {
                return Err(anyhow::Error::msg("Too much faces collected somehow"));
            }

            if let Some(couple) = palve_path.get_couple_to_fusion() {
                let mut paths_res = palve_path.fusion_couple(couple)?;
                palve_paths_external.append(&mut paths_res);
                continue;
            } else if let Some(ind_exp) = palve_path.get_ind_to_expand() {
                palve_path.expand_ind(ind_exp, &mut closing_faces)?;
                palve_paths_external.push(palve_path);
                continue;
            } else {
                return Ok(None);
            }
        } else {
            break;
        }
    }
    Ok(Some(closing_faces))
}

/// (Debug) Estimates Delaunay faces to add on mesh to close the given separation
pub fn collectable_closing_faces(
    skeleton_separation: &SkeletonSeparation,
    removed_faces: &Vec<usize>,
) -> Result<Vec<[usize; 3]>> {
    let mut unfaced_hedges = HashSet::new();
    for &ind_fac in removed_faces {
        let fac = skeleton_separation
            .skeleton_interface()
            .get_mesh()
            .get_face(ind_fac)?;
        let [hedg0, hedg1, hedg2] = fac.halfedges();
        unfaced_hedges.insert(hedg0.halfedge());
        unfaced_hedges.insert(hedg1.halfedge());
        unfaced_hedges.insert(hedg2.halfedge());
    }

    let mut palve_paths_external = {
        let ind_palves = skeleton_separation
            .external_path()
            .alveolae_path(&skeleton_separation.skeleton_interface())?;
        vec![MovableDelaunayPath::create(
            skeleton_separation.skeleton_interface(),
            ind_palves,
            &unfaced_hedges,
        )?]
    };
    let mut palve_paths_internal = Vec::new();
    for internal_path in skeleton_separation.internal_paths().iter() {
        let ind_palves = internal_path.alveolae_path(&skeleton_separation.skeleton_interface())?;
        palve_paths_internal.push(MovableDelaunayPath::create(
            skeleton_separation.skeleton_interface(),
            ind_palves,
            &unfaced_hedges,
        )?);
    }

    let mut closing_faces = Vec::new();
    loop {
        println!("new iter");
        for path in palve_paths_external.iter() {
            path.print();
            println!("");
        }
        println!("");

        if let Some(mut palve_path) = palve_paths_external.pop() {
            if palve_path.closed() {
                continue;
            }
            if closing_faces.len()
                > skeleton_separation
                    .skeleton_interface()
                    .get_mesh()
                    .get_nb_faces()
            {
                return Ok(closing_faces);
            }

            if let Some(couple) = palve_path.get_couple_to_fusion() {
                println!("fusion");
                let mut paths_res = palve_path.fusion_couple(couple)?;
                palve_paths_external.append(&mut paths_res);
                continue;
            } else if let Some(ind_exp) = palve_path.get_ind_to_expand() {
                println!("expand");
                palve_path.expand_ind(ind_exp, &mut closing_faces)?;
                palve_paths_external.push(palve_path);
                continue;
            } else {
                return Ok(closing_faces);
            }
        } else {
            break;
        }
    }
    Ok(closing_faces)
}

/// (Debug) Creates debug meshes associated to set of faces
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
            .get_face(ind_face)
            .unwrap()
            .vertices_inds();
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
