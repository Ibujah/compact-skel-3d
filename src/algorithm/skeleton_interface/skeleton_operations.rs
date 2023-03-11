use anyhow::Result;
use nalgebra::base::*;
use rand::Rng;
use std::collections::{HashMap, HashSet};

use crate::algorithm::skeleton_interface::{skeleton_path::SkeletonPath, SkeletonInterface3D};

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

pub fn first_alveola_in<'a, 'b, 'c>(
    skeleton_interface: &'c mut SkeletonInterface3D<'a, 'b>,
) -> Result<usize> {
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

fn follow_singular_path(skeleton_path: &mut SkeletonPath) -> Result<()> {
    loop {
        if let Some(ind_pedge) = skeleton_path.ind_last_partial_edge() {
            // print!("{}: last {} ", cpt, ind_pedge);
            let pedge = skeleton_path
                .skeleton_interface()
                .get_partial_edge_uncheck(ind_pedge);
            // print!(
            //     "{} -> {} (deg {}) ",
            //     pedge.partial_node_first().unwrap().node().ind(),
            //     pedge.partial_node_last().unwrap().node().ind(),
            //     pedge.edge().degree()
            // );
            if pedge.is_singular() {
                skeleton_path.append_last()?;
                // println!("append ");
            } else {
                skeleton_path.rotate_last()?;
                // println!("rotate ");
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
) -> HashSet<usize> {
    let mut pedges_set = HashSet::new();
    for &ind_alveola in current_sheet.iter() {
        for palve in skeleton_interface
            .get_alveola_uncheck(ind_alveola)
            .partial_alveolae()
        {
            for pedge in palve.partial_edges().iter() {
                let pedge_neigh = pedge.partial_edge_neighbor();
                if pedge_neigh.is_singular() {
                    pedges_set.insert(pedge_neigh.ind());
                }
            }
        }
    }
    pedges_set
}

pub fn extract_one_skeleton_path<'a, 'b, 'c>(
    skeleton_interface: &'c mut SkeletonInterface3D<'a, 'b>,
    pedges_set: &HashSet<usize>,
) -> Result<Option<SkeletonPath<'a, 'b, 'c>>> {
    for ind_pedge in pedges_set.iter() {
        let pedge = skeleton_interface.get_partial_edge_uncheck(*ind_pedge);
        if pedge.is_singular() {
            let mut skeleton_path = SkeletonPath::new(skeleton_interface, *ind_pedge);
            follow_singular_path(&mut skeleton_path)?;
            return Ok(Some(skeleton_path));
        }
    }

    Ok(None)
}
