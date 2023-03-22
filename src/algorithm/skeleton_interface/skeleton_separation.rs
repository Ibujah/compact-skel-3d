use anyhow::Result;
use nalgebra::base::*;
use std::collections::HashSet;

use super::skeleton_path::SkeletonPath;
use crate::algorithm::skeleton_interface::SkeletonInterface3D;
use crate::mesh3d::manifold_mesh3d::IterHalfEdge;

pub struct SkeletonSeparation<'a, 'b> {
    skeleton_interface: &'b mut SkeletonInterface3D<'a>,
    external_path: SkeletonPath,
    internal_paths: Vec<SkeletonPath>,
}

impl<'a, 'b> SkeletonSeparation<'a, 'b> {
    pub fn new(
        skeleton_interface: &'b mut SkeletonInterface3D<'a>,
        ind_pedge: usize,
    ) -> SkeletonSeparation<'a, 'b> {
        SkeletonSeparation {
            skeleton_interface,
            external_path: SkeletonPath::new(ind_pedge),
            internal_paths: Vec::new(),
        }
    }

    pub fn skeleton_interface(&self) -> &SkeletonInterface3D {
        self.skeleton_interface
    }

    pub fn external_path(&self) -> &SkeletonPath {
        &self.external_path
    }

    pub fn internal_paths(&self) -> &Vec<SkeletonPath> {
        &self.internal_paths
    }

    fn follow_external_path(&mut self) -> Result<()> {
        self.external_path
            .follow_singular_path(&mut self.skeleton_interface)
    }

    fn internal_partial_edges(&self) -> Vec<usize> {
        let vec_pedges_ext = self.external_path.ind_partial_edges();
        let mut vec_pedges_int = Vec::new();
        for &ind_pedge in vec_pedges_ext.iter() {
            let ind_pedge_opp = self
                .skeleton_interface
                .get_partial_edge_uncheck(ind_pedge)
                .partial_edge_opposite()
                .ind();
            if vec_pedges_ext
                .iter()
                .find(|&&ind| ind == ind_pedge_opp)
                .is_none()
            {
                vec_pedges_int.push(ind_pedge_opp);
            }
        }
        vec_pedges_int.sort();
        vec_pedges_int.dedup();
        vec_pedges_int
    }

    fn follow_internal_paths(&mut self) -> Result<()> {
        let mut vec_internal_pedges = self.internal_partial_edges();
        loop {
            if let Some(ind_pedge) = vec_internal_pedges.pop() {
                let mut skeleton_path_int = SkeletonPath::new(ind_pedge);
                skeleton_path_int.follow_singular_path(&mut self.skeleton_interface)?;
                for &ind_pedge_new in skeleton_path_int.ind_partial_edges().iter() {
                    if let Some(pos) = vec_internal_pedges
                        .iter()
                        .position(|&ind| ind == ind_pedge_new)
                    {
                        vec_internal_pedges.remove(pos);
                    }
                }
                self.internal_paths.push(skeleton_path_int);
            } else {
                break;
            }
        }
        Ok(())
    }

    pub fn follow_separation(&mut self) -> Result<()> {
        self.follow_external_path()?;
        self.follow_internal_paths()?;
        Ok(())
    }

    pub fn closable_path(&self) -> Result<bool> {
        self.external_path.closable_path(&self.skeleton_interface)
        // Ok(true)
    }

    pub fn collect_closing_faces(
        &self,
        removed_faces: &Vec<usize>,
    ) -> Result<Option<Vec<[usize; 3]>>> {
        let mut palve_paths = {
            let palve_path = self.external_path.alveolae_path(&self.skeleton_interface)?;
            vec![palve_path]
        };
        let mut fixed_internal = HashSet::new();
        for internal_path in self.internal_paths.iter() {
            let palve_internal = internal_path.alveolae_path(&self.skeleton_interface)?;
            for ind_palve in palve_internal {
                fixed_internal.insert(ind_palve);
            }
        }
        let mut unfaced_hedges = HashSet::new();
        for &ind_fac in removed_faces {
            let fac = self.skeleton_interface.get_mesh().get_face(ind_fac)?;
            let [hedg0, hedg1, hedg2] = fac.halfedges();
            unfaced_hedges.insert(hedg0.halfedge());
            unfaced_hedges.insert(hedg1.halfedge());
            unfaced_hedges.insert(hedg2.halfedge());
        }
        let mut faces = Vec::new();
        loop {
            // println!("new iter");
            // for path in palve_paths.iter() {
            //     for &ind_palve in path.iter() {
            //         let palve = self
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

            if let Some(palve_path) = palve_paths.pop() {
                if palve_path.is_empty() {
                    continue;
                }
                let nb_palve = palve_path.len();
                let mut operation_done = false;
                let mut all_fixed = true;
                for ind in 0..nb_palve {
                    let ind_palve = palve_path[ind];
                    if fixed_internal.contains(&ind) {
                        continue;
                    }
                    all_fixed = false;

                    let palve = self
                        .skeleton_interface
                        .get_partial_alveola_uncheck(ind_palve);
                    let seg = palve.alveola().delaunay_segment();
                    let seg = if seg[0] == palve.corner() {
                        seg
                    } else {
                        [seg[1], seg[0]]
                    };
                    if unfaced_hedges.contains(&seg) {
                        continue;
                    }

                    let ind_palve_opp = palve.partial_alveola_opposite().ind();
                    let opt_position = palve_path
                        .iter()
                        .position(|&ind_alve| ind_alve == ind_palve_opp);
                    if let Some(position) = opt_position {
                        let mut path1 = Vec::new();
                        let mut path2 = Vec::new();
                        let (ind_min, ind_max) = if ind < position {
                            (ind, position)
                        } else {
                            (position, ind)
                        };
                        for i in 0..nb_palve {
                            if i < ind_min || i > ind_max {
                                path1.push(palve_path[i]);
                            }
                            if i > ind_min && i < ind_max {
                                path2.push(palve_path[i]);
                            }
                        }
                        palve_paths.push(path1);
                        palve_paths.push(path2);
                        operation_done = true;
                        break;
                    }
                }
                if all_fixed {
                    operation_done = true
                }
                if !operation_done {
                    for ind in 0..nb_palve {
                        let ind_palve = palve_path[ind];

                        let palve = self
                            .skeleton_interface
                            .get_partial_alveola_uncheck(ind_palve);
                        let ind_prev = (ind - 1 + nb_palve) % nb_palve;
                        let ind_palve_prev = palve_path[ind_prev];
                        let palve_prev = self
                            .skeleton_interface
                            .get_partial_alveola_uncheck(ind_palve_prev);
                        let ind_vertex = palve_prev.corner();
                        for pedge in palve.partial_edges() {
                            if !pedge.edge().is_full() {
                                continue;
                            }
                            let tri = pedge.edge().delaunay_triangle();
                            if tri.contains(&ind_vertex) {
                                let seg = palve.alveola().delaunay_segment();
                                if ind_vertex == seg[0] || ind_vertex == seg[1] {
                                    continue;
                                }
                                let seg = if seg[0] == palve.corner() {
                                    seg
                                } else {
                                    [seg[1], seg[0]]
                                };
                                let mut tri = [seg[0], seg[1], ind_vertex];
                                tri.sort();
                                faces.push([seg[0], seg[1], ind_vertex]);

                                let seg2 = if seg[1] < ind_vertex {
                                    [seg[1], ind_vertex]
                                } else {
                                    [ind_vertex, seg[1]]
                                };
                                let &ind_alveola2 =
                                    self.skeleton_interface.del_seg.get(&seg2).ok_or(
                                        anyhow::Error::msg(format!(
                                            "Partial alveola ({}, {}) not in skeleton",
                                            seg2[0], seg2[1]
                                        )),
                                    )?;
                                let alve2 =
                                    self.skeleton_interface.get_alveola_uncheck(ind_alveola2);
                                let ind_palve2 =
                                    if alve2.partial_alveolae()[0].corner() == ind_vertex {
                                        alve2.partial_alveolae()[0].ind()
                                    } else {
                                        alve2.partial_alveolae()[1].ind()
                                    };
                                let mut palve_path_new = Vec::new();
                                for i in 0..palve_path.len() {
                                    if i != ind && i != ind_prev {
                                        palve_path_new.push(palve_path[i]);
                                    }
                                    if i == ind {
                                        palve_path_new.push(ind_palve2);
                                    }
                                }
                                palve_paths.push(palve_path_new);

                                operation_done = true;
                                break;
                            }
                        }
                        if operation_done {
                            break;
                        }
                    }
                }

                if !operation_done {
                    return Ok(None);
                }
            } else {
                break;
            }
            if faces.len() > self.skeleton_interface.get_mesh().get_nb_faces() {
                return Err(anyhow::Error::msg("Too much faces collected somehow"));
            }
        }
        Ok(Some(faces))
    }

    pub fn collectable_closing_faces(&self, removed_faces: &Vec<usize>) -> Result<Vec<[usize; 3]>> {
        let mut palve_paths = {
            let palve_path = self.external_path.alveolae_path(&self.skeleton_interface)?;
            vec![palve_path]
        };
        let mut fixed_internal = HashSet::new();
        for internal_path in self.internal_paths.iter() {
            let palve_internal = internal_path.alveolae_path(&self.skeleton_interface)?;
            for ind_palve in palve_internal {
                fixed_internal.insert(ind_palve);
            }
        }
        let mut unfaced_hedges = HashSet::new();
        for &ind_fac in removed_faces {
            let fac = self.skeleton_interface.get_mesh().get_face(ind_fac)?;
            let [hedg0, hedg1, hedg2] = fac.halfedges();
            unfaced_hedges.insert(hedg0.halfedge());
            unfaced_hedges.insert(hedg1.halfedge());
            unfaced_hedges.insert(hedg2.halfedge());
        }
        let mut faces = Vec::new();
        loop {
            println!("new iter");
            for path in palve_paths.iter() {
                for &ind_palve in path.iter() {
                    let palve = self
                        .skeleton_interface
                        .get_partial_alveola_uncheck(ind_palve);
                    let seg = palve.alveola().delaunay_segment();
                    let seg = if seg[0] == palve.corner() {
                        seg
                    } else {
                        [seg[1], seg[0]]
                    };
                    print!("({} -> {}), ", seg[0], seg[1]);
                }
                println!("");
            }
            println!("");

            if let Some(palve_path) = palve_paths.pop() {
                if palve_path.is_empty() {
                    continue;
                }
                let nb_palve = palve_path.len();
                let mut operation_done = false;
                let mut all_fixed = true;
                for ind in 0..nb_palve {
                    let ind_palve = palve_path[ind];
                    if fixed_internal.contains(&ind) {
                        continue;
                    }
                    all_fixed = false;

                    let palve = self
                        .skeleton_interface
                        .get_partial_alveola_uncheck(ind_palve);
                    let seg = palve.alveola().delaunay_segment();
                    let seg = if seg[0] == palve.corner() {
                        seg
                    } else {
                        [seg[1], seg[0]]
                    };
                    if unfaced_hedges.contains(&seg) {
                        continue;
                    }

                    let ind_palve_opp = palve.partial_alveola_opposite().ind();
                    let opt_position = palve_path.iter().position(|&ind| ind == ind_palve_opp);
                    if let Some(position) = opt_position {
                        let mut path1 = Vec::new();
                        let mut path2 = Vec::new();
                        for i in 0..nb_palve {
                            if i < ind || i > position {
                                path1.push(palve_path[i]);
                            }
                            if i > ind && i < position {
                                path2.push(palve_path[i]);
                            }
                        }
                        palve_paths.push(path1);
                        palve_paths.push(path2);
                        operation_done = true;
                        break;
                    }
                }
                if all_fixed {
                    operation_done = true
                }
                if !operation_done {
                    for ind in 0..nb_palve {
                        let ind_palve = palve_path[ind];

                        let palve = self
                            .skeleton_interface
                            .get_partial_alveola_uncheck(ind_palve);
                        let ind_prev = (ind - 1 + nb_palve) % nb_palve;
                        let ind_palve_prev = palve_path[ind_prev];
                        let palve_prev = self
                            .skeleton_interface
                            .get_partial_alveola_uncheck(ind_palve_prev);
                        let ind_vertex = palve_prev.corner();
                        for pedge in palve.partial_edges() {
                            if !pedge.edge().is_full() {
                                continue;
                            }
                            let tri = pedge.edge().delaunay_triangle();
                            if tri.contains(&ind_vertex) {
                                let seg = palve.alveola().delaunay_segment();
                                if ind_vertex == seg[0] || ind_vertex == seg[1] {
                                    continue;
                                }
                                let seg = if seg[0] == palve.corner() {
                                    seg
                                } else {
                                    [seg[1], seg[0]]
                                };
                                let mut tri = [seg[0], seg[1], ind_vertex];
                                tri.sort();
                                faces.push([seg[0], seg[1], ind_vertex]);

                                let seg2 = if seg[1] < ind_vertex {
                                    [seg[1], ind_vertex]
                                } else {
                                    [ind_vertex, seg[1]]
                                };
                                let &ind_alveola2 =
                                    self.skeleton_interface.del_seg.get(&seg2).ok_or(
                                        anyhow::Error::msg(
                                            "Second partial alveola not in skeleton",
                                        ),
                                    )?;
                                let alve2 =
                                    self.skeleton_interface.get_alveola_uncheck(ind_alveola2);
                                let ind_palve2 =
                                    if alve2.partial_alveolae()[0].corner() == ind_vertex {
                                        alve2.partial_alveolae()[0].ind()
                                    } else {
                                        alve2.partial_alveolae()[1].ind()
                                    };
                                let mut palve_path_new = Vec::new();
                                for i in 0..palve_path.len() {
                                    if i != ind && i != ind_prev {
                                        palve_path_new.push(palve_path[i]);
                                    }
                                    if i == ind {
                                        palve_path_new.push(ind_palve2);
                                    }
                                }
                                palve_paths.push(palve_path_new);

                                operation_done = true;
                                break;
                            }
                        }
                        if operation_done {
                            break;
                        }
                    }
                }

                if !operation_done {
                    return Ok(faces);
                }
            } else {
                break;
            }
            if faces.len() > self.skeleton_interface.get_mesh().get_nb_faces() {
                return Ok(faces);
            }
        }
        Ok(faces)
    }
}

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

    fn last_paths_fusion(
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

    fn expand_last_face(
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
        .basis_spheres(&skeleton_separation.skeleton_interface())?;
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

        if faces.len()
            > skeleton_separation
                .skeleton_interface()
                .get_mesh()
                .get_nb_faces()
                * 2
        {
            return Ok(None);
            // return Err(anyhow::Error::msg("Too much faces collected somehow"));
        }

        if last_hedge_deletion(&mut mesh_paths_external, &skeleton_separation)? {
            continue;
        }

        if last_paths_fusion(
            &mut mesh_paths_external,
            &mut mesh_paths_internal,
            &skeleton_separation,
        )? {
            continue;
        }

        if expand_last_face(
            &mut mesh_paths_external,
            &skeleton_separation,
            &center_mat,
            &radius_mat,
            epsilon,
            &mut faces,
        )? {
            continue;
        }
        return Ok(None);
    }

    Ok(Some(faces))
}
