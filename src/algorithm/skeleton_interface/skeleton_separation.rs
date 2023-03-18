use crate::algorithm::skeleton_interface::SkeletonInterface3D;
use anyhow::Result;

use super::skeleton_path::SkeletonPath;

pub struct SkeletonSeparation<'a, 'b> {
    skeleton_interface: &'b mut SkeletonInterface3D<'a>,
    external_path: SkeletonPath,
}

impl<'a, 'b> SkeletonSeparation<'a, 'b> {
    pub fn new(
        skeleton_interface: &'b mut SkeletonInterface3D<'a>,
        ind_pedge: usize,
    ) -> SkeletonSeparation<'a, 'b> {
        SkeletonSeparation {
            skeleton_interface,
            external_path: SkeletonPath::new(ind_pedge),
        }
    }

    pub fn skeleton_interface(&self) -> &SkeletonInterface3D {
        self.skeleton_interface
    }

    pub fn follow_external_path(&mut self) -> Result<()> {
        self.external_path
            .follow_singular_path(&mut self.skeleton_interface)
    }

    pub fn closable_path(&self) -> Result<bool> {
        self.external_path.closable_path(&self.skeleton_interface)
    }

    pub fn collect_mesh_faces_index(&self, epsilon: f32) -> Result<Option<Vec<usize>>> {
        let vec_bnd_path_vert = self.external_path.mesh_path(&self.skeleton_interface);
        let (center_mat, radius_mat) =
            self.external_path.basis_spheres(&self.skeleton_interface)?;
        let mut mesh_path_hedge = Vec::new();
        for ind1 in 0..vec_bnd_path_vert.len() {
            let ind2 = (ind1 + 1) % vec_bnd_path_vert.len();
            let ind_vertex1 = vec_bnd_path_vert[ind1];
            let ind_vertex2 = vec_bnd_path_vert[ind2];

            let hedge = self
                .skeleton_interface
                .get_mesh()
                .is_edge_in(ind_vertex1, ind_vertex2)
                .ok_or(anyhow::Error::msg("Part of the path not on the boundary"))?;
            mesh_path_hedge.push(hedge.ind());
        }

        let mut mesh_paths_hedge = vec![mesh_path_hedge];
        let mut faces = Vec::new();
        loop {
            // println!("new iter");
            // for path in mesh_paths_hedge.iter() {
            //     for &ind_he in path.iter() {
            //         let hedge = self.skeleton_interface.get_mesh().get_halfedge(ind_he)?;
            //         print!(
            //             "({} -> {}), ",
            //             hedge.first_vertex().ind(),
            //             hedge.last_vertex().ind()
            //         );
            //     }
            //     println!("");
            // }
            // println!("");

            if let Some(mut mesh_path_edge) = mesh_paths_hedge.pop() {
                if let Some(ind_hedge) = mesh_path_edge.pop() {
                    let hedge = self.skeleton_interface.get_mesh().get_halfedge(ind_hedge)?;
                    let ind_hedge_opp = hedge.opposite_halfedge().unwrap().ind();
                    let opt_position = mesh_path_edge.iter().position(|&ind| ind == ind_hedge_opp);
                    if let Some(position) = opt_position {
                        if position != 0 {
                            let mut path1 = vec![0 as usize; position];
                            path1.copy_from_slice(&mesh_path_edge[..position]);
                            mesh_paths_hedge.push(path1);
                        }
                        if position != mesh_path_edge.len() - 1 {
                            let mut path2 = vec![0 as usize; mesh_path_edge.len() - 1 - position];
                            path2.copy_from_slice(&mesh_path_edge[position + 1..]);
                            mesh_paths_hedge.push(path2);
                        }
                    } else {
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
                            return Ok(None);
                        }

                        let face = hedge.face().unwrap();

                        faces.push(face.ind());
                        let hedge_rep1 =
                            hedge.prev_halfedge().unwrap().opposite_halfedge().unwrap();
                        let hedge_rep2 =
                            hedge.next_halfedge().unwrap().opposite_halfedge().unwrap();
                        mesh_path_edge.push(hedge_rep1.ind());
                        mesh_path_edge.push(hedge_rep2.ind());
                        mesh_paths_hedge.push(mesh_path_edge);
                    }
                }
            } else {
                break;
            }
            if faces.len() > self.skeleton_interface.get_mesh().get_nb_faces() * 2 {
                return Err(anyhow::Error::msg("Too much faces collected somehow"));
            }
        }

        Ok(Some(faces))
    }

    pub fn collect_closing_faces(&self) -> Result<Option<Vec<[usize; 3]>>> {
        let mesh_path = self.external_path.mesh_path(&self.skeleton_interface);
        let mut palve_path = Vec::new();
        for ind1 in 0..mesh_path.len() {
            let ind2 = (ind1 + 1) % mesh_path.len();
            let ind_vertex1 = mesh_path[ind1];
            let ind_vertex2 = mesh_path[ind2];
            let seg = if ind_vertex1 < ind_vertex2 {
                [ind_vertex1, ind_vertex2]
            } else {
                [ind_vertex2, ind_vertex1]
            };

            let &ind_alveola = self
                .skeleton_interface
                .del_seg
                .get(&seg)
                .ok_or(anyhow::Error::msg("Alveola not in skeleton"))?;
            let alve = self.skeleton_interface.get_alveola_uncheck(ind_alveola);
            let ind_palve = if alve.partial_alveolae()[0].corner() == ind_vertex1 {
                alve.partial_alveolae()[0].ind()
            } else {
                alve.partial_alveolae()[1].ind()
            };
            palve_path.push(ind_palve);
        }

        let mut palve_paths = vec![palve_path];
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
                for ind in 0..nb_palve {
                    let ind_palve = palve_path[ind];

                    let palve = self
                        .skeleton_interface
                        .get_partial_alveola_uncheck(ind_palve);
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
}
