use crate::algorithm::skeleton_interface::SkeletonInterface3D;
use anyhow::Result;
use nalgebra::{MatrixXx1, MatrixXx3, Vector3};

use super::skeleton_interface::IterPartialEdge;
use crate::geometry::geometry_operations;

#[derive(Copy, Clone)]
pub enum PathPart {
    PartialNode(usize),
    PartialEdge(usize),
}

pub enum State {
    Computing,
    Closed,
}

pub struct SkeletonPath<'a, 'b> {
    skeleton_interface: &'b mut SkeletonInterface3D<'a>,
    components: Vec<PathPart>,
    opt_ind_pedge_last: Option<usize>,
}

impl<'a, 'b> SkeletonPath<'a, 'b> {
    pub fn new(
        skeleton_interface: &'b mut SkeletonInterface3D<'a>,
        ind_pedge: usize,
    ) -> SkeletonPath<'a, 'b> {
        SkeletonPath {
            skeleton_interface,
            components: Vec::new(),
            opt_ind_pedge_last: Some(ind_pedge),
        }
    }

    pub fn skeleton_interface(&self) -> &SkeletonInterface3D {
        self.skeleton_interface
    }

    pub fn mesh_path(&self) -> Vec<usize> {
        let mut path = Vec::new();
        for ind1 in 0..self.components.len() {
            let ind2 = (ind1 + 1) % self.components.len();
            match (self.components[ind1], self.components[ind2]) {
                (PathPart::PartialNode(ind_pnode), PathPart::PartialNode(_)) => {
                    path.push(
                        self.skeleton_interface
                            .get_partial_node_uncheck(ind_pnode)
                            .corner(),
                    );
                }
                (_, _) => (),
            }
        }
        path
    }

    pub fn append_last(&mut self) -> Result<State> {
        if let Some(ind_pedge_last) = self.opt_ind_pedge_last {
            let edge = self
                .skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_last)
                .edge();
            if !edge.is_computed() {
                let ind_edge = edge.ind();
                self.skeleton_interface.propagate_edge(ind_edge)?;
            }
            self.components.push(PathPart::PartialEdge(ind_pedge_last));

            let ind_pedge_next = self
                .skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_last)
                .partial_edge_next()
                .ok_or(anyhow::Error::msg("No next partial edge"))?
                .ind();

            self.opt_ind_pedge_last = Some(ind_pedge_next);

            self.check_loop()
        } else {
            Ok(State::Closed)
        }
    }

    fn check_loop(&mut self) -> Result<State> {
        if let Some(ind_pedge_last) = self.opt_ind_pedge_last {
            let part_first = self.components.first().unwrap();
            let looped = match part_first {
                &PathPart::PartialNode(_) => todo!(),
                &PathPart::PartialEdge(ind_pedge) => ind_pedge == ind_pedge_last,
            };

            if looped {
                self.opt_ind_pedge_last = None
            }
        }

        if self.opt_ind_pedge_last.is_none() {
            Ok(State::Closed)
        } else {
            Ok(State::Computing)
        }
    }

    pub fn rotate_last(&mut self) -> Result<State> {
        if let Some(ind_pedge_last) = self.opt_ind_pedge_last {
            let pedge = self
                .skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_last);
            let pedge_neigh = pedge.partial_edge_neighbor();

            let alve = pedge_neigh.partial_alveola().alveola();
            let pedge_next = if alve.is_full() {
                pedge_neigh
                    .partial_edge_next()
                    .ok_or(anyhow::Error::msg("No next partial edge"))?
            } else {
                pedge_neigh.partial_edge_opposite()
            };
            let ind_pnode = pedge
                .partial_node_first()
                .ok_or(anyhow::Error::msg("No first node"))?
                .ind();
            if let Some(&plast) = self.components.last() {
                if let PathPart::PartialNode(nod) = plast {
                    if nod != ind_pnode {
                        self.components.push(PathPart::PartialNode(ind_pnode));
                    }
                } else {
                    self.components.push(PathPart::PartialNode(ind_pnode));
                }
            } else {
                self.components.push(PathPart::PartialNode(ind_pnode));
            }
            let ind_pedge_new = pedge_next.ind();
            self.opt_ind_pedge_last = Some(ind_pedge_new);
            self.check_loop()
        } else {
            Ok(State::Closed)
        }
    }

    pub fn ind_last_partial_edge(&self) -> Option<usize> {
        self.opt_ind_pedge_last
    }

    pub fn last_partial_edge(&self) -> Option<IterPartialEdge> {
        if let Some(ind_pedge_last) = self.opt_ind_pedge_last {
            Some(
                self.skeleton_interface
                    .get_partial_edge_uncheck(ind_pedge_last),
            )
        } else {
            None
        }
    }

    pub fn ind_partial_edges(&self) -> Vec<usize> {
        self.components
            .iter()
            .filter_map(|pp| match pp {
                &PathPart::PartialEdge(ind_pedge) => Some(ind_pedge),
                &PathPart::PartialNode(_) => None,
            })
            .collect()
    }

    pub fn ind_alveolae(&self) -> Vec<usize> {
        self.components
            .iter()
            .filter_map(|pp| match pp {
                &PathPart::PartialEdge(ind_pedge) => Some(
                    self.skeleton_interface
                        .get_partial_edge_uncheck(ind_pedge)
                        .partial_alveola()
                        .alveola()
                        .ind(),
                ),
                &PathPart::PartialNode(_) => None,
            })
            .collect()
    }

    pub fn closable_path(&self) -> Result<bool> {
        let mut has_deg1 = false;
        for ind in 0..self.components.len() {
            let ind_next = (ind + 1) % self.components.len();
            match (self.components[ind], self.components[ind_next]) {
                (PathPart::PartialEdge(ind_pedge), PathPart::PartialNode(_)) => {
                    let pedge = self.skeleton_interface.get_partial_edge_uncheck(ind_pedge);
                    if pedge.partial_alveola().alveola().label().is_some() {
                        return Ok(false);
                    }
                    let deg = pedge
                        .partial_edge_next()
                        .ok_or(anyhow::Error::msg("Next edge does not exist"))?
                        .edge()
                        .degree();
                    if deg != 1 && deg != 2 {
                        return Ok(false);
                    }
                    if deg == 1 {
                        has_deg1 = true;
                    }
                }
                (_, _) => (),
            };
        }
        return Ok(has_deg1);
    }

    pub fn nodes(&self) -> Vec<usize> {
        let mut nodes: Vec<usize> = self
            .components
            .iter()
            .filter_map(|&cmp| match cmp {
                PathPart::PartialNode(ind_pnode) => Some(
                    self.skeleton_interface
                        .get_partial_node_uncheck(ind_pnode)
                        .node()
                        .ind(),
                ),
                PathPart::PartialEdge(_) => None,
            })
            .collect();
        nodes.sort();
        nodes
    }

    pub fn basis_spheres(&self) -> Result<(MatrixXx3<f32>, MatrixXx1<f32>)> {
        let ind_nodes = self.nodes();
        let mut center_mat = MatrixXx3::<f32>::zeros(ind_nodes.len());
        let mut radius_mat = MatrixXx1::<f32>::zeros(ind_nodes.len());

        for i in 0..ind_nodes.len() {
            let ind_node = ind_nodes[i];
            let tet_vert: Vec<Vector3<f32>> = self
                .skeleton_interface
                .get_node_uncheck(ind_node)
                .delaunay_tetrahedron()
                .iter()
                .map(|&ind| {
                    self.skeleton_interface
                        .get_mesh()
                        .get_vertex(ind)
                        .unwrap()
                        .vertex()
                })
                .collect();
            let (center, radius) = geometry_operations::center_and_radius(
                [tet_vert[0], tet_vert[1], tet_vert[2], tet_vert[3]],
                None,
            )
            .ok_or(anyhow::Error::msg("Could not find sphere center"))?;
            center_mat[(i, 0)] = center[0];
            center_mat[(i, 1)] = center[1];
            center_mat[(i, 2)] = center[2];
            radius_mat[i] = radius;
        }

        Ok((center_mat, radius_mat))
    }

    pub fn collect_mesh_faces_index(&self, epsilon: f32) -> Result<Option<Vec<usize>>> {
        let vec_bnd_path_vert = self.mesh_path();
        let (center_mat, radius_mat) = self.basis_spheres()?;
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
        let mesh_path = self.mesh_path();
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

    pub fn compute_debug_mesh(&mut self) -> Result<()> {
        for ind1 in 0..self.components.len() {
            if let PathPart::PartialNode(ind_pnode1) = self.components[ind1] {
                let ind2 = (ind1 + 1) % self.components.len();
                if let PathPart::PartialNode(ind_pnode2) = self.components[ind2] {
                    let pnode1 = self.skeleton_interface.get_partial_node_uncheck(ind_pnode1);
                    let pnode2 = self.skeleton_interface.get_partial_node_uncheck(ind_pnode2);
                    let vert1 = self
                        .skeleton_interface
                        .get_mesh()
                        .get_vertex(pnode1.corner())?
                        .vertex();
                    let vert2 = self
                        .skeleton_interface
                        .get_mesh()
                        .get_vertex(pnode2.corner())?
                        .vertex();
                    let tet_vert3: Vec<Vector3<f32>> = pnode1
                        .node()
                        .delaunay_tetrahedron()
                        .iter()
                        .map(|&ind| {
                            self.skeleton_interface
                                .get_mesh()
                                .get_vertex(ind)
                                .unwrap()
                                .vertex()
                        })
                        .collect();
                    let (vert3, _) = geometry_operations::center_and_radius(
                        [tet_vert3[0], tet_vert3[1], tet_vert3[2], tet_vert3[3]],
                        None,
                    )
                    .ok_or(anyhow::Error::msg("Could not find sphere center"))?;
                    let ind_vertex1 = self.skeleton_interface.debug_mesh.add_vertex(&vert1);
                    let ind_vertex2 = self.skeleton_interface.debug_mesh.add_vertex(&vert2);
                    let ind_vertex3 = self.skeleton_interface.debug_mesh.add_vertex(&vert3);
                    self.skeleton_interface.debug_mesh.add_face(
                        ind_vertex1,
                        ind_vertex2,
                        ind_vertex3,
                    )?;
                }
            }
            if let PathPart::PartialEdge(pedge) = self.components[ind1] {
                let pedge = self.skeleton_interface.get_partial_edge_uncheck(pedge);
                let corner1 = pedge.corner();
                let vert1 = self
                    .skeleton_interface
                    .get_mesh()
                    .get_vertex(corner1)?
                    .vertex();
                let tet_vert2: Vec<Vector3<f32>> = pedge
                    .partial_node_first()
                    .ok_or(anyhow::Error::msg("No first node"))?
                    .node()
                    .delaunay_tetrahedron()
                    .iter()
                    .map(|&ind| {
                        self.skeleton_interface
                            .get_mesh()
                            .get_vertex(ind)
                            .unwrap()
                            .vertex()
                    })
                    .collect();
                let (vert2, _) = geometry_operations::center_and_radius(
                    [tet_vert2[0], tet_vert2[1], tet_vert2[2], tet_vert2[3]],
                    None,
                )
                .ok_or(anyhow::Error::msg("Flat tetrahedron"))?;
                let tet_vert3: Vec<Vector3<f32>> = pedge
                    .partial_node_last()
                    .ok_or(anyhow::Error::msg("No first node"))?
                    .node()
                    .delaunay_tetrahedron()
                    .iter()
                    .map(|&ind| {
                        self.skeleton_interface
                            .get_mesh()
                            .get_vertex(ind)
                            .unwrap()
                            .vertex()
                    })
                    .collect();
                let (vert3, _) = geometry_operations::center_and_radius(
                    [tet_vert3[0], tet_vert3[1], tet_vert3[2], tet_vert3[3]],
                    None,
                )
                .ok_or(anyhow::Error::msg("Flat tetrahedron"))?;
                let ind_vertex1 = self.skeleton_interface.debug_mesh.add_vertex(&vert1);
                let ind_vertex2 = self.skeleton_interface.debug_mesh.add_vertex(&vert2);
                let ind_vertex3 = self.skeleton_interface.debug_mesh.add_vertex(&vert3);
                self.skeleton_interface.debug_mesh.add_face(
                    ind_vertex1,
                    ind_vertex2,
                    ind_vertex3,
                )?;
            }
        }

        let opt_faces = self.collect_closing_faces()?;

        if let Some(faces) = opt_faces {
            for [ind_vertex1, ind_vertex2, ind_vertex3] in faces {
                let vert1 = self
                    .skeleton_interface
                    .get_mesh()
                    .get_vertex(ind_vertex1)?
                    .vertex();
                let vert2 = self
                    .skeleton_interface
                    .get_mesh()
                    .get_vertex(ind_vertex2)?
                    .vertex();
                let vert3 = self
                    .skeleton_interface
                    .get_mesh()
                    .get_vertex(ind_vertex3)?
                    .vertex();
                let ind_v1 = self.skeleton_interface.debug_mesh.add_vertex(&vert1);
                let ind_v2 = self.skeleton_interface.debug_mesh.add_vertex(&vert2);
                let ind_v3 = self.skeleton_interface.debug_mesh.add_vertex(&vert3);
                self.skeleton_interface
                    .debug_mesh
                    .add_face(ind_v1, ind_v2, ind_v3)?;
            }
        }

        Ok(())
    }
}
