use crate::algorithm::skeleton_interface::SkeletonInterface3D;
use anyhow::Result;
use nalgebra::Vector3;

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
        self.components
            .iter()
            .map(|part| match part {
                &PathPart::PartialNode(pnode) => self
                    .skeleton_interface
                    .get_partial_node_uncheck(pnode)
                    .corner(),
                &PathPart::PartialEdge(pedge) => self
                    .skeleton_interface
                    .get_partial_edge_uncheck(pedge)
                    .corner(),
            })
            .collect()
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
            let ind_pedge_new = pedge_next.ind();
            self.components.push(PathPart::PartialNode(ind_pnode));
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
        if self.components.len() > 500 {
            return Ok(false);
        }
        for ind in 0..self.components.len() {
            let ind_next = (ind + 1) % self.components.len();
            match (self.components[ind], self.components[ind_next]) {
                (PathPart::PartialEdge(ind_pedge), PathPart::PartialNode(_)) => {
                    let pedge = self.skeleton_interface.get_partial_edge_uncheck(ind_pedge);
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
        // for ind in 0..self.components.len() {
        //     let ind_next = (ind + 1) % self.components.len();
        //     match (self.components[ind], self.components[ind_next]) {
        //         (PathPart::PartialEdge(ind_pedge), PathPart::PartialNode(_)) => {
        //             let edge = self
        //                 .skeleton_interface
        //                 .get_partial_edge_uncheck(ind_pedge)
        //                 .partial_edge_next()
        //                 .ok_or(anyhow::Error::msg("Next edge does not exist"))?
        //                 .edge();
        //             if edge.degree() < 3 {
        //                 let [ind_vertex1, ind_vertex2, ind_vertex3] = edge.delaunay_triangle();
        //                 let vert1 = self
        //                     .skeleton_interface
        //                     .get_mesh()
        //                     .get_vertex(ind_vertex1)?
        //                     .vertex();
        //                 let vert2 = self
        //                     .skeleton_interface
        //                     .get_mesh()
        //                     .get_vertex(ind_vertex2)?
        //                     .vertex();
        //                 let vert3 = self
        //                     .skeleton_interface
        //                     .get_mesh()
        //                     .get_vertex(ind_vertex3)?
        //                     .vertex();
        //                 let ind_v1 = self.skeleton_interface.debug_mesh.add_vertex(&vert1);
        //                 let ind_v2 = self.skeleton_interface.debug_mesh.add_vertex(&vert2);
        //                 let ind_v3 = self.skeleton_interface.debug_mesh.add_vertex(&vert3);
        //                 self.skeleton_interface
        //                     .debug_mesh
        //                     .add_face(ind_v1, ind_v2, ind_v3)?;
        //             }
        //         }
        //         (_, _) => (),
        //     };
        // }

        let faces = self.collect_faces()?;

        for ind_face in faces {
            let [ind_vertex1, ind_vertex2, ind_vertex3] = self
                .skeleton_interface
                .get_mesh()
                .get_face(ind_face)?
                .vertices_inds();
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

        Ok(())
    }

    pub fn collect_faces(&self) -> Result<Vec<usize>> {
        let vec_bnd_path_vert = self.mesh_path();
        let mut mesh_path_hedge = Vec::new();
        for ind1 in 0..vec_bnd_path_vert.len() {
            let ind2 = (ind1 + 1) % vec_bnd_path_vert.len();
            let ind_vertex1 = vec_bnd_path_vert[ind1];
            let ind_vertex2 = vec_bnd_path_vert[ind2];

            if ind_vertex1 != ind_vertex2 {
                let res_ind_edge = self
                    .skeleton_interface
                    .get_mesh()
                    .is_edge_in(ind_vertex1, ind_vertex2)
                    .ok_or(anyhow::Error::msg("Part of the path not on the boundary"));
                match res_ind_edge {
                    Err(e) => {
                        if self
                            .skeleton_interface
                            .closing_mesh
                            .is_edge_in(ind_vertex1, ind_vertex2)
                            .is_some()
                        {
                            println!("In new mesh");
                        } else {
                            println!("Not in new mesh");
                        }
                        return Err(e);
                    }
                    Ok(hedge) => mesh_path_hedge.push(hedge.ind()),
                }
            }
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
            if faces.len() > self.skeleton_interface.get_mesh().get_nb_faces() {
                println!("Too much faces collected somehow");
                faces.sort();
                faces.dedup();
                return Ok(faces);
            }
        }

        Ok(faces)
    }

    pub fn close_path(&mut self) -> Result<()> {
        for ind in 0..self.components.len() {
            let ind_next = (ind + 1) % self.components.len();
            match (self.components[ind], self.components[ind_next]) {
                (PathPart::PartialEdge(ind_pedge), PathPart::PartialNode(_)) => {
                    let pedge = self.skeleton_interface.get_partial_edge_uncheck(ind_pedge);
                    let [ind_vertex1, ind_vertex2, ind_vertex3] = pedge
                        .partial_edge_next()
                        .ok_or(anyhow::Error::msg("Next edge does not exist"))?
                        .edge()
                        .delaunay_triangle();
                    self.skeleton_interface
                        .close_face(ind_vertex1, ind_vertex2, ind_vertex3)?
                }
                (_, _) => (),
            };
        }
        Ok(())
    }
}
