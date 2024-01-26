use anyhow::Result;
use simple_delaunay_lib::delaunay_3d::delaunay_struct_3d::DelaunayStructure3D;
use simple_delaunay_lib::delaunay_3d::simplicial_struct_3d::Node;
use std::collections::{HashMap, HashSet};

use crate::mesh3d::mesh_operations;
use crate::mesh3d::{manifold_mesh3d, ManifoldMesh3D};

pub type Edge = [usize; 2];
pub type Triangle = [usize; 3];
pub type Tetrahedron = [usize; 4];

/// Delaunay structure
pub struct DelaunayInterface<'a> {
    mesh: &'a mut ManifoldMesh3D,
    del_struct: DelaunayStructure3D,

    vertex_edges: Vec<Vec<(usize, usize)>>,

    non_del_edges: Vec<usize>,
    non_del_faces: Vec<usize>,

    initial_vertices_number: usize,
}

impl<'a> DelaunayInterface<'a> {
    fn generate_struct(&mut self) -> Result<()> {
        let mut points = Vec::new();
        for v in self.mesh.vertex_indices() {
            let vert = self.mesh.get_vertex(v)?.vertex();
            points.push([vert[0], vert[1], vert[2]]);
            self.vertex_edges.push(Vec::new());
        }
        self.del_struct.insert_vertices(&points, true)?;

        for ind_tet in 0..self.del_struct.get_simplicial().get_nb_tetrahedra() {
            let tetra = self.del_struct.get_simplicial().get_tetrahedron(ind_tet)?;

            for tri in tetra.halftriangles() {
                let hes = tri.halfedges();
                for i in 0..3 {
                    if let (Node::Value(i1), Node::Value(_)) =
                        (hes[i].first_node(), hes[i].last_node())
                    {
                        self.vertex_edges[i1].push((tri.ind(), i));
                    }
                }
            }
        }

        Ok(())
    }

    fn insert_vertex(&mut self, ind_vertex: usize, near_to: usize) -> Result<()> {
        let vert = self.mesh.get_vertex(ind_vertex)?.vertex();
        self.del_struct.insert_vertex(
            [vert[0] as f64, vert[1] as f64, vert[2] as f64],
            Some(near_to),
        )?;
        self.vertex_edges.push(Vec::new());

        let tet_update = self
            .del_struct
            .get_simplicial()
            .get_tetrahedra_containing(&Node::Value(self.vertex_edges.len() - 1));

        let mut vert_to_check = HashSet::new();
        for tetra in tet_update {
            for tri in tetra.halftriangles() {
                let hes = tri.halfedges();
                for i in 0..3 {
                    if let (Node::Value(i1), Node::Value(i2)) =
                        (hes[i].first_node(), hes[i].last_node())
                    {
                        self.vertex_edges[i1].push((tri.ind(), i));
                        vert_to_check.insert(i1);
                        vert_to_check.insert(i2);
                    }
                }
            }
        }

        for iv in vert_to_check {
            self.vertex_edges[iv] = self.vertex_edges[iv]
                .iter()
                .filter_map(|&(it, i)| {
                    if let Ok(tri) = self.del_struct.get_simplicial().get_halftriangle(it) {
                        if tri.halfedges()[i].first_node().equals(&Node::Value(iv)) {
                            Some((it, i))
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                })
                .collect();
        }

        Ok(())
    }

    fn fill_non_del(&mut self) -> () {
        self.non_del_edges.clear();
        self.non_del_faces.clear();

        for (&ind_fac, _) in self.mesh.faces() {
            let face = self.mesh.get_face(ind_fac).unwrap();
            let face_vert = face.vertices_inds();
            if !self.is_face_in(&face_vert) {
                self.non_del_faces.push(ind_fac);
                for ind_he in face.face_halfedges() {
                    let edge = self.mesh.get_halfedge(ind_he).unwrap();
                    let edge_vert = edge.halfedge();
                    if !self.is_edge_in(&edge_vert) {
                        self.non_del_edges.push(ind_he);
                    };
                }
            };
        }
    }

    /// Creates Delaunay structure from mesh
    pub fn from_mesh(mesh: &'a mut ManifoldMesh3D) -> Result<DelaunayInterface<'a>> {
        let initial_vertices_number = mesh.get_nb_vertices();
        let mut deltet = DelaunayInterface {
            mesh,
            del_struct: DelaunayStructure3D::new(),
            vertex_edges: Vec::new(),
            non_del_edges: Vec::new(),
            non_del_faces: Vec::new(),
            initial_vertices_number,
        };

        deltet.generate_struct()?;

        deltet.fill_non_del();

        Ok(deltet)
    }

    /// Mesh getter
    pub fn get_mesh(&self) -> &ManifoldMesh3D {
        self.mesh
    }

    /// Tetrahedra set getter
    pub fn get_faces(&self) -> HashMap<Triangle, Vec<Tetrahedron>> {
        let mut face_set = HashMap::new();
        for ind_tet in 0..self.del_struct.get_simplicial().get_nb_tetrahedra() {
            let tetra = self
                .del_struct
                .get_simplicial()
                .get_tetrahedron(ind_tet)
                .unwrap();
            if let [Node::Value(i1), Node::Value(i2), Node::Value(i3), Node::Value(i4)] =
                tetra.nodes()
            {
                let mut tetra_ind = [i1, i2, i3, i4];
                tetra_ind.sort();

                face_set
                    .entry([tetra_ind[0], tetra_ind[1], tetra_ind[2]])
                    .or_insert(Vec::new())
                    .push(tetra_ind);
                face_set
                    .entry([tetra_ind[0], tetra_ind[1], tetra_ind[3]])
                    .or_insert(Vec::new())
                    .push(tetra_ind);
                face_set
                    .entry([tetra_ind[0], tetra_ind[2], tetra_ind[3]])
                    .or_insert(Vec::new())
                    .push(tetra_ind);
                face_set
                    .entry([tetra_ind[1], tetra_ind[2], tetra_ind[3]])
                    .or_insert(Vec::new())
                    .push(tetra_ind);
            }
        }
        face_set
    }

    /// Checks if vertex was an original mesh vertex
    pub fn is_original_vertex(&self, ind_vertex: usize) -> bool {
        ind_vertex < self.initial_vertices_number
    }

    /// Checks if edge is in Delaunay
    fn is_edge_in(&self, edge: &Edge) -> bool {
        self.vertex_edges[edge[0]]
            .iter()
            .position(|&(it, i)| {
                if let Ok(tri) = self.del_struct.get_simplicial().get_halftriangle(it) {
                    tri.halfedges()[i].last_node().equals(&Node::Value(edge[1]))
                } else {
                    false
                }
            })
            .is_some()
    }

    /// Checks if face is in Delaunay
    fn is_face_in(&self, face: &Triangle) -> bool {
        self.vertex_edges[face[0]]
            .iter()
            .position(|&(it, i)| {
                if let Ok(tri) = self.del_struct.get_simplicial().get_halftriangle(it) {
                    let he = tri.halfedges()[i];
                    if he.last_node().equals(&Node::Value(face[1])) {
                        he.next().last_node().equals(&Node::Value(face[2]))
                    } else {
                        false
                    }
                } else {
                    false
                }
            })
            .is_some()
    }

    /// Count number of non Delaunay halfedges
    pub fn count_non_del_halfedges(&mut self) -> usize {
        if self.non_del_edges.len() == 0 {
            self.fill_non_del();
        }
        self.non_del_edges.len()
    }

    /// Count number of non Delaunay faces
    pub fn count_non_del_faces(&mut self) -> usize {
        if self.non_del_faces.len() == 0 {
            self.fill_non_del();
        }
        self.non_del_faces.len()
    }

    /// Gets first globally non Delaunay halfedge, starting from a shift
    pub fn get_non_del_halfedge(&mut self) -> Result<Option<manifold_mesh3d::IterHalfEdge>> {
        let mut length: Vec<(usize, f64)> = self
            .non_del_edges
            .iter()
            .map(|&ind_he| {
                if let Ok(edge) = self.mesh.get_halfedge(ind_he) {
                    (
                        ind_he,
                        (edge.first_vertex().vertex() - edge.last_vertex().vertex()).norm_squared(),
                    )
                } else {
                    (ind_he, 0.0)
                }
            })
            .collect();
        length.sort_by(|(_, l1), (_, l2)| l1.partial_cmp(l2).unwrap());
        self.non_del_edges = length.iter().map(|&(ind_he, _)| ind_he).collect();

        loop {
            if let Some(ind_he) = self.non_del_edges.pop() {
                if let Ok(edge) = self.mesh.get_halfedge(ind_he) {
                    let seg = edge.halfedge();
                    if !self.is_edge_in(&seg) {
                        let he = self.mesh.get_halfedge(ind_he)?;
                        break Ok(Some(he));
                    };
                }
            } else {
                break Ok(None);
            }
        }
    }

    /// Gets first globally non Delaunay face, starting from a shift
    pub fn get_non_del_face(&mut self) -> Result<Option<manifold_mesh3d::IterFace>> {
        loop {
            if let Some(ind_fac) = self.non_del_faces.pop() {
                if let Ok(face) = self.mesh.get_face(ind_fac) {
                    let face_vert = face.vertices_inds();
                    if !self.is_face_in(&face_vert) {
                        break Ok(Some(face));
                    };
                }
            } else {
                break Ok(None);
            }
        }
    }

    /// Flips given halfedge
    pub fn flip_halfedge(&mut self, ind_halfedge: usize) -> Result<bool> {
        mesh_operations::flip_halfedge(self.mesh, ind_halfedge)
    }

    /// Splits given halfedge
    pub fn split_halfedge(
        &mut self,
        vert: &manifold_mesh3d::Vertex,
        ind_halfedge: usize,
    ) -> Result<()> {
        let ind_near_vert = self.mesh.get_halfedge(ind_halfedge)?.first_vertex().ind();
        let (it, _) = self.vertex_edges[ind_near_vert][0];
        let ind_tet = self
            .del_struct
            .get_simplicial()
            .get_halftriangle(it)?
            .tetrahedron()
            .ind();

        log::debug!("he to split {}", ind_halfedge);
        let ind_vertex = mesh_operations::split_halfedge(self.mesh, vert, ind_halfedge)?;
        log::debug!(
            "Added {} , ({}, {}, {})",
            ind_vertex,
            vert[0],
            vert[1],
            vert[2]
        );
        self.insert_vertex(ind_vertex, ind_tet)
    }

    /// Splits given face
    pub fn split_face(&mut self, vert: &manifold_mesh3d::Vertex, ind_face: usize) -> Result<()> {
        let ind_near_vert = self.mesh.get_face(ind_face)?.vertices_inds()[0];
        let (it, _) = self.vertex_edges[ind_near_vert][0];
        let ind_tet = self
            .del_struct
            .get_simplicial()
            .get_halftriangle(it)?
            .tetrahedron()
            .ind();
        let ind_vertex = mesh_operations::split_face(self.mesh, vert, ind_face)?;
        self.insert_vertex(ind_vertex, ind_tet)
    }
}
