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

    initial_vertices_number: usize,
}

impl<'a> DelaunayInterface<'a> {
    fn generate_struct(&mut self) -> Result<()> {
        let mut points = Vec::new();
        for v in self.mesh.vertex_indices() {
            let vert = self.mesh.get_vertex(v)?.vertex();
            points.push([vert[0] as f64, vert[1] as f64, vert[2] as f64]);
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

    fn insert_vertex(&mut self, ind_vertex: usize) -> Result<()> {
        let vert = self.mesh.get_vertex(ind_vertex)?.vertex();
        self.del_struct
            .insert_vertex([vert[0] as f64, vert[1] as f64, vert[2] as f64], None)?;

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
                    if self
                        .del_struct
                        .get_simplicial()
                        .get_halftriangle(it)
                        .unwrap()
                        .halfedges()[i]
                        .first_node()
                        .equals(&Node::Value(iv))
                    {
                        Some((it, i))
                    } else {
                        None
                    }
                })
                .collect();
        }

        Ok(())
    }

    /// Creates Delaunay structure from mesh
    pub fn from_mesh(mesh: &'a mut ManifoldMesh3D) -> Result<DelaunayInterface<'a>> {
        let initial_vertices_number = mesh.get_nb_vertices();
        let mut deltet = DelaunayInterface {
            mesh,
            del_struct: DelaunayStructure3D::new(),
            vertex_edges: Vec::new(),
            initial_vertices_number,
        };

        deltet.generate_struct()?;

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
    pub fn is_edge_in(&self, edge: &Edge) -> bool {
        self.vertex_edges[edge[0]]
            .iter()
            .position(|&(it, i)| {
                self.del_struct
                    .get_simplicial()
                    .get_halftriangle(it)
                    .unwrap()
                    .halfedges()[i]
                    .last_node()
                    .equals(&Node::Value(edge[1]))
            })
            .is_some()
    }

    /// Checks if face is in Delaunay
    pub fn is_face_in(&self, face: &Triangle) -> bool {
        self.vertex_edges[face[0]]
            .iter()
            .position(|&(it, i)| {
                let he = self
                    .del_struct
                    .get_simplicial()
                    .get_halftriangle(it)
                    .unwrap()
                    .halfedges()[i];
                if he.last_node().equals(&Node::Value(face[1])) {
                    he.next().last_node().equals(&Node::Value(face[2]))
                } else {
                    false
                }
            })
            .is_some()
    }

    /// Count number of non Delaunay halfedges
    pub fn count_non_del_halfedges(&self) -> Result<usize> {
        let mut nb_non_del = 0;
        for (_, &he) in self.mesh.halfedges() {
            // let he = self.mesh.get_halfedge(i)?.halfedge();
            nb_non_del = nb_non_del + if self.is_edge_in(&he) { 0 } else { 1 };
        }
        Ok(nb_non_del)
    }

    /// Count number of non Delaunay faces
    pub fn count_non_del_faces(&self) -> Result<usize> {
        let mut nb_non_del = 0;
        for (&i, _) in self.mesh.faces() {
            let face = self.mesh.get_face(i).unwrap().vertices_inds();
            nb_non_del = nb_non_del + if self.is_face_in(&face) { 0 } else { 1 };
        }
        Ok(nb_non_del)
    }

    fn get_opposite_angle(&self, halfedge: manifold_mesh3d::IterHalfEdge) -> Result<f32> {
        let vert1 = halfedge.first_vertex().vertex();
        let vert2 = halfedge.last_vertex().vertex();
        let vert3 = halfedge
            .next_halfedge()
            .ok_or(anyhow::Error::msg("get_opposite_angle(): No next halfedge"))?
            .last_vertex()
            .vertex();

        let vec31 = vert1 - vert3;
        let vec32 = vert2 - vert3;

        let vx = vec31.dot(&vec32);
        let vy = vec31.cross(&vec32).norm();
        let angle = vy.atan2(vx);
        Ok(angle)
    }

    /// Gets first locally non Delaunay halfedge, starting from a shift
    pub fn get_local_non_del_halfedge(
        &self,
        shift: Option<usize>,
    ) -> Result<Option<manifold_mesh3d::IterHalfEdge>> {
        let shift = shift.unwrap_or(0) % self.mesh.get_nb_halfedges();
        let mut i = 0;
        for (&ind_he, seg) in self.mesh.halfedges() {
            if i <= shift {
                i = i + 1;
                continue;
            }
            if !self.is_edge_in(&seg) {
                let he = self.mesh.get_halfedge(ind_he)?;
                let angle1 = self.get_opposite_angle(he)?;
                let angle2 = self.get_opposite_angle(he.opposite_halfedge().ok_or(
                    anyhow::Error::msg("get_opposite_angle(): No opposite halfedge"),
                )?)?;

                if angle1 + angle2 >= std::f32::consts::PI {
                    return Ok(Some(he));
                }
            };
        }
        Ok(None)
    }

    /// Gets first globally non Delaunay halfedge, starting from a shift
    pub fn get_non_del_halfedge(
        &self,
        shift: Option<usize>,
    ) -> Result<Option<manifold_mesh3d::IterHalfEdge>> {
        let shift = shift.unwrap_or(0) % self.mesh.get_nb_halfedges();
        let mut i = 0;
        for (&ind_he, seg) in self.mesh.halfedges() {
            if i <= shift {
                i = i + 1;
                continue;
            }
            if !self.is_edge_in(&seg) {
                let he = self.mesh.get_halfedge(ind_he)?;
                return Ok(Some(he));
            };
        }
        Ok(None)
    }

    /// Gets first globally non Delaunay face, starting from a shift
    pub fn get_non_del_face(
        &self,
        shift: Option<usize>,
    ) -> Result<Option<manifold_mesh3d::IterFace>> {
        let shift = shift.unwrap_or(0) % self.mesh.get_nb_halfedges();
        let mut i = 0;
        for (&ind_face, _) in self.mesh.faces() {
            if i <= shift {
                i = i + 1;
                continue;
            }
            let face = self.mesh.get_face(ind_face)?;
            let face_ind = face.vertices_inds();
            if !self.is_face_in(&face_ind) {
                return Ok(Some(face));
            };
        }
        Ok(None)
    }

    /// Gets all non delaunay halfedges
    pub fn get_all_non_del_halfedge(&self) -> Result<Vec<usize>> {
        let mut non_del = Vec::new();

        for (&ind_he, seg) in self.mesh.halfedges() {
            if !self.is_edge_in(&seg) {
                non_del.push(ind_he);
            };
        }

        Ok(non_del)
    }

    /// Gets all non delaunay faces
    pub fn get_all_non_del_face(&self) -> Result<Vec<usize>> {
        let mut non_del = Vec::new();

        for (&ind_fac, _) in self.mesh.faces() {
            let face = self.mesh.get_face(ind_fac).unwrap().vertices_inds();
            if !self.is_face_in(&face) {
                non_del.push(ind_fac);
            };
        }

        Ok(non_del)
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
        let ind_vertex = mesh_operations::split_halfedge(self.mesh, vert, ind_halfedge)?;
        self.insert_vertex(ind_vertex)
    }

    /// Splits given face
    pub fn split_face(&mut self, vert: &manifold_mesh3d::Vertex, ind_face: usize) -> Result<()> {
        let ind_vertex = mesh_operations::split_face(self.mesh, vert, ind_face)?;
        self.insert_vertex(ind_vertex)
    }
}
