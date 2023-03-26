use anyhow::Result;
use std::collections::{HashMap, HashSet};
use tritet::{StrError, Tetgen};

use crate::mesh3d::mesh_operations;
use crate::mesh3d::{manifold_mesh3d, ManifoldMesh3D};

pub type Edge = [usize; 2];
pub type Triangle = [usize; 3];
pub type Tetrahedron = [usize; 4];

/// Delaunay structure
pub struct DelaunayInterface<'a> {
    mesh: &'a mut ManifoldMesh3D,

    edges: HashSet<Edge>,
    faces: HashMap<Triangle, Vec<Tetrahedron>>,
    tetras: HashSet<Tetrahedron>,

    initial_vertices_number: usize,
}

fn to_anyhow(err: StrError) -> anyhow::Error {
    anyhow::Error::msg(err.to_string())
}

impl<'a> DelaunayInterface<'a> {
    fn insert_tetra(&mut self, tetra: &mut Tetrahedron) -> () {
        tetra.sort();

        self.edges.insert([tetra[0], tetra[1]]);
        self.edges.insert([tetra[0], tetra[2]]);
        self.edges.insert([tetra[0], tetra[3]]);
        self.edges.insert([tetra[1], tetra[2]]);
        self.edges.insert([tetra[1], tetra[3]]);
        self.edges.insert([tetra[2], tetra[3]]);

        self.faces
            .entry([tetra[0], tetra[1], tetra[2]])
            .or_insert(Vec::new())
            .push(*tetra);
        self.faces
            .entry([tetra[0], tetra[1], tetra[3]])
            .or_insert(Vec::new())
            .push(*tetra);
        self.faces
            .entry([tetra[0], tetra[2], tetra[3]])
            .or_insert(Vec::new())
            .push(*tetra);
        self.faces
            .entry([tetra[1], tetra[2], tetra[3]])
            .or_insert(Vec::new())
            .push(*tetra);

        self.tetras.insert([tetra[0], tetra[1], tetra[2], tetra[3]]);
    }

    fn generate_struct(&mut self) -> Result<()> {
        let mut tetgen = Tetgen::new(
            self.mesh.get_nb_vertices(),
            Some(vec![3; self.mesh.get_nb_faces()]),
            None,
            None,
        )
        .map_err(to_anyhow)?;

        for v in self.mesh.vertex_indices() {
            // let mut i = 0;
            // for (&v, vert) in self.mesh.vertices() {
            //     println!("{}:{} ({}, {}, {})", i, v, vert[0], vert[1], vert[2]);
            //     i = i + 1;
            let vert = self.mesh.get_vertex(v)?.vertex();
            if v >= self.mesh.get_nb_vertices() {
                return Err(anyhow::Error::msg(
                    "generate_struct(): Vertex index over vertex number, not currently handled",
                ));
            }
            tetgen
                .set_point(v, vert[0] as f64, vert[1] as f64, vert[2] as f64)
                .map_err(to_anyhow)?;
        }

        tetgen.generate_delaunay(false).map_err(to_anyhow)?;

        for t in 0..tetgen.ntet() {
            let mut tetra = [0; 4];
            for m in 0..4 {
                tetra[m] = tetgen.tet_node(t, m);
            }

            self.insert_tetra(&mut tetra);
        }

        Ok(())
    }

    fn recompute_struct(&mut self) -> Result<()> {
        self.edges = HashSet::new();
        self.faces = HashMap::new();
        self.tetras = HashSet::new();

        self.generate_struct()
    }

    /// Creates Delaunay structure from mesh
    pub fn from_mesh(mesh: &'a mut ManifoldMesh3D) -> Result<DelaunayInterface<'a>> {
        let initial_vertices_number = mesh.get_nb_vertices();
        let mut deltet = DelaunayInterface {
            mesh,
            edges: HashSet::new(),
            faces: HashMap::new(),
            tetras: HashSet::new(),
            initial_vertices_number,
        };

        deltet.generate_struct()?;

        Ok(deltet)
    }

    /// Mesh getter
    pub fn get_mesh(&self) -> &ManifoldMesh3D {
        self.mesh
    }

    /// Edge set getter
    pub fn get_edges(&self) -> &HashSet<Edge> {
        &self.edges
    }

    /// Face map getter
    pub fn get_faces(&self) -> &HashMap<Triangle, Vec<Tetrahedron>> {
        &self.faces
    }

    /// Tetrahedra set getter
    pub fn get_tetrahedra(&self) -> &HashSet<Tetrahedron> {
        &self.tetras
    }

    /// Gets tetrahedra surrounding a given triangle
    pub fn get_tetrahedra_from_triangle(&self, del_tri: Triangle) -> Result<Vec<Tetrahedron>> {
        let vec = self
            .faces
            .get(&del_tri)
            .ok_or(anyhow::Error::msg("Triangle does not exist"))?
            .iter()
            .map(|&x| x)
            .collect();
        Ok(vec)
    }

    /// Checks if vertex was an original mesh vertex
    pub fn is_original_vertex(&self, ind_vertex: usize) -> bool {
        ind_vertex < self.initial_vertices_number
    }

    /// Checks if edge is in Delaunay
    pub fn is_edge_in(&self, edge: &Edge) -> bool {
        let mut edge_sort = [edge[0], edge[1]];
        edge_sort.sort();
        self.edges.contains(&edge_sort)
    }

    /// Checks if face is in Delaunay
    pub fn is_face_in(&self, face: &Triangle) -> bool {
        let mut face_sort = [face[0], face[1], face[2]];
        face_sort.sort();
        self.faces.contains_key(&face_sort)
    }

    /// Checks if tetrahedron is in Delaunay
    pub fn is_tetra_in(&self, tetra: &Tetrahedron) -> bool {
        let mut tetra_sort = [tetra[0], tetra[1], tetra[2], tetra[3]];
        tetra_sort.sort();
        self.tetras.contains(&tetra_sort)
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
        mesh_operations::split_halfedge(self.mesh, vert, ind_halfedge)?;
        self.recompute_struct()
    }

    /// Splits given face
    pub fn split_face(&mut self, vert: &manifold_mesh3d::Vertex, ind_face: usize) -> Result<()> {
        mesh_operations::split_face(self.mesh, vert, ind_face)?;
        self.recompute_struct()
    }
}
