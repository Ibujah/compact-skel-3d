use anyhow::Result;
use nalgebra::base::*;

/// Mesh vertex
pub type Vertex = Vector3<f32>;
/// Mesh edge
pub type Edge = [usize; 2];
/// Mesh face (array of vertex indices)
pub type Face = [usize; 3];

#[derive(Clone)]
/// Generic non manifold Mesh
pub struct GenericMesh3D {
    pub(super) vertices: Vec<Vertex>,
    pub(super) edges: Vec<Edge>,
    pub(super) faces: Vec<Face>,

    pub(super) map_vert_edg: Vec<Vec<usize>>,
    pub(super) map_edg_face: Vec<Vec<usize>>,
}

impl GenericMesh3D {
    /// Generic mesh constructor
    pub fn new() -> GenericMesh3D {
        GenericMesh3D {
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: Vec::new(),

            map_vert_edg: Vec::new(),
            map_edg_face: Vec::new(),
        }
    }

    /// Adds a vertex to the mesh
    pub fn add_vertex(&mut self, point: &Vector3<f32>) -> usize {
        self.vertices.push(*point);
        self.map_vert_edg.push(Vec::new());
        self.vertices.len() - 1
    }

    fn get_vertex_uncheck(&self, ind_vertex: usize) -> Vertex {
        self.vertices[ind_vertex]
    }

    /// Vertex getter
    pub fn get_vertex(&self, ind_vertex: usize) -> Result<Vertex> {
        if ind_vertex >= self.vertices.len() {
            return Err(anyhow::Error::msg("get_vertex(): Index out of bounds"));
        }

        Ok(self.get_vertex_uncheck(ind_vertex))
    }

    /// Gets number of vertices
    pub fn get_nb_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Adds an edge to the mesh
    pub fn add_edge(&mut self, ind_vertex1: usize, ind_vertex2: usize) -> Result<usize> {
        let edge = if ind_vertex1 < ind_vertex2 {
            [ind_vertex1, ind_vertex2]
        } else {
            [ind_vertex2, ind_vertex1]
        };
        let [ind_vertex1, ind_vertex2] = edge;
        if ind_vertex2 >= self.map_vert_edg.len() {
            return Err(anyhow::Error::msg("add_edge(): Vertex index out of bounds"));
        }

        let ind_edge = self.map_vert_edg[ind_vertex1]
            .iter()
            .fold(None, |res, &ind_edg| {
                if self.edges[ind_edg] == edge {
                    Some(ind_edg)
                } else {
                    res
                }
            });

        match ind_edge {
            Some(ind) => Ok(ind),
            None => {
                self.edges.push(edge);
                self.map_edg_face.push(Vec::new());
                self.map_vert_edg[ind_vertex1].push(self.edges.len() - 1);
                self.map_vert_edg[ind_vertex2].push(self.edges.len() - 1);
                Ok(self.edges.len() - 1)
            }
        }
    }

    fn get_edge_uncheck(&self, ind_edge: usize) -> Edge {
        self.edges[ind_edge]
    }

    /// Edge getter
    pub fn get_edge(&self, ind_edge: usize) -> Result<Edge> {
        if ind_edge >= self.edges.len() {
            return Err(anyhow::Error::msg("get_edge(): Index out of bounds"));
        }
        Ok(self.get_edge_uncheck(ind_edge))
    }

    /// Get number of edges
    pub fn get_nb_edges(&self) -> usize {
        self.edges.len()
    }

    /// Adds a face to the mesh
    pub fn add_face(
        &mut self,
        ind_vertex1: usize,
        ind_vertex2: usize,
        ind_vertex3: usize,
    ) -> Result<usize> {
        let mut face = [ind_vertex1, ind_vertex2, ind_vertex3];
        face.sort();
        let [ind_vertex1, ind_vertex2, ind_vertex3] = face;
        if ind_vertex3 >= self.map_vert_edg.len() {
            return Err(anyhow::Error::msg("add_face(): Vertex index out of bounds"));
        }
        let ind_edge1 = self.add_edge(ind_vertex1, ind_vertex2)?;

        let opt_fac = self.map_edg_face[ind_edge1]
            .iter()
            .find(|&ind_face| self.faces[*ind_face] == face);
        if let Some(&ind_fac) = opt_fac {
            return Ok(ind_fac);
        }

        let ind_edge2 = self.add_edge(ind_vertex2, ind_vertex3)?;
        let ind_edge3 = self.add_edge(ind_vertex1, ind_vertex3)?;

        self.faces.push([ind_vertex1, ind_vertex2, ind_vertex3]);
        let ind_face = self.faces.len() - 1;

        self.map_edg_face[ind_edge1].push(ind_face);
        self.map_edg_face[ind_edge2].push(ind_face);
        self.map_edg_face[ind_edge3].push(ind_face);

        Ok(ind_face)
    }

    fn get_face_uncheck(&self, ind_face: usize) -> Face {
        self.faces[ind_face]
    }

    /// Face getter
    pub fn get_face(&self, ind_face: usize) -> Result<Face> {
        if ind_face >= self.faces.len() {
            return Err(anyhow::Error::msg("get_face(): Index out of bounds"));
        }
        Ok(self.get_face_uncheck(ind_face))
    }

    /// Gets number of faces
    pub fn get_nb_faces(&self) -> usize {
        self.faces.len()
    }

    /// Checks if an edge is in the mesh
    ///
    /// Returns edge index if found
    pub fn is_edge_in(&self, ind_vertex1: usize, ind_vertex2: usize) -> Option<usize> {
        let edge = if ind_vertex1 < ind_vertex2 {
            [ind_vertex1, ind_vertex2]
        } else {
            [ind_vertex2, ind_vertex1]
        };
        let [ind_vertex1, ind_vertex2] = edge;
        if ind_vertex2 >= self.vertices.len() {
            return None;
        } else {
            let opt_ind_edg = self.map_vert_edg[ind_vertex1]
                .iter()
                .find(|&&ind_edg| self.edges[ind_edg] == edge);
            if let Some(&ind_edg) = opt_ind_edg {
                return Some(ind_edg);
            }
        }
        None
    }

    /// Checks if a face is in the mesh
    ///
    /// Returns face index if found
    pub fn is_face_in(
        &self,
        ind_vertex1: usize,
        ind_vertex2: usize,
        ind_vertex3: usize,
    ) -> Option<usize> {
        let mut face = [ind_vertex1, ind_vertex2, ind_vertex3];
        face.sort();
        let [ind_vertex1, ind_vertex2, ind_vertex3] = face;
        if ind_vertex3 >= self.vertices.len() {
            return None;
        } else {
            let opt_ind_edg = self.is_edge_in(ind_vertex1, ind_vertex2);
            if let Some(ind_edg) = opt_ind_edg {
                let opt_fac = self.map_edg_face[ind_edg]
                    .iter()
                    .find(|&&ind_face| self.faces[ind_face] == face);
                if let Some(&ind_fac) = opt_fac {
                    return Some(ind_fac);
                }
            }
        }
        None
    }
}
