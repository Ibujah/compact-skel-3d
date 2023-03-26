use anyhow::Result;
use nalgebra::base::*;
use std::collections::HashMap;

/// Mesh vertex
pub type Vertex = Vector3<f32>;
/// Mesh halfedge
pub type HalfEdge = [usize; 2];
/// Mesh face (array of halfedges)
pub type FaceHalfedges = [usize; 3];

#[derive(Clone)]
/// Manifold mesh
pub struct ManifoldMesh3D {
    pub(super) vertices: HashMap<usize, Vertex>,
    pub(super) halfedges: HashMap<usize, HalfEdge>,
    pub(super) faces: HashMap<usize, FaceHalfedges>,
    pub(super) last_ind_vert: usize,
    pub(super) last_ind_hedge: usize,
    pub(super) last_ind_face: usize,

    pub(super) map_vert_hedg: HashMap<usize, Vec<usize>>,
    pub(super) map_hedg_face: HashMap<usize, usize>,
    pub(super) map_hedg_opp: HashMap<usize, usize>,
    pub(super) map_hedg_next: HashMap<usize, usize>,
    pub(super) map_hedg_prev: HashMap<usize, usize>,
}

#[derive(Copy, Clone)]
/// Vertex iterator
pub struct IterVertex<'a> {
    mesh: &'a ManifoldMesh3D,
    ind_vertex: usize,
}

#[derive(Copy, Clone)]
/// Halfedge iterator
pub struct IterHalfEdge<'a> {
    mesh: &'a ManifoldMesh3D,
    ind_halfedge: usize,
}

#[derive(Copy, Clone)]
/// Face iterator
pub struct IterFace<'a> {
    mesh: &'a ManifoldMesh3D,
    ind_face: usize,
}

impl ManifoldMesh3D {
    /// Manifold mesh constructor
    pub fn new() -> ManifoldMesh3D {
        ManifoldMesh3D {
            vertices: HashMap::new(),
            halfedges: HashMap::new(),
            faces: HashMap::new(),
            last_ind_vert: 0,
            last_ind_hedge: 0,
            last_ind_face: 0,

            map_vert_hedg: HashMap::new(),
            map_hedg_face: HashMap::new(),
            map_hedg_opp: HashMap::new(),
            map_hedg_next: HashMap::new(),
            map_hedg_prev: HashMap::new(),
        }
    }

    /// Adds a vertex to th mesh
    pub fn add_vertex(&mut self, point: &Vector3<f32>) -> usize {
        self.vertices.insert(self.last_ind_vert, *point);
        self.map_vert_hedg.insert(self.last_ind_vert, Vec::new());
        self.last_ind_vert = self.last_ind_vert + 1;
        self.last_ind_vert - 1
    }

    fn get_vertex_uncheck(&self, ind_vertex: usize) -> IterVertex {
        IterVertex {
            mesh: self,
            ind_vertex,
        }
    }

    /// Vertex getter
    pub fn get_vertex(&self, ind_vertex: usize) -> Result<IterVertex> {
        if !self.vertices.contains_key(&ind_vertex) {
            return Err(anyhow::Error::msg("get_vertex(): Index out of bounds"));
        }

        Ok(self.get_vertex_uncheck(ind_vertex))
    }

    /// Gets number of vertices
    pub fn get_nb_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Gets list of vertex index
    pub fn vertex_indices(&self) -> Vec<usize> {
        let mut inds: Vec<usize> = self.vertices.iter().map(|(&ind, _)| ind).collect();
        inds.sort();
        inds
    }

    /// Gets vertex map
    pub fn vertices(&self) -> &HashMap<usize, Vertex> {
        &self.vertices
    }

    fn add_halfedge_uncheck(&mut self, ind_vertex1: usize, ind_vertex2: usize) -> usize {
        self.halfedges
            .insert(self.last_ind_hedge, [ind_vertex1, ind_vertex2]);
        self.map_vert_hedg
            .get_mut(&ind_vertex1)
            .unwrap()
            .push(self.last_ind_hedge);

        if let Some(&ind_opp) = self
            .map_vert_hedg
            .get(&ind_vertex2)
            .unwrap()
            .iter()
            .find(|ind_he| self.halfedges.get(ind_he).unwrap()[1] == ind_vertex1)
        {
            self.map_hedg_opp.insert(ind_opp, self.last_ind_hedge);
            self.map_hedg_opp.insert(self.last_ind_hedge, ind_opp);
        }

        self.last_ind_hedge = self.last_ind_hedge + 1;
        self.last_ind_hedge - 1
    }

    fn get_halfedge_uncheck(&self, ind_halfedge: usize) -> IterHalfEdge {
        IterHalfEdge {
            mesh: self,
            ind_halfedge,
        }
    }

    /// Halfedge getter
    pub fn get_halfedge(&self, ind_halfedge: usize) -> Result<IterHalfEdge> {
        if !self.halfedges.contains_key(&ind_halfedge) {
            return Err(anyhow::Error::msg("get_halfedge(): Index out of bounds"));
        }
        Ok(self.get_halfedge_uncheck(ind_halfedge))
    }

    /// Gets number of halfedges
    pub fn get_nb_halfedges(&self) -> usize {
        self.halfedges.len()
    }

    /// Gets halfedge map
    pub fn halfedges(&self) -> &HashMap<usize, HalfEdge> {
        &self.halfedges
    }

    /// Adds a face and associated halfedges
    pub fn add_face(
        &mut self,
        ind_vertex1: usize,
        ind_vertex2: usize,
        ind_vertex3: usize,
    ) -> Result<usize> {
        if !self.vertices.contains_key(&ind_vertex1)
            || !self.vertices.contains_key(&ind_vertex2)
            || !self.vertices.contains_key(&ind_vertex3)
        {
            return Err(anyhow::Error::msg("add_face(): Index out of bounds"));
        }

        if self
            .map_vert_hedg
            .get(&ind_vertex1)
            .unwrap()
            .iter()
            .find(|ind_he| self.halfedges.get(ind_he).unwrap()[1] == ind_vertex2)
            .is_some()
        {
            return Err(anyhow::Error::msg("add_face(): halfedge already exists"));
        }

        if self
            .map_vert_hedg
            .get(&ind_vertex2)
            .unwrap()
            .iter()
            .find(|ind_he| self.halfedges.get(ind_he).unwrap()[1] == ind_vertex3)
            .is_some()
        {
            return Err(anyhow::Error::msg("add_face(): halfedge already exists"));
        }

        if self
            .map_vert_hedg
            .get(&ind_vertex3)
            .unwrap()
            .iter()
            .find(|ind_he| self.halfedges.get(ind_he).unwrap()[1] == ind_vertex1)
            .is_some()
        {
            return Err(anyhow::Error::msg("add_face(): halfedge already exists"));
        }

        let ind_halfedge1 = self.add_halfedge_uncheck(ind_vertex1, ind_vertex2);
        let ind_halfedge2 = self.add_halfedge_uncheck(ind_vertex2, ind_vertex3);
        let ind_halfedge3 = self.add_halfedge_uncheck(ind_vertex3, ind_vertex1);

        self.faces.insert(
            self.last_ind_face,
            [ind_halfedge1, ind_halfedge2, ind_halfedge3],
        );

        self.map_hedg_face.insert(ind_halfedge1, self.last_ind_face);
        self.map_hedg_face.insert(ind_halfedge2, self.last_ind_face);
        self.map_hedg_face.insert(ind_halfedge3, self.last_ind_face);

        self.map_hedg_next.insert(ind_halfedge1, ind_halfedge2);
        self.map_hedg_next.insert(ind_halfedge2, ind_halfedge3);
        self.map_hedg_next.insert(ind_halfedge3, ind_halfedge1);

        self.map_hedg_prev.insert(ind_halfedge1, ind_halfedge3);
        self.map_hedg_prev.insert(ind_halfedge2, ind_halfedge1);
        self.map_hedg_prev.insert(ind_halfedge3, ind_halfedge2);

        self.last_ind_face = self.last_ind_face + 1;
        Ok(self.last_ind_face - 1)
    }

    /// Removes a face and associated halfedges
    pub fn remove_face(&mut self, ind_face: usize) -> Result<()> {
        let [ind_he1, ind_he2, ind_he3] = self
            .faces
            .remove(&ind_face)
            .ok_or(anyhow::Error::msg("remove_face() face does not exist"))?;

        self.map_hedg_face.remove(&ind_he1);
        self.map_hedg_face.remove(&ind_he2);
        self.map_hedg_face.remove(&ind_he3);

        self.map_hedg_next.remove(&ind_he1);
        self.map_hedg_next.remove(&ind_he2);
        self.map_hedg_next.remove(&ind_he3);

        self.map_hedg_prev.remove(&ind_he1);
        self.map_hedg_prev.remove(&ind_he2);
        self.map_hedg_prev.remove(&ind_he3);

        if let Some(ind_he1_opp) = self.map_hedg_opp.remove(&ind_he1) {
            self.map_hedg_opp.remove(&ind_he1_opp).unwrap();
        }
        if let Some(ind_he2_opp) = self.map_hedg_opp.remove(&ind_he2) {
            self.map_hedg_opp.remove(&ind_he2_opp).unwrap();
        }
        if let Some(ind_he3_opp) = self.map_hedg_opp.remove(&ind_he3) {
            self.map_hedg_opp.remove(&ind_he3_opp).unwrap();
        }

        let [ind_v1, _] = self.halfedges.remove(&ind_he1).unwrap();
        let [ind_v2, _] = self.halfedges.remove(&ind_he2).unwrap();
        let [ind_v3, _] = self.halfedges.remove(&ind_he3).unwrap();

        self.map_vert_hedg
            .get_mut(&ind_v1)
            .unwrap()
            .retain(|&ind| ind != ind_he1);
        self.map_vert_hedg
            .get_mut(&ind_v2)
            .unwrap()
            .retain(|&ind| ind != ind_he2);
        self.map_vert_hedg
            .get_mut(&ind_v3)
            .unwrap()
            .retain(|&ind| ind != ind_he3);

        Ok(())
    }

    fn get_face_uncheck(&self, ind_face: usize) -> IterFace {
        IterFace {
            mesh: self,
            ind_face,
        }
    }

    /// Face getter
    pub fn get_face(&self, ind_face: usize) -> Result<IterFace> {
        if !self.faces.contains_key(&ind_face) {
            return Err(anyhow::Error::msg("get_face(): Index out of bounds"));
        }
        Ok(self.get_face_uncheck(ind_face))
    }

    /// gets number of faces
    pub fn get_nb_faces(&self) -> usize {
        self.faces.len()
    }

    /// Gets face map
    pub fn faces(&self) -> &HashMap<usize, FaceHalfedges> {
        &self.faces
    }

    /// Checks if an edge is in the mesh
    ///
    /// Returns halfedge iterator if found
    pub fn is_edge_in(&self, ind_vertex1: usize, ind_vertex2: usize) -> Option<IterHalfEdge> {
        if !self.vertices.contains_key(&ind_vertex1) || !self.vertices.contains_key(&ind_vertex2) {
            return None;
        } else {
            let vertex1 = self.get_vertex_uncheck(ind_vertex1);
            for he in vertex1.halfedges() {
                if he.last_vertex().ind() == ind_vertex2 {
                    return Some(he);
                }
            }
        }
        None
    }

    /// Checks if a face is in the mesh
    ///
    /// Returns face iterator if found
    pub fn is_face_in(
        &self,
        ind_vertex1: usize,
        ind_vertex2: usize,
        ind_vertex3: usize,
    ) -> Option<IterFace> {
        if !self.vertices.contains_key(&ind_vertex3) {
            return None;
        } else {
            let opt_he = self.is_edge_in(ind_vertex1, ind_vertex2);
            if let Some(he) = opt_he {
                if let Some(he_next) = he.next_halfedge() {
                    if he_next.last_vertex().ind() == ind_vertex3 {
                        if let Some(face) = he.face() {
                            return Some(face);
                        }
                    }
                }
                if let Some(he_opp) = he.opposite_halfedge() {
                    if let Some(he_next) = he_opp.next_halfedge() {
                        if he_next.last_vertex().ind() == ind_vertex3 {
                            if let Some(face) = he_next.face() {
                                return Some(face);
                            }
                        }
                    }
                }
            }
        }

        None
    }

    fn check_face(&self, ind_face: usize) -> Result<()> {
        let face = self.get_face(ind_face)?;
        // check edges existence
        for hedg in face.halfedges() {
            let face_comp = hedg.face().ok_or(anyhow::Error::msg(
                "check_face(): Halfedge should be linked to a face",
            ))?;
            if face_comp.ind() != face.ind() {
                return Err(anyhow::Error::msg(
                    "check_face(): HalfEdge linked to wrong face",
                ));
            }
        }
        Ok(())
    }

    fn check_halfedge(&self, ind_hedge: usize) -> Result<()> {
        let halfedge = self.get_halfedge(ind_hedge)?;

        let face = halfedge.face().ok_or(anyhow::Error::msg(
            "check_halfedge(): Halfedge should be linked to a face",
        ))?;

        let halfedge_next = halfedge.next_halfedge().ok_or(anyhow::Error::msg(
            "check_halfedge(): Halfedge should have next halfedge",
        ))?;

        let halfedge_prev = halfedge.prev_halfedge().ok_or(anyhow::Error::msg(
            "check_halfedge(): Halfedge should have previous halfedge",
        ))?;

        let face_next = halfedge_next.face().ok_or(anyhow::Error::msg(
            "check_halfedge(): Next halfedge should be linked to face",
        ))?;

        let face_prev = halfedge_prev.face().ok_or(anyhow::Error::msg(
            "check_halfedge(): Previous halfedge should be linked to face",
        ))?;

        if halfedge.last_vertex().ind() != halfedge_next.first_vertex().ind() {
            return Err(anyhow::Error::msg(
                "check_halfedge(): Next halfedge not starting with last vertex",
            ));
        }
        if halfedge.first_vertex().ind() != halfedge_prev.last_vertex().ind() {
            return Err(anyhow::Error::msg(
                "check_halfedge(): Previous halfedge not ending with first vertex",
            ));
        }

        if face.ind() != face_next.ind() {
            return Err(anyhow::Error::msg(
                "check_halfedge(): Next halfedge not on the same face",
            ));
        }
        if face.ind() != face_prev.ind() {
            return Err(anyhow::Error::msg(
                "check_halfedge(): Previous halfedge not on the same face",
            ));
        }

        // check opposite
        let halfedge_opp = halfedge.opposite_halfedge().ok_or(anyhow::Error::msg(
            "check_halfedge(): Halfedge should have opposite halfedge",
        ))?;

        if halfedge.first_vertex().ind() != halfedge_opp.last_vertex().ind() {
            return Err(anyhow::Error::msg(
                "check_halfedge(): Opposite halfedge not starting with last vertex",
            ));
        }
        if halfedge.last_vertex().ind() != halfedge_opp.first_vertex().ind() {
            return Err(anyhow::Error::msg(
                "check_halfedge(): Opposite halfedge not ending with first vertex",
            ));
        }

        // check vertices
        let neigh_hedges = halfedge.first_vertex().halfedges();

        let is_in = neigh_hedges.iter().fold(false, |res, &iterhedge| {
            res || iterhedge.ind() == halfedge.ind()
        });

        if !is_in {
            return Err(anyhow::Error::msg(
                "check_halfedge(): Halfedge not in vertex",
            ));
        }

        Ok(())
    }

    fn check_vertex(&self, ind_vertex: usize) -> Result<()> {
        let vertex = self.get_vertex(ind_vertex)?;

        for he in vertex.halfedges() {
            if he.first_vertex().ind() != ind_vertex {
                return Err(anyhow::Error::msg(
                    "check_vertex(): Vertex contains non coherent halfedge",
                ));
            }
        }

        Ok(())
    }

    /// Checks integrity of the mesh
    pub fn check_mesh(&self) -> Result<()> {
        for (&f, _) in self.faces.iter() {
            self.check_face(f)?;
        }

        for (&e, _) in self.halfedges.iter() {
            self.check_halfedge(e)?;
        }

        for (&v, _) in self.vertices.iter() {
            self.check_vertex(v)?;
        }

        Ok(())
    }
}

impl<'a> IterVertex<'a> {
    /// Gets vertex coordinates
    pub fn vertex(&self) -> Vertex {
        *self.mesh.vertices.get(&self.ind_vertex).unwrap()
    }

    /// Gets vertex index
    pub fn ind(&self) -> usize {
        self.ind_vertex
    }

    /// Gets list of halfedges starting at this vertex
    pub fn halfedges(&self) -> Vec<IterHalfEdge<'a>> {
        let vec_he = self
            .mesh
            .map_vert_hedg
            .get(&self.ind_vertex)
            .unwrap()
            .iter()
            .fold(Vec::new(), |mut v, &x| {
                v.push(IterHalfEdge {
                    mesh: self.mesh,
                    ind_halfedge: x,
                });
                v
            });
        vec_he
    }
}

impl<'a> IterHalfEdge<'a> {
    /// Gets halfedge (array of vertex indices)
    pub fn halfedge(&self) -> HalfEdge {
        *self.mesh.halfedges.get(&self.ind_halfedge).unwrap()
    }

    /// Gets halfedge index
    pub fn ind(&self) -> usize {
        self.ind_halfedge
    }

    /// First vertex iterator
    pub fn first_vertex(&self) -> IterVertex<'a> {
        IterVertex {
            mesh: self.mesh,
            ind_vertex: self.halfedge()[0],
        }
    }

    /// Last vertex iterator
    pub fn last_vertex(&self) -> IterVertex<'a> {
        IterVertex {
            mesh: self.mesh,
            ind_vertex: self.halfedge()[1],
        }
    }

    /// Next halfedge on same face
    pub fn next_halfedge(&self) -> Option<IterHalfEdge<'a>> {
        if let Some(&ind_next) = self.mesh.map_hedg_next.get(&self.ind_halfedge) {
            Some(IterHalfEdge {
                mesh: self.mesh,
                ind_halfedge: ind_next,
            })
        } else {
            None
        }
    }

    /// Previous halfedge on same face
    pub fn prev_halfedge(&self) -> Option<IterHalfEdge<'a>> {
        if let Some(&ind_prev) = self.mesh.map_hedg_prev.get(&self.ind_halfedge) {
            Some(IterHalfEdge {
                mesh: self.mesh,
                ind_halfedge: ind_prev,
            })
        } else {
            None
        }
    }

    /// Opposite halfedge: Same vertices in opposite order (on neighbor face)
    pub fn opposite_halfedge(&self) -> Option<IterHalfEdge<'a>> {
        if let Some(&ind_opp) = self.mesh.map_hedg_opp.get(&self.ind_halfedge) {
            Some(IterHalfEdge {
                mesh: self.mesh,
                ind_halfedge: ind_opp,
            })
        } else {
            None
        }
    }

    /// Face containing halfedge
    pub fn face(&self) -> Option<IterFace<'a>> {
        if let Some(&ind_face) = self.mesh.map_hedg_face.get(&self.ind_halfedge) {
            Some(IterFace {
                mesh: self.mesh,
                ind_face,
            })
        } else {
            None
        }
    }
}

impl<'a> IterFace<'a> {
    /// Gets face (array of halfedge indices)
    pub fn face_halfedges(&self) -> FaceHalfedges {
        *self.mesh.faces.get(&self.ind_face).unwrap()
    }

    /// Gets face index
    pub fn ind(&self) -> usize {
        self.ind_face
    }

    /// Surrounding halfedges (array of halfedge iterators)
    pub fn halfedges(&self) -> [IterHalfEdge<'a>; 3] {
        let &face = self.mesh.faces.get(&self.ind_face).unwrap();

        [
            IterHalfEdge {
                mesh: self.mesh,
                ind_halfedge: face[0],
            },
            IterHalfEdge {
                mesh: self.mesh,
                ind_halfedge: face[1],
            },
            IterHalfEdge {
                mesh: self.mesh,
                ind_halfedge: face[2],
            },
        ]
    }

    /// Surrouding vertices (array of vertex iterators)
    pub fn vertices(&self) -> [IterVertex<'a>; 3] {
        let he = self.halfedges();

        [
            he[0].first_vertex(),
            he[1].first_vertex(),
            he[2].first_vertex(),
        ]
    }

    /// Surrouding vertices (array of vertex indices)
    pub fn vertices_inds(&self) -> [usize; 3] {
        let ve = self.vertices();
        [ve[0].ind(), ve[1].ind(), ve[2].ind()]
    }
}
