use anyhow::Result;
use nalgebra::base::*;

/// Mesh vertex
pub type Vertex = Vector3<f64>;
/// Mesh halfedge
pub type HalfEdge = [usize; 2];
/// Mesh face (array of halfedges)
pub type FaceHalfedges = [usize; 3];

#[derive(Clone)]
/// Manifold mesh
pub struct ManifoldMesh3D {
    pub(super) vertices: Vec<Vertex>,
    pub(super) hedg_vert_inds: Vec<usize>, // [ind_f1_v1, ind_f1_v2, ind_f1_v3, ind_f2_v1, ind_f2_v2, ind_f2_v3, ...]
    pub(super) face_groups: Vec<Option<usize>>,

    pub(super) vert_hedg: Vec<Vec<usize>>,
    pub(super) hedg_opp: Vec<Option<usize>>,
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
            vertices: Vec::new(),

            hedg_vert_inds: Vec::new(),
            face_groups: Vec::new(),

            vert_hedg: Vec::new(),
            hedg_opp: Vec::new(),
        }
    }

    /// Adds a vertex to th mesh
    pub fn add_vertex(&mut self, point: &Vector3<f64>) -> usize {
        self.vertices.push(*point);
        self.vert_hedg.push(Vec::new());

        self.vertices.len() - 1
    }

    fn get_vertex_uncheck(&self, ind_vertex: usize) -> IterVertex {
        IterVertex {
            mesh: self,
            ind_vertex,
        }
    }

    /// Vertex getter
    pub fn get_vertex(&self, ind_vertex: usize) -> Result<IterVertex> {
        if ind_vertex >= self.vertices.len() {
            return Err(anyhow::Error::msg("get_vertex(): Index out of bounds"));
        }

        Ok(self.get_vertex_uncheck(ind_vertex))
    }

    /// Gets number of vertices
    pub fn get_nb_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Gets vertex map
    pub fn vertices(&self) -> &Vec<Vertex> {
        &self.vertices
    }

    fn get_halfedge_uncheck(&self, ind_halfedge: usize) -> IterHalfEdge {
        IterHalfEdge {
            mesh: self,
            ind_halfedge,
        }
    }

    /// Halfedge getter
    pub fn get_halfedge(&self, ind_halfedge: usize) -> Result<IterHalfEdge> {
        if ind_halfedge >= self.hedg_vert_inds.len() {
            return Err(anyhow::Error::msg("get_halfedge(): Index out of bounds"));
        }
        Ok(self.get_halfedge_uncheck(ind_halfedge))
    }

    /// Gets number of halfedges
    pub fn get_nb_halfedges(&self) -> usize {
        self.hedg_vert_inds.len()
    }

    /// Adds a face and associated halfedges
    pub fn add_face(
        &mut self,
        ind_vertex1: usize,
        ind_vertex2: usize,
        ind_vertex3: usize,
    ) -> Result<usize> {
        if ind_vertex1 >= self.vertices.len()
            || ind_vertex2 >= self.vertices.len()
            || ind_vertex3 >= self.vertices.len()
        {
            return Err(anyhow::Error::msg("add_face(): Index out of bounds"));
        }

        if self.vert_hedg[ind_vertex1]
            .iter()
            .map(|&ind_he| {
                let ind_v1 = ind_he % 3;
                let ind_f = ind_he - ind_v1;
                let ind_v2 = (ind_v1 + 1) % 3;
                ind_f + ind_v2
            })
            .any(|ind_v2| self.hedg_vert_inds[ind_v2] == ind_vertex2)
        {
            return Err(anyhow::Error::msg("add_face(): halfedge already exists"));
        }

        if self.vert_hedg[ind_vertex2]
            .iter()
            .map(|&ind_he| {
                let ind_v1 = ind_he % 3;
                let ind_f = ind_he - ind_v1;
                let ind_v2 = (ind_v1 + 1) % 3;
                ind_f + ind_v2
            })
            .any(|ind_v2| self.hedg_vert_inds[ind_v2] == ind_vertex3)
        {
            return Err(anyhow::Error::msg("add_face(): halfedge already exists"));
        }

        if self.vert_hedg[ind_vertex3]
            .iter()
            .map(|&ind_he| {
                let ind_v1 = ind_he % 3;
                let ind_f = ind_he - ind_v1;
                let ind_v2 = (ind_v1 + 1) % 3;
                ind_f + ind_v2
            })
            .any(|ind_v2| self.hedg_vert_inds[ind_v2] == ind_vertex1)
        {
            return Err(anyhow::Error::msg("add_face(): halfedge already exists"));
        }

        let ind_he_21_opt = if let Some(he21) = self.is_edge_in(ind_vertex2, ind_vertex1) {
            if he21.opposite_halfedge().is_some() {
                return Err(anyhow::Error::msg(
                    "add_face(): adding face removes manifoldness",
                ));
            }
            Some(he21.ind())
        } else {
            None
        };
        let ind_he_32_opt = if let Some(he32) = self.is_edge_in(ind_vertex3, ind_vertex2) {
            if he32.opposite_halfedge().is_some() {
                return Err(anyhow::Error::msg(
                    "add_face(): adding face removes manifoldness",
                ));
            }
            Some(he32.ind())
        } else {
            None
        };
        let ind_he_13_opt = if let Some(he13) = self.is_edge_in(ind_vertex1, ind_vertex3) {
            if he13.opposite_halfedge().is_some() {
                return Err(anyhow::Error::msg(
                    "add_face(): adding face removes manifoldness",
                ));
            }
            Some(he13.ind())
        } else {
            None
        };

        let ind_he_12 = self.hedg_vert_inds.len();
        self.hedg_vert_inds.push(ind_vertex1);
        let ind_he_23 = self.hedg_vert_inds.len();
        self.hedg_vert_inds.push(ind_vertex2);
        let ind_he_31 = self.hedg_vert_inds.len();
        self.hedg_vert_inds.push(ind_vertex3);

        self.face_groups.push(None);

        self.vert_hedg[ind_vertex1].push(ind_he_12);
        self.vert_hedg[ind_vertex2].push(ind_he_23);
        self.vert_hedg[ind_vertex3].push(ind_he_31);

        self.hedg_opp.push(None);
        self.hedg_opp.push(None);
        self.hedg_opp.push(None);

        if let Some(ind_he_21) = ind_he_21_opt {
            self.hedg_opp[ind_he_12] = Some(ind_he_21);
            self.hedg_opp[ind_he_21] = Some(ind_he_12);
        }
        if let Some(ind_he_32) = ind_he_32_opt {
            self.hedg_opp[ind_he_23] = Some(ind_he_32);
            self.hedg_opp[ind_he_32] = Some(ind_he_23);
        }
        if let Some(ind_he_13) = ind_he_13_opt {
            self.hedg_opp[ind_he_31] = Some(ind_he_13);
            self.hedg_opp[ind_he_13] = Some(ind_he_31);
        }

        Ok((self.hedg_vert_inds.len() - 3) / 3)
    }

    /// removes a face and associated halfedges
    pub fn remove_face(
        &mut self,
        ind_vertex1: usize,
        ind_vertex2: usize,
        ind_vertex3: usize,
    ) -> Result<()> {
        let face = self
            .is_face_in(ind_vertex1, ind_vertex2, ind_vertex3)
            .ok_or(anyhow::Error::msg("Face does not exist in mesh"))?;

        let ind_face = face.ind();

        // Previous face
        let ind_he_ab = ind_face * 3;
        let ind_he_bc = ind_face * 3 + 1;
        let ind_he_ca = ind_face * 3 + 2;

        let ind_va = self.hedg_vert_inds[ind_he_ab];
        let ind_vb = self.hedg_vert_inds[ind_he_bc];
        let ind_vc = self.hedg_vert_inds[ind_he_ca];

        // Delete removed halfedges from vertex-halfedge mappings
        self.vert_hedg[ind_va].retain(|&ind_he| ind_he != ind_he_ab);
        self.vert_hedg[ind_vb].retain(|&ind_he| ind_he != ind_he_bc);
        self.vert_hedg[ind_vc].retain(|&ind_he| ind_he != ind_he_ca);

        // Deletes opposite halfedges from mappings
        let ind_he_ba_opt = self.hedg_opp[ind_he_ab];
        let ind_he_cb_opt = self.hedg_opp[ind_he_bc];
        let ind_he_ac_opt = self.hedg_opp[ind_he_ca];
        if let Some(ind_he_ba) = ind_he_ba_opt {
            self.hedg_opp[ind_he_ba] = None;
        }
        if let Some(ind_he_cb) = ind_he_cb_opt {
            self.hedg_opp[ind_he_cb] = None;
        }
        if let Some(ind_he_ac) = ind_he_ac_opt {
            self.hedg_opp[ind_he_ac] = None;
        }

        if ind_face < self.get_nb_faces() - 1 {
            // Get last face, halfedges and vertices
            let ind_fa_123 = self.hedg_vert_inds.len() - 3;

            let ind_he_12 = ind_fa_123;
            let ind_he_23 = ind_fa_123 + 1;
            let ind_he_31 = ind_fa_123 + 2;

            let ind_v1 = self.hedg_vert_inds[ind_he_12];
            let ind_v2 = self.hedg_vert_inds[ind_he_23];
            let ind_v3 = self.hedg_vert_inds[ind_he_31];

            // Delete removed halfedges from vertex-halfedge mappings
            self.vert_hedg[ind_v1].retain(|&ind_he| ind_he != ind_he_12);
            self.vert_hedg[ind_v2].retain(|&ind_he| ind_he != ind_he_23);
            self.vert_hedg[ind_v3].retain(|&ind_he| ind_he != ind_he_31);

            let ind_he_21_opt = self.hedg_opp[ind_he_12];
            let ind_he_32_opt = self.hedg_opp[ind_he_23];
            let ind_he_13_opt = self.hedg_opp[ind_he_31];

            // New halfedges indices
            let ind_he_12 = ind_he_ab;
            let ind_he_23 = ind_he_bc;
            let ind_he_31 = ind_he_ca;

            // Update the halfedge-vertex mappings
            self.hedg_vert_inds[ind_he_12] = ind_v1;
            self.hedg_vert_inds[ind_he_23] = ind_v2;
            self.hedg_vert_inds[ind_he_31] = ind_v3;

            // Update opposite halfedges mappings
            self.hedg_opp[ind_he_12] = ind_he_21_opt;
            self.hedg_opp[ind_he_23] = ind_he_32_opt;
            self.hedg_opp[ind_he_31] = ind_he_13_opt;
            if let Some(ind_he_21) = ind_he_21_opt {
                self.hedg_opp[ind_he_21] = Some(ind_he_12);
            }
            if let Some(ind_he_32) = ind_he_32_opt {
                self.hedg_opp[ind_he_32] = Some(ind_he_23);
            }
            if let Some(ind_he_13) = ind_he_13_opt {
                self.hedg_opp[ind_he_13] = Some(ind_he_31);
            }
            // Insert added halfedges to vertex-halfedge mappings
            self.vert_hedg[ind_v1].push(ind_he_12);
            self.vert_hedg[ind_v2].push(ind_he_23);
            self.vert_hedg[ind_v3].push(ind_he_31);
        }
        // Pop unused data
        self.hedg_vert_inds.pop();
        self.hedg_vert_inds.pop();
        self.hedg_vert_inds.pop();

        self.hedg_opp.pop();
        self.hedg_opp.pop();
        self.hedg_opp.pop();

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
        if ind_face * 3 >= self.hedg_vert_inds.len() {
            return Err(anyhow::Error::msg("get_face(): Index out of bounds"));
        }
        Ok(self.get_face_uncheck(ind_face))
    }

    /// gets number of faces
    pub fn get_nb_faces(&self) -> usize {
        (self.hedg_vert_inds.len() + 1) / 3
    }

    /// Checks if an edge is in the mesh
    ///
    /// Returns halfedge iterator if found
    pub fn is_edge_in(&self, ind_vertex1: usize, ind_vertex2: usize) -> Option<IterHalfEdge> {
        if ind_vertex1 >= self.vertices.len() || ind_vertex2 >= self.vertices.len() {
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
        if ind_vertex3 >= self.vertices.len() {
            return None;
        } else {
            let opt_he = self.is_edge_in(ind_vertex1, ind_vertex2);
            if let Some(he) = opt_he {
                if he.next_halfedge().last_vertex().ind() == ind_vertex3 {
                    return Some(he.face());
                }
                if let Some(he_opp) = he.opposite_halfedge() {
                    if he_opp.next_halfedge().last_vertex().ind() == ind_vertex3 {
                        return Some(he_opp.face());
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
            let face_comp = hedg.face();
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

        let face = halfedge.face();

        let halfedge_next = halfedge.next_halfedge();

        let halfedge_prev = halfedge.prev_halfedge();

        let face_next = halfedge_next.face();

        let face_prev = halfedge_prev.face();

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
        if let Some(halfedge_opp) = halfedge.opposite_halfedge() {
            //    .ok_or(anyhow::Error::msg("check_halfedge(): no opposite halfedge"))?;

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
        }

        // check vertices
        let neigh_hedges = halfedge.first_vertex().halfedges();

        let is_in = neigh_hedges
            .iter()
            .any(|&iterhedge| iterhedge.ind() == halfedge.ind());

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
        for f in 0..self.get_nb_faces() {
            self.check_face(f)?;
        }

        for e in 0..self.get_nb_halfedges() {
            self.check_halfedge(e)?;
        }

        for v in 0..self.get_nb_vertices() {
            self.check_vertex(v)?;
        }

        Ok(())
    }

    /// Checks if edge is too sharp
    fn check_too_sharp_edge(&self, ind_edge: usize) -> bool {
        let he = self.get_halfedge(ind_edge).unwrap();
        let face = he.face();
        let face_opp = he.opposite_halfedge().unwrap().face();

        let normal = face.normal();
        let normal_opp = face_opp.normal();

        let cos_ang = normal.dot(&normal_opp);

        cos_ang < -0.866
    }

    /// Checks if face has self intersection with edges
    fn check_face_self_inter(&self, ind_face: usize) -> bool {
        let face = self.get_face(ind_face).unwrap();

        let [ind_v1, ind_v2, ind_v3] = face.vertices_inds();
        let vert1 = self.get_vertex(ind_v1).unwrap().vertex();
        let vert2 = self.get_vertex(ind_v2).unwrap().vertex();
        let vert3 = self.get_vertex(ind_v3).unwrap().vertex();

        let vec1 = vert2 - vert1;
        let vec2 = vert3 - vert1;
        let normal = vec1.cross(&vec2).normalize();

        for ind_vert in 0..self.get_nb_vertices() {
            if ind_vert == ind_v1 || ind_vert == ind_v2 || ind_vert == ind_v3 {
                continue;
            }

            let vertex = self.get_vertex(ind_vert).unwrap();
            let coo = vertex.vertex();

            // checks if vertex is above the face
            let above = (coo - vert1).dot(&normal) > 0.0;
            if !above {
                continue;
            }

            // checks for all edges if other extremity is below
            for edg in vertex.halfedges() {
                let vertex_ext2 = edg.last_vertex();
                let ind_vert_ext2 = vertex_ext2.ind();
                if ind_vert_ext2 == ind_v1 || ind_vert_ext2 == ind_v2 || ind_vert_ext2 == ind_v3 {
                    continue;
                }
                let coo_ext2 = vertex_ext2.vertex();

                // checks if vertex is below the face
                let below = (coo_ext2 - vert1).dot(&normal) < 0.0;
                if !below {
                    continue;
                }

                // computes intersection between edge and plane
                let line_vec = (coo_ext2 - coo).normalize();
                let d = (vert1 - coo).dot(&normal) / line_vec.dot(&normal);
                let point = coo + d * line_vec;

                // checks if intersection is in the triangle
                let vert1_mov = vert1 - point;
                let vert2_mov = vert2 - point;
                let vert3_mov = vert3 - point;

                let n1 = vert2_mov.cross(&vert3_mov);
                let n2 = vert3_mov.cross(&vert1_mov);
                let n3 = vert1_mov.cross(&vert2_mov);

                let in_trangle = n1.dot(&n2) > 0.0 && n1.dot(&n3) > 0.0;
                if in_trangle {
                    println!("{}", (coo - vert1).dot(&normal));
                    println!("{}", (point - vert1).dot(&normal));
                    println!("{}", (coo_ext2 - vert1).dot(&normal));

                    print!("face ({} {} {})  ", ind_v1, ind_v2, ind_v3);
                    print!("[{}, {}, {}], ", vert1[0], vert1[1], vert1[2]);
                    print!("[{}, {}, {}], ", vert2[0], vert2[1], vert2[2]);
                    println!("[{}, {}, {}]", vert3[0], vert3[1], vert3[2]);

                    print!("edge ({} {})  ", ind_vert, ind_vert_ext2);
                    print!("[{}, {}, {}], ", coo[0], coo[1], coo[2]);
                    println!("[{}, {}, {}]", coo_ext2[0], coo_ext2[1], coo_ext2[2]);

                    print!("face_mov  ");
                    print!("[{}, {}, {}], ", vert1_mov[0], vert1_mov[1], vert1_mov[2]);
                    print!("[{}, {}, {}], ", vert2_mov[0], vert2_mov[1], vert2_mov[2]);
                    println!("[{}, {}, {}]", vert3_mov[0], vert3_mov[1], vert3_mov[2]);

                    return true;
                }
            }
        }

        false
    }

    /// Checks self intersections of the mesh
    pub fn has_self_intersection(&self) -> bool {
        for f in 0..self.get_nb_faces() {
            if self.check_face_self_inter(f) {
                return true;
            }
        }

        false
    }

    /// Checks sharp edges of the mesh
    pub fn has_sharp_edges(&self) -> bool {
        for e in 0..self.get_nb_halfedges() {
            let [ind1, ind2] = self.get_halfedge(e).unwrap().halfedge();
            if ind1 > ind2 {
                continue;
            }
            if self.check_too_sharp_edge(e) {
                return true;
            }
        }

        false
    }

    /// Assign a group to a face
    pub fn set_face_in_group(&mut self, ind_face: usize, group: usize) {
        self.face_groups[ind_face] = Some(group);
    }
}

impl<'a> IterVertex<'a> {
    /// Gets vertex coordinates
    pub fn vertex(&self) -> Vertex {
        self.mesh.vertices[self.ind_vertex]
    }

    /// Gets vertex index
    pub fn ind(&self) -> usize {
        self.ind_vertex
    }

    /// Gets list of halfedges starting at this vertex
    pub fn halfedges(&self) -> Vec<IterHalfEdge<'a>> {
        self.mesh.vert_hedg[self.ind_vertex]
            .iter()
            .map(|&x| IterHalfEdge {
                mesh: self.mesh,
                ind_halfedge: x,
            })
            .collect()
    }
}

impl<'a> IterHalfEdge<'a> {
    fn ind_first_vertex(&self) -> usize {
        self.mesh.hedg_vert_inds[self.ind_halfedge]
    }

    fn ind_next_halfedge(&self) -> usize {
        let ind_curr = self.ind_halfedge;
        let num_curr = ind_curr % 3;
        let ind_f = ind_curr - num_curr;
        let num_next = (num_curr + 1) % 3;
        ind_f + num_next
    }

    fn ind_last_vertex(&self) -> usize {
        self.mesh.hedg_vert_inds[self.ind_next_halfedge()]
    }

    /// Gets halfedge (array of vertex indices)
    pub fn halfedge(&self) -> HalfEdge {
        [self.ind_first_vertex(), self.ind_last_vertex()]
    }

    /// Gets halfedge index
    pub fn ind(&self) -> usize {
        self.ind_halfedge
    }

    /// First vertex iterator
    pub fn first_vertex(&self) -> IterVertex<'a> {
        IterVertex {
            mesh: self.mesh,
            ind_vertex: self.ind_first_vertex(),
        }
    }

    /// Last vertex iterator
    pub fn last_vertex(&self) -> IterVertex<'a> {
        IterVertex {
            mesh: self.mesh,
            ind_vertex: self.ind_last_vertex(),
        }
    }

    /// Next halfedge on same face
    pub fn next_halfedge(&self) -> IterHalfEdge<'a> {
        IterHalfEdge {
            mesh: self.mesh,
            ind_halfedge: self.ind_next_halfedge(),
        }
    }

    /// Previous halfedge on same face
    pub fn prev_halfedge(&self) -> IterHalfEdge<'a> {
        let ind_curr = self.ind_halfedge;
        let num_curr = ind_curr % 3;
        let ind_f = ind_curr - num_curr;
        let num_prev = (num_curr + 2) % 3;
        let ind_prev = ind_f + num_prev;
        IterHalfEdge {
            mesh: self.mesh,
            ind_halfedge: ind_prev,
        }
    }

    /// Opposite halfedge: Same vertices in opposite order (on neighbor face)
    pub fn opposite_halfedge(&self) -> Option<IterHalfEdge<'a>> {
        if let Some(ind_opp) = self.mesh.hedg_opp[self.ind_halfedge] {
            Some(IterHalfEdge {
                mesh: self.mesh,
                ind_halfedge: ind_opp,
            })
        } else {
            None
        }
    }

    /// Face containing halfedge
    pub fn face(&self) -> IterFace<'a> {
        let ind_face = self.ind_halfedge / 3;
        IterFace {
            mesh: self.mesh,
            ind_face,
        }
    }
}

impl<'a> IterFace<'a> {
    /// Gets face index
    pub fn ind(&self) -> usize {
        self.ind_face
    }

    /// Gets face (array of halfedge indices)
    pub fn halfedges_inds(&self) -> FaceHalfedges {
        [
            self.ind_face * 3,
            self.ind_face * 3 + 1,
            self.ind_face * 3 + 2,
        ]
    }

    /// Surrounding halfedges (array of halfedge iterators)
    pub fn halfedges(&self) -> [IterHalfEdge<'a>; 3] {
        let face_he = self.halfedges_inds();

        [
            IterHalfEdge {
                mesh: self.mesh,
                ind_halfedge: face_he[0],
            },
            IterHalfEdge {
                mesh: self.mesh,
                ind_halfedge: face_he[1],
            },
            IterHalfEdge {
                mesh: self.mesh,
                ind_halfedge: face_he[2],
            },
        ]
    }

    /// Surrouding vertices (array of vertex indices)
    pub fn vertices_inds(&self) -> [usize; 3] {
        let face_he = self.halfedges_inds();
        [
            self.mesh.hedg_vert_inds[face_he[0]],
            self.mesh.hedg_vert_inds[face_he[1]],
            self.mesh.hedg_vert_inds[face_he[2]],
        ]
    }

    /// Surrouding vertices (array of vertex iterators)
    pub fn vertices(&self) -> [IterVertex<'a>; 3] {
        let face_ve = self.vertices_inds();

        [
            IterVertex {
                mesh: self.mesh,
                ind_vertex: face_ve[0],
            },
            IterVertex {
                mesh: self.mesh,
                ind_vertex: face_ve[1],
            },
            IterVertex {
                mesh: self.mesh,
                ind_vertex: face_ve[2],
            },
        ]
    }

    /// Face normal
    pub fn normal(&self) -> Vector3<f64> {
        let [v1, v2, v3] = self.vertices();
        let vert1 = v1.vertex();
        let vert2 = v2.vertex();
        let vert3 = v3.vertex();

        let vec1 = vert2 - vert1;
        let vec2 = vert3 - vert2;

        vec1.cross(&vec2).normalize()
    }
}
