use anyhow::Result;
use nalgebra::base::*;

pub type Vertex = Vector3<f32>;
pub type HalfEdge = [usize; 2];
pub type FaceHalfedges = [usize; 3];
pub type FaceVertices = [usize; 3];

pub struct Mesh3D {
    pub(super) vertices: Vec<Vertex>,
    pub(super) halfedges: Vec<HalfEdge>,
    pub(super) faces: Vec<FaceHalfedges>,

    pub(super) map_vert_hedg: Vec<Vec<usize>>,
    pub(super) map_hedg_face: Vec<Option<usize>>,
    pub(super) map_hedg_opp: Vec<Option<usize>>,
    pub(super) map_hedg_next: Vec<Option<usize>>,
    pub(super) map_hedg_prev: Vec<Option<usize>>,
}

#[derive(Copy, Clone)]
pub struct IterVertex<'a> {
    mesh: &'a Mesh3D,
    ind_vertex: usize,
}

#[derive(Copy, Clone)]
pub struct IterHalfEdge<'a> {
    mesh: &'a Mesh3D,
    ind_halfedge: usize,
}

#[derive(Copy, Clone)]
pub struct IterFace<'a> {
    mesh: &'a Mesh3D,
    ind_face: usize,
}

impl Mesh3D {
    pub fn new() -> Mesh3D {
        Mesh3D {
            vertices: Vec::new(),
            halfedges: Vec::new(),
            faces: Vec::new(),

            map_vert_hedg: Vec::new(),
            map_hedg_face: Vec::new(),
            map_hedg_opp: Vec::new(),
            map_hedg_next: Vec::new(),
            map_hedg_prev: Vec::new(),
        }
    }

    pub fn add_vertex(&mut self, point: &Vector3<f32>) -> usize {
        self.vertices.push(*point);
        self.map_vert_hedg.push(Vec::new());
        self.vertices.len() - 1
    }

    fn get_vertex_uncheck(&self, ind_vertex: usize) -> IterVertex {
        IterVertex {
            mesh: self,
            ind_vertex,
        }
    }

    pub fn get_vertex(&self, ind_vertex: usize) -> Result<IterVertex> {
        if ind_vertex >= self.vertices.len() {
            return Err(anyhow::Error::msg("get_vertex(): Index out of bounds"));
        }

        Ok(self.get_vertex_uncheck(ind_vertex))
    }

    pub fn get_nb_vertices(&self) -> usize {
        self.vertices.len()
    }

    pub fn add_halfedge(&mut self, ind_vertex1: usize, ind_vertex2: usize) -> Result<usize> {
        if ind_vertex1 >= self.map_vert_hedg.len() || ind_vertex2 >= self.map_vert_hedg.len() {
            return Err(anyhow::Error::msg(
                "add_halfedge(): Vertex index out of bounds",
            ));
        }
        let ind_halfedge = self.map_vert_hedg[ind_vertex1]
            .iter()
            .fold(None, |res, &ind_he| {
                if self.halfedges[ind_he][1] == ind_vertex2 {
                    Some(ind_he)
                } else {
                    res
                }
            });

        match ind_halfedge {
            Some(ind) => Ok(ind),
            None => {
                self.halfedges.push([ind_vertex1, ind_vertex2]);
                self.map_hedg_face.push(None);
                self.map_hedg_prev.push(None);
                self.map_hedg_next.push(None);
                self.map_hedg_opp.push(None);
                self.map_vert_hedg[ind_vertex1].push(self.halfedges.len() - 1);
                Ok(self.halfedges.len() - 1)
            }
        }
    }

    fn get_halfedge_uncheck(&self, ind_halfedge: usize) -> IterHalfEdge {
        IterHalfEdge {
            mesh: self,
            ind_halfedge,
        }
    }

    pub fn get_halfedge(&self, ind_halfedge: usize) -> Result<IterHalfEdge> {
        if ind_halfedge >= self.halfedges.len() {
            return Err(anyhow::Error::msg("get_halfedge(): Index out of bounds"));
        }
        Ok(self.get_halfedge_uncheck(ind_halfedge))
    }

    pub fn get_nb_halfedges(&self) -> usize {
        self.halfedges.len()
    }

    pub(super) fn fill_face(
        &mut self,
        ind_face: usize,
        ind_halfedge1: usize,
        ind_halfedge2: usize,
        ind_halfedge3: usize,
        ind_halfedge1_opp: usize,
        ind_halfedge2_opp: usize,
        ind_halfedge3_opp: usize,
    ) -> () {
        self.map_hedg_face[ind_halfedge1] = Some(ind_face);
        self.map_hedg_face[ind_halfedge2] = Some(ind_face);
        self.map_hedg_face[ind_halfedge3] = Some(ind_face);

        self.map_hedg_next[ind_halfedge1] = Some(ind_halfedge2);
        self.map_hedg_next[ind_halfedge2] = Some(ind_halfedge3);
        self.map_hedg_next[ind_halfedge3] = Some(ind_halfedge1);

        self.map_hedg_prev[ind_halfedge1] = Some(ind_halfedge3);
        self.map_hedg_prev[ind_halfedge2] = Some(ind_halfedge1);
        self.map_hedg_prev[ind_halfedge3] = Some(ind_halfedge2);

        self.map_hedg_opp[ind_halfedge1] = Some(ind_halfedge1_opp);
        self.map_hedg_opp[ind_halfedge2] = Some(ind_halfedge2_opp);
        self.map_hedg_opp[ind_halfedge3] = Some(ind_halfedge3_opp);

        self.map_hedg_opp[ind_halfedge1_opp] = Some(ind_halfedge1);
        self.map_hedg_opp[ind_halfedge2_opp] = Some(ind_halfedge2);
        self.map_hedg_opp[ind_halfedge3_opp] = Some(ind_halfedge3);
    }

    pub fn add_face(
        &mut self,
        ind_vertex1: usize,
        ind_vertex2: usize,
        ind_vertex3: usize,
    ) -> Result<usize> {
        let ind_halfedge1 = self.add_halfedge(ind_vertex1, ind_vertex2)?;
        let ind_halfedge2 = self.add_halfedge(ind_vertex2, ind_vertex3)?;
        let ind_halfedge3 = self.add_halfedge(ind_vertex3, ind_vertex1)?;
        let ind_halfedge1_opp = self.add_halfedge(ind_vertex2, ind_vertex1)?;
        let ind_halfedge2_opp = self.add_halfedge(ind_vertex3, ind_vertex2)?;
        let ind_halfedge3_opp = self.add_halfedge(ind_vertex1, ind_vertex3)?;

        let ind_f1 = self.map_hedg_face[ind_halfedge1];
        let ind_f2 = self.map_hedg_face[ind_halfedge2];
        let ind_f3 = self.map_hedg_face[ind_halfedge3];

        if ind_f1 != ind_f2 || ind_f1 != ind_f3 {
            return Err(anyhow::Error::msg(
                "add_face(): Error in faces organisation",
            ));
        }

        if let Some(fac) = ind_f1 {
            return Ok(fac);
        }

        self.faces
            .push([ind_halfedge1, ind_halfedge2, ind_halfedge3]);
        let ind_face = self.faces.len() - 1;

        self.fill_face(
            ind_face,
            ind_halfedge1,
            ind_halfedge2,
            ind_halfedge3,
            ind_halfedge1_opp,
            ind_halfedge2_opp,
            ind_halfedge3_opp,
        );

        Ok(ind_face)
    }

    fn get_face_uncheck(&self, ind_face: usize) -> IterFace {
        IterFace {
            mesh: self,
            ind_face,
        }
    }

    pub fn get_face(&self, ind_face: usize) -> Result<IterFace> {
        if ind_face >= self.faces.len() {
            return Err(anyhow::Error::msg("get_face(): Index out of bounds"));
        }
        Ok(self.get_face_uncheck(ind_face))
    }

    pub fn get_nb_faces(&self) -> usize {
        self.faces.len()
    }

    pub fn get_face_vertices(&self, ind_face: usize) -> Result<FaceVertices> {
        if ind_face >= self.faces.len() {
            return Err(anyhow::Error::msg("get_face(): Index out of bounds"));
        }

        let face_he = self.faces[ind_face];
        let he1 = self.halfedges[face_he[0]];
        let he2 = self.halfedges[face_he[1]];
        let he3 = self.halfedges[face_he[2]];

        let face_v = [he1[0], he2[0], he3[0]];

        Ok(face_v)
    }

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

    pub fn check_mesh(&self) -> Result<()> {
        for f in 0..self.faces.len() {
            self.check_face(f)?;
        }

        for e in 0..self.halfedges.len() {
            self.check_halfedge(e)?;
        }

        for v in 0..self.vertices.len() {
            self.check_vertex(v)?;
        }

        Ok(())
    }
}

impl<'a> IterVertex<'a> {
    pub fn vertex(&self) -> Vertex {
        self.mesh.vertices[self.ind_vertex]
    }

    pub fn ind(&self) -> usize {
        self.ind_vertex
    }

    pub fn halfedges(&self) -> Vec<IterHalfEdge<'a>> {
        let vec_he =
            self.mesh.map_vert_hedg[self.ind_vertex]
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
    pub fn halfedge(&self) -> HalfEdge {
        self.mesh.halfedges[self.ind_halfedge]
    }

    pub fn ind(&self) -> usize {
        self.ind_halfedge
    }

    pub fn first_vertex(&self) -> IterVertex<'a> {
        IterVertex {
            mesh: self.mesh,
            ind_vertex: self.halfedge()[0],
        }
    }

    pub fn last_vertex(&self) -> IterVertex<'a> {
        IterVertex {
            mesh: self.mesh,
            ind_vertex: self.halfedge()[1],
        }
    }

    pub fn next_halfedge(&self) -> Option<IterHalfEdge<'a>> {
        if let Some(ind_next) = self.mesh.map_hedg_next[self.ind_halfedge] {
            Some(IterHalfEdge {
                mesh: self.mesh,
                ind_halfedge: ind_next,
            })
        } else {
            None
        }
    }

    pub fn prev_halfedge(&self) -> Option<IterHalfEdge<'a>> {
        if let Some(ind_prev) = self.mesh.map_hedg_prev[self.ind_halfedge] {
            Some(IterHalfEdge {
                mesh: self.mesh,
                ind_halfedge: ind_prev,
            })
        } else {
            None
        }
    }

    pub fn opposite_halfedge(&self) -> Option<IterHalfEdge<'a>> {
        if let Some(ind_opp) = self.mesh.map_hedg_opp[self.ind_halfedge] {
            Some(IterHalfEdge {
                mesh: self.mesh,
                ind_halfedge: ind_opp,
            })
        } else {
            None
        }
    }

    pub fn face(&self) -> Option<IterFace<'a>> {
        if let Some(ind_face) = self.mesh.map_hedg_face[self.ind_halfedge] {
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
    pub fn face_halfedges(&self) -> FaceHalfedges {
        self.mesh.faces[self.ind_face]
    }

    pub fn halfedges(&self) -> [IterHalfEdge<'a>; 3] {
        let face = self.mesh.faces[self.ind_face];

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

    pub fn vertices(&self) -> [IterVertex<'a>; 3] {
        let he = self.halfedges();

        [
            he[0].first_vertex(),
            he[1].first_vertex(),
            he[2].first_vertex(),
        ]
    }

    pub fn vertices_inds(&self) -> [usize; 3] {
        let ve = self.vertices();
        [ve[0].ind(), ve[1].ind(), ve[2].ind()]
    }

    pub fn ind(&self) -> usize {
        self.ind_face
    }
}
