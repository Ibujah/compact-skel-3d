use anyhow::Result;
use nalgebra::base::*;

pub type Vertex = Vector3<f32>;
pub type HalfEdge = [usize; 2];
pub type FaceHalfedges = [usize; 3];
pub type FaceVertices = [usize; 3];

pub struct TopoMesh {
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
    topomesh: &'a TopoMesh,
    ind_vertex: usize,
}

#[derive(Copy, Clone)]
pub struct IterHalfEdge<'a> {
    topomesh: &'a TopoMesh,
    ind_halfedge: usize,
}

#[derive(Copy, Clone)]
pub struct IterFace<'a> {
    topomesh: &'a TopoMesh,
    ind_face: usize,
}

impl TopoMesh {
    pub fn init() -> TopoMesh {
        TopoMesh {
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

    pub fn get_vertex(&self, ind_vertex: usize) -> Result<IterVertex> {
        if ind_vertex >= self.vertices.len() {
            return Err(anyhow::Error::msg("get_vertex(): Index out of bounds"));
        }

        Ok(IterVertex {
            topomesh: self,
            ind_vertex,
        })
    }

    pub fn get_nb_vertices(&self) -> usize {
        self.vertices.len()
    }

    pub fn add_halfedge(&mut self, ind_vertex1: usize, ind_vertex2: usize) -> Result<usize> {
        if ind_vertex1 >= self.map_vert_hedg.len() || ind_vertex2 >= self.map_vert_hedg.len() {
            return Err(anyhow::Error::msg("add_halfedge(): Vertex index out of bounds"));
        }
        let ind_halfedge = 
            self.map_vert_hedg[ind_vertex1]
            .iter()
            .fold(None,
                  |res, &ind_he|
                  {
                      if self.halfedges[ind_he][1] == ind_vertex2 {
                        Some(ind_he)
                      }
                      else {
                        res
                      }
                  }
            );

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

    pub fn get_halfedge(&self, ind_halfedge: usize) -> Result<IterHalfEdge> {
        if ind_halfedge >= self.halfedges.len() {
            return Err(anyhow::Error::msg("get_halfedge(): Index out of bounds"));
        }
        Ok(IterHalfEdge {
            topomesh: self,
            ind_halfedge,
        })
    }

    pub fn get_nb_halfedges(&self) -> usize {
        self.halfedges.len()
    }
    
    fn fill_face(&mut self, ind_face: usize, ind_halfedge1: usize, ind_halfedge2: usize, ind_halfedge3: usize, ind_halfedge1_opp: usize, ind_halfedge2_opp: usize, ind_halfedge3_opp: usize) -> () {
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

    pub fn add_face(&mut self, ind_vertex1: usize, ind_vertex2: usize, ind_vertex3: usize) -> Result<usize> {
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
            return Err(anyhow::Error::msg("add_face(): Error in faces organisation"));
        }

        if let Some(fac) = ind_f1 {
            return Ok(fac)
        }

        self.faces.push([ind_halfedge1, ind_halfedge2, ind_halfedge3]);
        let ind_face = self.faces.len() - 1;

        self.fill_face(ind_face, ind_halfedge1, ind_halfedge2, ind_halfedge3, ind_halfedge1_opp, ind_halfedge2_opp, ind_halfedge3_opp);
        
        Ok(ind_face)
    }

    pub fn get_face(&self, ind_face: usize) -> Result<IterFace> {
        if ind_face >= self.faces.len() {
            return Err(anyhow::Error::msg("get_face(): Index out of bounds"));
        }
        Ok(IterFace {
            topomesh: self,
            ind_face,
        })
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

    pub fn can_flip_halfedge(&self, ind_halfedge: usize) -> Result<bool> {
        let halfedge = self.get_halfedge(ind_halfedge)?;
        
        let opp_vert1 = 
            halfedge
            .next_halfedge()?
            .last_vertex();
        let opp_vert2 = 
            halfedge
            .opposite_halfedge()?
            .next_halfedge()?
            .last_vertex();
        
        let edg_found = 
            opp_vert1
            .halfedges()
            .iter()
            .fold(false,
                  |res, &he| {
                      res || he.last_vertex().ind() == opp_vert2.ind()
                  });

        Ok(!edg_found)
    }

    pub fn flip_halfedge(&mut self, ind_halfedge: usize) -> Result<bool> {
        if ind_halfedge >= self.halfedges.len() {
            return Err(anyhow::Error::msg("flip_halfedge(): Index out of bounds"));
        }
        let flippable = self.can_flip_halfedge(ind_halfedge)?;

        if !flippable {
            return Ok(false)
        }
        
        // direct order
        //     1             1   
        //   / | \         /   \ 
        //  4  |  3  -->  4 --- 3
        //   \ | /         \   / 
        //     2             2   

        let ind_he_12 = ind_halfedge;
        let ind_fac_123 = self.map_hedg_face[ind_he_12].ok_or(anyhow::Error::msg("flip_halfedge(): Empty face"))?;
        let ind_he_23 = self.map_hedg_next[ind_he_12].ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
        let ind_he_31 = self.map_hedg_prev[ind_he_12].ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;

        let ind_he_21 = self.map_hedg_opp[ind_he_12].ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
        let ind_fac_214 = self.map_hedg_face[ind_he_21].ok_or(anyhow::Error::msg("flip_halfedge(): Empty face"))?;
        let ind_he_14 = self.map_hedg_next[ind_he_21].ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
        let ind_he_42 = self.map_hedg_prev[ind_he_21].ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
        
        let ind_v1 = self.halfedges[ind_he_12][0];
        let ind_v2 = self.halfedges[ind_he_12][1];
        let ind_v3 = self.halfedges[ind_he_23][1];
        let ind_v4 = self.halfedges[ind_he_14][1];
        
        // rename edges and faces
        let ind_fac_143 = ind_fac_123;
        let ind_fac_234 = ind_fac_214;

        let ind_he_43 = ind_he_12;
        let ind_he_34 = ind_he_21;

        // update faces
        self.faces[ind_fac_143] = [ind_he_14, ind_he_43, ind_he_31];
        self.faces[ind_fac_234] = [ind_he_23, ind_he_34, ind_he_42];

        // update edges
        self.halfedges[ind_he_43] = [ind_v4, ind_v3];
        self.halfedges[ind_he_34] = [ind_v3, ind_v4];
        
        // update connectivity
        self.map_hedg_face[ind_he_14] = Some(ind_fac_143);
        self.map_hedg_face[ind_he_23] = Some(ind_fac_234);

        self.map_hedg_next[ind_he_43] = Some(ind_he_31);
        self.map_hedg_next[ind_he_31] = Some(ind_he_14);
        self.map_hedg_next[ind_he_14] = Some(ind_he_43);

        self.map_hedg_prev[ind_he_43] = Some(ind_he_14);
        self.map_hedg_prev[ind_he_14] = Some(ind_he_31);
        self.map_hedg_prev[ind_he_31] = Some(ind_he_43);

        self.map_hedg_next[ind_he_34] = Some(ind_he_42);
        self.map_hedg_next[ind_he_42] = Some(ind_he_23);
        self.map_hedg_next[ind_he_23] = Some(ind_he_34);

        self.map_hedg_prev[ind_he_34] = Some(ind_he_23);
        self.map_hedg_prev[ind_he_23] = Some(ind_he_42);
        self.map_hedg_prev[ind_he_42] = Some(ind_he_34);
        
        self.map_vert_hedg[ind_v1].retain(|&ind| ind != ind_he_12);
        self.map_vert_hedg[ind_v2].retain(|&ind| ind != ind_he_21);
        self.map_vert_hedg[ind_v3].push(ind_he_34);
        self.map_vert_hedg[ind_v4].push(ind_he_43);
        
        Ok(true)
    }

    pub fn split_halfedge(&mut self, vert: &Vertex, ind_halfedge: usize) -> Result<()> {
        if ind_halfedge >= self.halfedges.len() {
            return Err(anyhow::Error::msg("split_halfedge(): Index out of bounds"));
        }
        
        // direct order
        //     1             1   
        //   / | \         / | \ 
        //  4  |  3  -->  4 -5- 3
        //   \ | /         \ | / 
        //     2             2   

        let ind_he_12 = ind_halfedge;
        let ind_fac_123 = self.map_hedg_face[ind_he_12].ok_or(anyhow::Error::msg("flip_halfedge(): Empty face"))?;
        let ind_he_23 = self.map_hedg_next[ind_he_12].ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
        let ind_he_31 = self.map_hedg_prev[ind_he_12].ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;

        let ind_he_21 = self.map_hedg_opp[ind_he_12].ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
        let ind_fac_214 = self.map_hedg_face[ind_he_21].ok_or(anyhow::Error::msg("flip_halfedge(): Empty face"))?;
        let ind_he_14 = self.map_hedg_next[ind_he_21].ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
        let ind_he_42 = self.map_hedg_prev[ind_he_21].ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;

        let ind_he_13 = self.map_hedg_opp[ind_he_31].ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
        let ind_he_32 = self.map_hedg_opp[ind_he_23].ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
        let ind_he_24 = self.map_hedg_opp[ind_he_42].ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
        let ind_he_41 = self.map_hedg_opp[ind_he_14].ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
        
        let ind_v1 = self.halfedges[ind_he_12][0];
        let ind_v2 = self.halfedges[ind_he_12][1];
        let ind_v3 = self.halfedges[ind_he_23][1];
        let ind_v4 = self.halfedges[ind_he_14][1];
        let ind_v5 = self.add_vertex(vert);

        self.map_vert_hedg[ind_v1].retain(|&ind| ind != ind_he_12);
        self.map_vert_hedg[ind_v2].retain(|&ind| ind != ind_he_21);

        // rename edges
        let ind_he_52 = ind_he_12;
        let ind_he_51 = ind_he_21;
        
        // reinit edges
        self.halfedges[ind_he_51] = [ind_v5, ind_v1];
        self.halfedges[ind_he_52] = [ind_v5, ind_v2];

        self.map_hedg_face[ind_he_51] = None;
        self.map_hedg_face[ind_he_52] = None;

        self.map_hedg_prev[ind_he_51] = None;
        self.map_hedg_prev[ind_he_52] = None;

        self.map_hedg_next[ind_he_51] = None;
        self.map_hedg_next[ind_he_52] = None;

        self.map_hedg_opp[ind_he_51] = None;
        self.map_hedg_opp[ind_he_52] = None;

        self.map_vert_hedg[ind_v5].push(ind_he_51);
        self.map_vert_hedg[ind_v5].push(ind_he_52);
        
        // add new halfedges
        let ind_he_53 = self.add_halfedge(ind_v5, ind_v3)?;
        let ind_he_54 = self.add_halfedge(ind_v5, ind_v4)?;
        let ind_he_15 = self.add_halfedge(ind_v1, ind_v5)?;
        let ind_he_25 = self.add_halfedge(ind_v2, ind_v5)?;
        let ind_he_35 = self.add_halfedge(ind_v3, ind_v5)?;
        let ind_he_45 = self.add_halfedge(ind_v4, ind_v5)?;

        // rename faces
        let ind_fac_235 = ind_fac_123;
        let ind_fac_145 = ind_fac_214;

        // reinit faces
        self.faces[ind_fac_235] = [ind_he_23, ind_he_35, ind_he_52];
        self.faces[ind_fac_145] = [ind_he_14, ind_he_45, ind_he_51];

        // add new faces
        self.faces.push([ind_he_15, ind_he_53, ind_he_31]);
        let ind_fac_153 = self.faces.len() - 1;

        self.faces.push([ind_he_25, ind_he_54, ind_he_42]);
        let ind_fac_254 = self.faces.len() - 1;
        
        self.fill_face(ind_fac_235, 
                       ind_he_23, ind_he_35, ind_he_52, 
                       ind_he_32, ind_he_53, ind_he_25);
        self.fill_face(ind_fac_145, 
                       ind_he_14, ind_he_45, ind_he_51, 
                       ind_he_41, ind_he_54, ind_he_15);
        self.fill_face(ind_fac_153, 
                       ind_he_15, ind_he_53, ind_he_31, 
                       ind_he_51, ind_he_35, ind_he_13);
        self.fill_face(ind_fac_254, 
                       ind_he_25, ind_he_54, ind_he_42, 
                       ind_he_52, ind_he_45, ind_he_24);

        Ok(())
    }

    pub fn split_face(&mut self, vert: &Vertex, ind_face: usize) -> Result<()> {
        if ind_face >= self.faces.len() {
            return Err(anyhow::Error::msg("split_face(): Index out of bounds"));
        }
        
        // direct order
        //      1                1       
        //    /   \            / | \     
        //   /     \   -->    /  4  \    
        //  /       \        / /   \ \   
        // 2 ------- 3      2 ------- 3 
        
        let ind_fac_123 = ind_face;
        let [ind_he_12, ind_he_23, ind_he_31] = self.get_face(ind_face)?.face_halfedges();
        let ind_he_21 = self.map_hedg_opp[ind_he_12].ok_or(anyhow::Error::msg("split_face(): Empty halfedge"))?;
        let ind_he_32 = self.map_hedg_opp[ind_he_23].ok_or(anyhow::Error::msg("split_face(): Empty halfedge"))?;
        let ind_he_13 = self.map_hedg_opp[ind_he_31].ok_or(anyhow::Error::msg("split_face(): Empty halfedge"))?;
        
        let ind_v1 = self.halfedges[ind_he_12][0];
        let ind_v2 = self.halfedges[ind_he_23][0];
        let ind_v3 = self.halfedges[ind_he_31][0];
        let ind_v4 = self.add_vertex(vert);

        // rename face
        let ind_fac_124 = ind_fac_123;

        // add halfedges
        let ind_he_14 = self.add_halfedge(ind_v1, ind_v4)?;
        let ind_he_41 = self.add_halfedge(ind_v4, ind_v1)?;
        let ind_he_24 = self.add_halfedge(ind_v2, ind_v4)?;
        let ind_he_42 = self.add_halfedge(ind_v4, ind_v2)?;
        let ind_he_34 = self.add_halfedge(ind_v3, ind_v4)?;
        let ind_he_43 = self.add_halfedge(ind_v4, ind_v3)?;
        
        // reinit faces
        self.faces[ind_fac_124] = [ind_he_12, ind_he_24, ind_he_41];

        // add faces
        self.faces.push([ind_he_23, ind_he_34, ind_he_42]);
        let ind_fac_234 = self.faces.len() - 1;
        self.faces.push([ind_he_14, ind_he_43, ind_he_31]);
        let ind_fac_143 = self.faces.len() - 1;
        
        // fill faces
        self.fill_face(ind_fac_124, 
                       ind_he_12, ind_he_24, ind_he_41, 
                       ind_he_21, ind_he_42, ind_he_14);

        self.fill_face(ind_fac_234, 
                       ind_he_23, ind_he_34, ind_he_42, 
                       ind_he_32, ind_he_43, ind_he_24);

        self.fill_face(ind_fac_143, 
                       ind_he_14, ind_he_43, ind_he_31, 
                       ind_he_41, ind_he_34, ind_he_13);

        Ok(())
    }
    
    fn check_face(&self, ind_face: usize) -> Result<()> {
        let face = self.get_face(ind_face)?;
        // check edges existence
        for hedg in face.halfedges() {
            let face_comp = hedg.face()?;
            if face_comp.ind() != face.ind() {
                return Err(anyhow::Error::msg("check_face(): HalfEdge linked to wrong face"));
            }
        }
        Ok(())
    }
    
    fn check_halfedge(&self, ind_hedge: usize) -> Result<()> {
        let halfedge = self.get_halfedge(ind_hedge)?;
        
        let face = halfedge.face()?;

        let halfedge_next = halfedge.next_halfedge()?;
        let halfedge_prev = halfedge.prev_halfedge()?;

        let face_next = halfedge_next.face()?;
        let face_prev = halfedge_prev.face()?;

        if halfedge.last_vertex().ind() != halfedge_next.first_vertex().ind() {
            return Err(anyhow::Error::msg("check_halfedge(): Next halfedge not starting with last vertex"));
        }
        if halfedge.first_vertex().ind() != halfedge_prev.last_vertex().ind() {
            return Err(anyhow::Error::msg("check_halfedge(): Previous halfedge not ending with first vertex"));
        }

        if face.ind() != face_next.ind() {
            return Err(anyhow::Error::msg("check_halfedge(): Next halfedge not on the same face"));
        }
        if face.ind() != face_prev.ind() {
            return Err(anyhow::Error::msg("check_halfedge(): Previous halfedge not on the same face"));
        }

        // check opposite
        let halfedge_opp = halfedge.opposite_halfedge()?;
        if halfedge.first_vertex().ind() != halfedge_opp.last_vertex().ind() {
            return Err(anyhow::Error::msg("check_halfedge(): Opposite halfedge not starting with last vertex"));
        }
        if halfedge.last_vertex().ind() != halfedge_opp.first_vertex().ind() {
            return Err(anyhow::Error::msg("check_halfedge(): Opposite halfedge not ending with first vertex"));
        }

        // check vertices
        let neigh_hedges = halfedge.first_vertex().halfedges();
        
        let is_in = neigh_hedges
            .iter()
            .fold(false,
                  |res, &iterhedge| {
                      res || iterhedge.ind() == halfedge.ind()
                  });

        if !is_in {
            return Err(anyhow::Error::msg("check_halfedge(): Halfedge not in vertex"));
        }

        Ok(())
    }

    fn check_vertex(&self, ind_vertex: usize) -> Result<()> {
        let vertex = self.get_vertex(ind_vertex)?;
        
        for he in vertex.halfedges() {
            if he.first_vertex().ind() != ind_vertex {
                return Err(anyhow::Error::msg("check_vertex(): Vertex contains non coherent halfedge"));
            }
        }

        Ok(())
    }

    pub fn check_topo_mesh(&self) -> Result<()> {
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
        self.topomesh.vertices[self.ind_vertex]
    }

    pub fn ind(&self) -> usize {
        self.ind_vertex
    }

    pub fn halfedges(&self) -> Vec<IterHalfEdge<'a>> {
        let vec_he =
            self.topomesh.map_vert_hedg[self.ind_vertex] 
            .iter()
            .fold(Vec::new(), 
                  |mut v, &x| {v.push(
                          IterHalfEdge{
                              topomesh: self.topomesh,
                              ind_halfedge: x,
                          }); v});
        vec_he
    }
}

impl<'a> IterHalfEdge<'a> {
    pub fn halfedge(&self) -> HalfEdge {
        self.topomesh.halfedges[self.ind_halfedge]
    }

    pub fn ind(&self) -> usize {
        self.ind_halfedge
    }

    pub fn first_vertex(&self) -> IterVertex<'a> {
        IterVertex {
            topomesh: self.topomesh,
            ind_vertex: self.halfedge()[0],
        }
    }

    pub fn last_vertex(&self) -> IterVertex<'a> {
        IterVertex {
            topomesh: self.topomesh,
            ind_vertex: self.halfedge()[1],
        }
    }

    pub fn next_halfedge(&self) -> Result<IterHalfEdge<'a>> {
        let ind_next = self.topomesh.map_hedg_next[self.ind_halfedge]
            .ok_or(anyhow::Error::msg("next_halfedge(): No next halfedge"))?;
        Ok(IterHalfEdge {
            topomesh: self.topomesh,
            ind_halfedge: ind_next,
        })
    }

    pub fn prev_halfedge(&self) -> Result<IterHalfEdge<'a>> {
        let ind_prev = self.topomesh.map_hedg_prev[self.ind_halfedge]
            .ok_or(anyhow::Error::msg("prev_halfedge(): No previous halfedge"))?;
        Ok(IterHalfEdge {
            topomesh: self.topomesh,
            ind_halfedge: ind_prev,
        })
    }

    pub fn opposite_halfedge(&self) -> Result<IterHalfEdge<'a>> {
        let ind_opp = self.topomesh.map_hedg_opp[self.ind_halfedge]
            .ok_or(anyhow::Error::msg("opposite_halfedge(): No opposite halfedge"))?;
        Ok(IterHalfEdge {
            topomesh: self.topomesh,
            ind_halfedge: ind_opp,
        })
    }

    pub fn face(&self) -> Result<IterFace<'a>> {
        let ind_face = self.topomesh.map_hedg_face[self.ind_halfedge]
            .ok_or(anyhow::Error::msg("face(): No associated face"))?;
        Ok(IterFace{
            topomesh: self.topomesh,
            ind_face,
        })
    }
}

impl<'a> IterFace<'a> {
    
    pub fn face_halfedges(&self) -> FaceHalfedges {
        self.topomesh.faces[self.ind_face]
    }

    pub fn halfedges(&self) -> [IterHalfEdge<'a>; 3] {
        let face = self.topomesh.faces[self.ind_face];

        [
            IterHalfEdge{
                topomesh: self.topomesh,
                ind_halfedge: face[0],
            },
            IterHalfEdge{
                topomesh: self.topomesh,
                ind_halfedge: face[1],
            },
            IterHalfEdge{
                topomesh: self.topomesh,
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

    pub fn vertices_inds(&self) -> [usize; 3]{
        let ve = self.vertices();
        [
            ve[0].ind(),
            ve[1].ind(),
            ve[2].ind(),
        ]
    }

    pub fn ind(&self) -> usize {
        self.ind_face
    }
}
