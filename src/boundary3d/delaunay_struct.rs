use anyhow::Result;
use tritet::{StrError, Tetgen};
use std::collections::{HashSet, HashMap};

use crate::boundary3d::topo_mesh::{self, TopoMesh, HalfEdge, FaceVertices};

pub type Edge = [usize; 2];
pub type Triangle = [usize; 3];
pub type Tetrahedra = [usize; 4];

pub struct DelaunayStruct<'a>{
    topomesh: &'a mut TopoMesh,

    edges: HashSet<Edge>,
    faces: HashMap<Triangle, Vec<Tetrahedra> >,
    tetras: HashSet<Tetrahedra>,
    
    initial_vertices_number: usize,
}

fn to_anyhow(err: StrError) -> anyhow::Error
{
    anyhow::Error::msg(err.to_string())
}

impl<'a> DelaunayStruct<'a> {
    
    fn insert_tetra(&mut self, tetra: &mut Tetrahedra) -> () {
        tetra.sort();

        self.edges.insert([tetra[0], tetra[1]]);
        self.edges.insert([tetra[0], tetra[2]]);
        self.edges.insert([tetra[0], tetra[3]]);
        self.edges.insert([tetra[1], tetra[2]]);
        self.edges.insert([tetra[1], tetra[3]]);
        self.edges.insert([tetra[2], tetra[3]]);
        
        self.faces.entry([tetra[0], tetra[1], tetra[2]]).or_insert(Vec::new()).push(*tetra);
        self.faces.entry([tetra[0], tetra[1], tetra[3]]).or_insert(Vec::new()).push(*tetra);
        self.faces.entry([tetra[0], tetra[2], tetra[3]]).or_insert(Vec::new()).push(*tetra);
        self.faces.entry([tetra[1], tetra[2], tetra[3]]).or_insert(Vec::new()).push(*tetra);

        self.tetras.insert([tetra[0], tetra[1], tetra[2], tetra[3]]);
    }
    
    fn generate_struct(&mut self) -> Result<()> {
        let mut tetgen = 
            Tetgen::new(self.topomesh.get_nb_vertices(), 
                        Some(vec![3; self.topomesh.get_nb_faces()]), 
                        None, 
                        None)
            .map_err(to_anyhow)?;

        for v in 0..self.topomesh.get_nb_vertices() {
            let vert = self.topomesh.get_vertex(v)?;

            tetgen
                .set_point(v, vert[0] as f64, vert[1] as f64, vert[2] as f64)
                .map_err(to_anyhow)?;
        }

        tetgen.generate_delaunay(false)
            .map_err(to_anyhow)?;

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

    pub fn init(topomesh: &'a mut TopoMesh) -> Result<DelaunayStruct<'a>> {
        let initial_vertices_number = topomesh.get_nb_vertices();
        let mut deltet = DelaunayStruct { 
            topomesh,
            edges: HashSet::new(),
            faces: HashMap::new(),
            tetras: HashSet::new(),
            initial_vertices_number,
        };

        deltet.generate_struct()?;

        Ok(deltet)
    }
    
    pub fn get_topomesh(&self) -> &TopoMesh {
        self.topomesh
    }

    pub fn is_original_vertex(&self, ind_vertex: usize) -> bool {
        ind_vertex < self.initial_vertices_number
    }

    pub fn is_edge_in(&self, edge: &Edge) -> bool {
        let mut edge_sort = [edge[0], edge[1]];
        edge_sort.sort();
        self.edges.contains(&edge_sort)
    }

    pub fn is_face_in(&self, face: &Triangle) -> bool {
        let mut face_sort = [face[0], face[1], face[2]];
        face_sort.sort();
        self.faces.contains_key(&face_sort)
    }

    pub fn is_tetra_in(&self, tetra: &Tetrahedra) -> bool {
        let mut tetra_sort = [tetra[0], tetra[1], tetra[2], tetra[3]];
        tetra_sort.sort();
        self.tetras.contains(&tetra_sort)
    }

    pub fn count_non_del_halfedges(&self) -> Result<usize> {
        let mut nb_non_del = 0;
        for i in 0..self.topomesh.get_nb_halfedges() {
            let he = self.topomesh.get_halfedge(i)?;
            nb_non_del = nb_non_del + if self.is_edge_in(&he) {0} else {1};
        }
        Ok(nb_non_del)
    }

    pub fn count_non_del_faces(&self) -> Result<usize> {
        let mut nb_non_del = 0;
        for i in 0..self.topomesh.get_nb_faces() {
            let face = self.topomesh.get_face_vertices(i)?;
            nb_non_del = nb_non_del + if self.is_face_in(&face) {0} else {1};
        }
        Ok(nb_non_del)
    }
    
    fn get_opposite_angle(&self, ind_halfedge: usize) -> Result<f32> {
        let ind_next = self.topomesh.get_next_halfedge(ind_halfedge)?;

        let ind_vert1 = self.topomesh.get_halfedge(ind_halfedge)?[0];
        let ind_vert2 = self.topomesh.get_halfedge(ind_halfedge)?[1];
        let ind_vert3 = self.topomesh.get_halfedge(ind_next)?[1];

        let vert1 = self.topomesh.get_vertex(ind_vert1)?;
        let vert2 = self.topomesh.get_vertex(ind_vert2)?;
        let vert3 = self.topomesh.get_vertex(ind_vert3)?;

        let vec31 = vert1 - vert3;
        let vec32 = vert2 - vert3;
        
        let vx = vec31.dot(&vec32);
        let vy = vec31.cross(&vec32).norm();
        let angle = vy.atan2(vx);
        Ok(angle)
    }

    pub fn get_local_non_del_halfedge(&self, shift: Option<usize>) -> Result<Option<(usize, HalfEdge)>>{
        let shift = shift.unwrap_or(0);
        for i in 0..self.topomesh.get_nb_halfedges() {
            let ind_he = (i+shift)%self.topomesh.get_nb_halfedges();
            let he = self.topomesh.get_halfedge(ind_he)?;
            if !self.is_edge_in(&he) {
                let angle1 = self.get_opposite_angle(ind_he)?;
                let ind_opp = self.topomesh.get_opposite_halfedge(ind_he)?;
                let angle2 = self.get_opposite_angle(ind_opp)?;

                if angle1 + angle2 >= std::f32::consts::PI {
                    return Ok(Some((ind_he, he)));
                }
            };
        }
        Ok(None)
    }
    
    pub fn get_non_del_halfedge(&self, shift: Option<usize>) -> Result<Option<(usize, HalfEdge)>>{
        let shift = shift.unwrap_or(0);
        for i in 0..self.topomesh.get_nb_halfedges() {
            let ind_he = (i+shift)%self.topomesh.get_nb_halfedges();
            let he = self.topomesh.get_halfedge(ind_he)?;
            if !self.is_edge_in(&he) {return Ok(Some((ind_he, he)))};
        }
        Ok(None)
    }
    
    pub fn get_non_del_face(&self, shift: Option<usize>) -> Result<Option<(usize, FaceVertices)>>{
        let shift = shift.unwrap_or(0);
        for i in 0..self.topomesh.get_nb_faces() {
            let ind_face = (i+shift)%self.topomesh.get_nb_faces();
            let face = self.topomesh.get_face_vertices(ind_face)?;
            if !self.is_face_in(&face) {return Ok(Some((ind_face, face)))};
        }
        Ok(None)
    }
    
    pub fn get_all_non_del_halfedge(&self) -> Result<Vec<usize>>{
        let mut non_del = Vec::new();
        
        for i in 0..self.topomesh.get_nb_halfedges() {
            let he = self.topomesh.get_halfedge(i)?;
            if !self.is_edge_in(&he) {non_del.push(i);};
        }

        Ok(non_del)
    }
    
    pub fn get_all_non_del_face(&self) -> Result<Vec<usize>>{
        let mut non_del = Vec::new();
        
        for i in 0..self.topomesh.get_nb_faces() {
            let face = self.topomesh.get_face_vertices(i)?;
            if !self.is_face_in(&face) {non_del.push(i);};
        }

        Ok(non_del)
    }
    
    pub fn flip_halfedge(&mut self, ind_halfedge: usize) -> Result<bool>{
        self.topomesh.flip_halfedge(ind_halfedge)
    }

    pub fn split_halfedge(&mut self, vert: &topo_mesh::Vertex, ind_halfedge: usize) -> Result<()> {
        self.topomesh.split_halfedge(vert, ind_halfedge)?;
        self.recompute_struct()
    }

    pub fn split_face(&mut self, vert: &topo_mesh::Vertex, ind_face: usize) -> Result<()> {
        self.topomesh.split_face(vert, ind_face)?;
        self.recompute_struct()
    }
}

