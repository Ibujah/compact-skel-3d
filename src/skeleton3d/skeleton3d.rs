use nalgebra::base::*;
use std::collections::HashMap;

#[derive(Copy, Clone)]
pub struct Sphere {
    pub center: Vector3<f32>,
    pub radius: f32,
}

pub struct Skeleton3D {
    nodes: HashMap<usize, Sphere>,
    edges: HashMap<usize, [usize; 2]>,    // connects two nodes
    alveolae: HashMap<usize, Vec<usize>>, // ordered list of edges

    labels: HashMap<usize, Option<usize>>, // alveolae labels
}

impl Skeleton3D {
    pub fn new() -> Skeleton3D {
        Skeleton3D {
            nodes: HashMap::new(),
            edges: HashMap::new(),
            alveolae: HashMap::new(),
            labels: HashMap::new(),
        }
    }

    pub fn add_node(&mut self, ind_node: usize, boundary_points: [Vector3<f32>; 4]) -> () {
        if !self.nodes.contains_key(&ind_node) {
            let sphere = Sphere {
                center: Vector3::new(0.0, 0.0, 0.0),
                radius: 0.0,
            };
            self.nodes.insert(ind_node, sphere);
        }
    }

    pub fn add_edge(&mut self, ind_edge: usize, ind_nodes: [usize; 2]) -> () {
        if !self.edges.contains_key(&ind_edge) {
            self.edges.insert(ind_edge, ind_nodes);
        }
    }

    pub fn add_alveola(&mut self, ind_alveola: usize, ind_edges: Vec<usize>) -> () {
        if !self.alveolae.contains_key(&ind_alveola) {
            self.alveolae.insert(ind_alveola, ind_edges);
            self.labels.insert(ind_alveola, None);
        }
    }

    pub fn set_label(&mut self, ind_alveola: usize, label: usize) -> Option<usize> {
        if let Some(l) = self.labels.get_mut(&ind_alveola) {
            let prev = *l;
            *l = Some(label);
            return prev;
        }
        return None;
    }
}
