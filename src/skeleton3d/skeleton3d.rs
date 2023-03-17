use anyhow::Result;
use nalgebra::base::*;
use std::collections::HashMap;

use crate::geometry::geometry_operations;

#[derive(Copy, Clone)]
pub struct Sphere {
    pub center: Vector3<f32>,
    pub radius: f32,
}

#[derive(Clone)]
pub struct Skeleton3D {
    pub(super) nodes: HashMap<usize, Sphere>,
    pub(super) edges: HashMap<usize, [usize; 2]>, // connects two nodes
    pub(super) alveolae: HashMap<usize, Vec<usize>>, // ordered list of nodes

    pub(super) labels: HashMap<usize, Option<usize>>, // alveolae labels
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

    pub fn add_node(&mut self, ind_node: usize, boundary_points: [Vector3<f32>; 4]) -> Result<()> {
        if !self.nodes.contains_key(&ind_node) {
            let (center, radius) = geometry_operations::center_and_radius(boundary_points, None)
                .ok_or(anyhow::Error::msg("Flat tetrahedron"))?;
            let sphere = Sphere { center, radius };
            self.nodes.insert(ind_node, sphere);
        }
        Ok(())
    }

    pub fn add_edge(&mut self, ind_edge: usize, ind_nodes: [usize; 2]) -> () {
        if !self.edges.contains_key(&ind_edge) {
            self.edges.insert(ind_edge, ind_nodes);
        }
    }

    pub fn add_alveola(&mut self, ind_alveola: usize, ind_nodes: Vec<usize>) -> () {
        if !self.alveolae.contains_key(&ind_alveola) {
            self.alveolae.insert(ind_alveola, ind_nodes);
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
