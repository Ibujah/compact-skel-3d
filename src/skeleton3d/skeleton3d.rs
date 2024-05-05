use anyhow::Result;
use nalgebra::base::*;
use std::collections::HashMap;

use crate::geometry::geometry_operations;

#[derive(Copy, Clone)]
/// Sphere
pub struct Sphere {
    /// Sphere center
    pub center: Vector3<f64>,
    /// sphere radius
    pub radius: f64,
}

#[derive(Clone)]
/// 3D Skeleton structure
pub struct Skeleton3D {
    pub(super) nodes: HashMap<usize, Sphere>,
    pub(super) boundary_inds: HashMap<usize, [usize; 4]>,
    pub(super) lone_edges: HashMap<usize, [usize; 2]>, // connects two nodes
    pub(super) edges_on_alv: HashMap<usize, [usize; 2]>, // connects two nodes
    pub(super) edges_on_alv_rev: HashMap<[usize; 2], usize>, // connects two nodes
    pub(super) edges_alveolae: HashMap<usize, Vec<usize>>, // alveolae containing edges
    pub(super) alveolae: HashMap<usize, Vec<usize>>,   // ordered list of nodes
    pub(super) alveolae_edges: HashMap<usize, Vec<usize>>, // edges composing alveola

    pub(super) labels: HashMap<usize, Option<usize>>, // alveolae labels
}

impl Skeleton3D {
    /// Skeleton 3D constructor
    pub fn new() -> Skeleton3D {
        Skeleton3D {
            nodes: HashMap::new(),
            boundary_inds: HashMap::new(),
            lone_edges: HashMap::new(),
            edges_on_alv: HashMap::new(),
            edges_on_alv_rev: HashMap::new(),
            edges_alveolae: HashMap::new(),
            alveolae: HashMap::new(),
            alveolae_edges: HashMap::new(),
            labels: HashMap::new(),
        }
    }

    /// Adds a node to the skeleton
    pub fn add_node(&mut self, ind_node: usize, boundary_points: [Vector3<f64>; 4], boundary_inds: [usize; 4]) -> Result<()> {
        if !self.nodes.contains_key(&ind_node) {
            let (center, radius) = geometry_operations::center_and_radius(boundary_points, None)
                .ok_or(anyhow::Error::msg("Flat tetrahedron"))?;
            let sphere = Sphere { center, radius };
            self.nodes.insert(ind_node, sphere);
            self.boundary_inds.insert(ind_node, boundary_inds);
        }
        Ok(())
    }

    /// Adds a node to the skeleton
    pub fn add_sphere(&mut self, ind_node: usize, sphere: Sphere) -> Result<()> {
        if !self.nodes.contains_key(&ind_node) {
            self.nodes.insert(ind_node, sphere);
        }
        Ok(())
    }

    /// Get nodes hashmap
    pub fn get_nodes(&self) -> &HashMap<usize, Sphere> {
        &self.nodes
    }

    /// Get lone edges hashmap
    pub fn get_lone_edges(&self) -> &HashMap<usize, [usize; 2]> {
        &self.lone_edges
    }

    /// Get edges on alveola hashmap
    pub fn get_edges_on_alv(&self) -> &HashMap<usize, [usize; 2]> {
        &self.edges_on_alv
    }

    /// Get alveolae per edges hashmap
    pub fn get_edges_alveolae(&self) -> &HashMap<usize, Vec<usize>> {
        &self.edges_alveolae
    }

    /// Get nodes per alveolae hashmap
    pub fn get_alveolae(&self) -> &HashMap<usize, Vec<usize>> {
        &self.alveolae
    }

    /// Get edges per alveolae hashmap
    pub fn get_alveolae_edges(&self) -> &HashMap<usize, Vec<usize>> {
        &self.alveolae_edges
    }

    /// Get labels per alveolae hashmap
    pub fn get_labels(&self) -> &HashMap<usize, Option<usize>> {
        &self.labels
    }

    /// Adds an edge to the skeleton
    pub fn add_lone_edge(&mut self, ind_edge: usize, ind_nodes: [usize; 2]) -> () {
        if !self.lone_edges.contains_key(&ind_edge) {
            self.lone_edges.insert(ind_edge, ind_nodes);
            self.edges_alveolae.insert(ind_edge, Vec::new());
        }
    }

    /// Adds an edge to the skeleton
    fn add_alv_edge(&mut self, ind_nodes: [usize; 2]) -> usize {
        if let Some(&ind_edge) = self.edges_on_alv_rev.get(&ind_nodes) {
            ind_edge
        } else {
            let ind_edge = self.edges_on_alv.len();
            self.edges_on_alv.insert(ind_edge, ind_nodes);
            self.edges_on_alv_rev.insert(ind_nodes, ind_edge);
            self.edges_alveolae.insert(ind_edge, Vec::new());
            ind_edge
        }
    }

    /// Adds an alveola to the skeleton
    pub fn add_alveola(&mut self, ind_alveola: usize, ind_nodes: Vec<usize>) -> () {
        if !self.alveolae.contains_key(&ind_alveola) {
            let mut vec_edg = Vec::new();
            for i in 0..ind_nodes.len() {
                let ind1 = ind_nodes[i];
                let ind2 = ind_nodes[(i + 1) % ind_nodes.len()];
                let nods = if ind1 < ind2 {
                    [ind1, ind2]
                } else {
                    [ind2, ind1]
                };
                let ind_edg = self.add_alv_edge(nods);

                vec_edg.push(ind_edg);
                self.edges_alveolae
                    .get_mut(&ind_edg)
                    .unwrap()
                    .push(ind_alveola);
            }
            self.alveolae.insert(ind_alveola, ind_nodes);
            self.alveolae_edges.insert(ind_alveola, vec_edg);
            self.labels.insert(ind_alveola, None);
        }
    }

    /// Assignate a label to a given alveola
    pub fn set_label(&mut self, ind_alveola: usize, label: usize) -> Option<usize> {
        if let Some(l) = self.labels.get_mut(&ind_alveola) {
            let prev = *l;
            *l = Some(label);
            return prev;
        }
        return None;
    }
}
