use std::collections::HashSet;

use anyhow::Result;

use super::skeleton_interface;
use super::SkeletonInterface3D;
use super::SkeletonSingularPath;

pub enum State {
    Computing,
    Closed,
}

pub struct SkeletonBoundaryPath {
    components: Vec<usize>,
    opt_ind_pedge_first: Option<usize>,
    opt_ind_pedge_last: Option<usize>,
}

impl SkeletonBoundaryPath {
    pub fn create(
        ind_pedge: usize,
        skeleton_interface: &mut SkeletonInterface3D,
    ) -> Result<SkeletonBoundaryPath> {
        let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge);
        if !pedge.is_boundary() {
            return Err(anyhow::Error::msg("Partial edge should be boundary"));
        }
        let edge = pedge.edge();
        if !edge.is_computed() {
            let ind_edge = edge.ind();
            skeleton_interface.propagate_edge(ind_edge)?;
        }
        let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge);

        let ind_pedge_next = pedge.partial_edge_next().unwrap().ind();
        let ind_pedge_prev = pedge.partial_edge_prev().unwrap().ind();

        Ok(SkeletonBoundaryPath {
            components: vec![ind_pedge],
            opt_ind_pedge_first: Some(ind_pedge_prev),
            opt_ind_pedge_last: Some(ind_pedge_next),
        })
    }

    pub fn is_valid(&self, skeleton_interface: &SkeletonInterface3D) -> bool {
        for &ind_pedge in self.components.iter() {
            let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge);
            if pedge.edge().degree() != 1 || pedge.partial_alveola().alveola().label().is_none() {
                return false;
            }
        }
        true
    }

    pub fn components(&self) -> &Vec<usize> {
        &self.components
    }

    pub fn append_last(&mut self, skeleton_interface: &mut SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_last) = self.opt_ind_pedge_last {
            if let Some(&ind_pedge_first) = self.components.first() {
                if ind_pedge_last == ind_pedge_first {
                    self.opt_ind_pedge_last = None;
                    self.opt_ind_pedge_first = None;
                    return Ok(State::Closed);
                }
            }
            let edge = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_last)
                .edge();
            if !edge.is_computed() {
                let ind_edge = edge.ind();
                skeleton_interface.propagate_edge(ind_edge)?;
            }
            self.components.push(ind_pedge_last);

            let ind_pedge_next = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_last)
                .partial_edge_next()
                .unwrap()
                .ind();

            self.opt_ind_pedge_last = Some(ind_pedge_next);

            self.check_reach_end(&skeleton_interface)
        } else {
            Ok(State::Closed)
        }
    }

    pub fn rotate_last(&mut self, skeleton_interface: &SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_last) = self.opt_ind_pedge_last {
            let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge_last);
            let pedge_neigh = pedge.partial_edge_neighbor();

            let alve = pedge_neigh.partial_alveola().alveola();
            let pedge_next = if alve.is_full() {
                pedge_neigh.partial_edge_next().unwrap()
            } else {
                pedge_neigh.partial_edge_opposite()
            };
            let ind_pedge_new = pedge_next.ind();
            self.opt_ind_pedge_last = Some(ind_pedge_new);
            self.check_reach_end(skeleton_interface)
        } else {
            Ok(State::Closed)
        }
    }

    fn check_reach_end(&mut self, skeleton_interface: &SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_last) = self.opt_ind_pedge_last {
            let ind_pedge_opp = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_last)
                .partial_edge_opposite()
                .ind();

            if let Some(&ind_pedge) = self.components.last() {
                if ind_pedge == ind_pedge_opp {
                    self.opt_ind_pedge_last = None;
                }
            }
        }

        if self.opt_ind_pedge_last.is_none() {
            Ok(State::Closed)
        } else {
            Ok(State::Computing)
        }
    }

    pub fn append_first(&mut self, skeleton_interface: &mut SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_first) = self.opt_ind_pedge_first {
            let edge = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_first)
                .edge();
            if !edge.is_computed() {
                let ind_edge = edge.ind();
                skeleton_interface.propagate_edge(ind_edge)?;
            }
            self.components.insert(0, ind_pedge_first);

            let ind_pedge_prev = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_first)
                .partial_edge_prev()
                .unwrap()
                .ind();

            self.opt_ind_pedge_first = Some(ind_pedge_prev);

            self.check_reach_beg(&skeleton_interface)
        } else {
            Ok(State::Closed)
        }
    }

    pub fn rotate_first(&mut self, skeleton_interface: &SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_first) = self.opt_ind_pedge_first {
            let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge_first);
            let pedge_neigh = pedge.partial_edge_neighbor();

            let alve = pedge_neigh.partial_alveola().alveola();
            let pedge_prev = if alve.is_full() {
                pedge_neigh.partial_edge_prev().unwrap()
            } else {
                pedge_neigh.partial_edge_opposite()
            };
            let ind_pedge_new = pedge_prev.ind();
            self.opt_ind_pedge_first = Some(ind_pedge_new);
            self.check_reach_beg(skeleton_interface)
        } else {
            Ok(State::Closed)
        }
    }

    fn check_reach_beg(&mut self, skeleton_interface: &SkeletonInterface3D) -> Result<State> {
        if let Some(ind_pedge_first) = self.opt_ind_pedge_first {
            let ind_pedge_opp = skeleton_interface
                .get_partial_edge_uncheck(ind_pedge_first)
                .partial_edge_opposite()
                .ind();

            if let Some(&ind_pedge) = self.components.first() {
                if ind_pedge_opp == ind_pedge {
                    self.opt_ind_pedge_first = None;
                }
            }
        }

        if self.opt_ind_pedge_first.is_none() {
            Ok(State::Closed)
        } else {
            Ok(State::Computing)
        }
    }

    pub fn ind_last_partial_edge(&self) -> Option<usize> {
        self.opt_ind_pedge_last
    }

    pub fn ind_first_partial_edge(&self) -> Option<usize> {
        self.opt_ind_pedge_first
    }

    pub fn follow_boundary_path(
        &mut self,
        skeleton_interface: &mut SkeletonInterface3D,
    ) -> Result<()> {
        loop {
            if let Some(ind_pedge_last) = self.ind_last_partial_edge() {
                let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge_last);
                if pedge.is_boundary() {
                    self.append_last(skeleton_interface)?;
                } else {
                    self.rotate_last(skeleton_interface)?;
                }
            } else {
                break;
            }
        }
        loop {
            if let Some(ind_pedge_first) = self.ind_first_partial_edge() {
                let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge_first);
                if pedge.is_boundary() {
                    self.append_first(skeleton_interface)?;
                } else {
                    self.rotate_first(skeleton_interface)?;
                }
            } else {
                break;
            }
        }
        Ok(())
    }

    pub fn most_salient(&self, skeleton_interface: &SkeletonInterface3D) -> Result<Option<usize>> {
        let mut sin_min = 1.0;
        let mut ind_min = None;
        // for i in 0..self.components.len() {
        //     let pedge_cur = skeleton_interface.get_partial_edge_uncheck(self.components[i]);
        //     println!(
        //         "({}, {}) {} -> {}",
        //         pedge_cur.partial_alveola().alveola().label().unwrap_or(0),
        //         pedge_cur.edge().degree(),
        //         pedge_cur.partial_node_first().unwrap().node().ind(),
        //         pedge_cur.partial_node_last().unwrap().node().ind(),
        //     );
        // }
        // println!("");
        for i in 0..(self.components.len() - 1) {
            let pedge_cur = skeleton_interface.get_partial_edge_uncheck(self.components[i]);
            if pedge_cur.edge().degree() != 1 {
                continue;
            }
            let pedge_aft = skeleton_interface.get_partial_edge_uncheck(self.components[i + 1]);
            let pedge_nex = pedge_cur.partial_edge_next().unwrap();
            let (vert1, _) = pedge_cur
                .partial_node_first()
                .unwrap()
                .node()
                .center_and_radius()?;
            let (vert2, _) = pedge_cur
                .partial_node_last()
                .unwrap()
                .node()
                .center_and_radius()?;
            let (vert3, _) = pedge_aft
                .partial_node_last()
                .unwrap()
                .node()
                .center_and_radius()?;
            let (vert4, _) = pedge_nex
                .partial_node_last()
                .unwrap()
                .node()
                .center_and_radius()?;

            let vec_cur = (vert2 - vert1).normalize();
            let vec_nex = vert4 - vert2;
            let nor = vec_cur.cross(&vec_nex).normalize();
            let vec_aft = (vert3 - vert2).normalize();
            let vec_cross = vec_cur.cross(&vec_aft);
            if vec_cross.dot(&nor) > 0.0 {
                let sin_ang = vec_cross.norm();
                if sin_ang < sin_min {
                    ind_min = Some(i);
                    sin_min = sin_ang;
                }
            }
        }
        Ok(ind_min)
    }

    pub fn separate_path_at(
        &self,
        ind_in_path: usize,
        skeleton_interface: &SkeletonInterface3D,
    ) -> Result<Vec<SkeletonBoundaryPath>> {
        let mut vec1 = Vec::new();
        for i in 0..(ind_in_path + 1) {
            let ind_pedge = self.components[i];
            let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge);
            if pedge.edge().degree() != 1 || pedge.partial_alveola().alveola().label().is_none() {
                break;
            }
            vec1.push(ind_pedge);
        }
        let mut vec2 = Vec::new();
        for i in (ind_in_path + 1)..self.components.len() {
            let ind_pedge = self.components[i];
            let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge);
            if pedge.edge().degree() != 1 || pedge.partial_alveola().alveola().label().is_none() {
                break;
            }
            vec2.push(ind_pedge);
        }

        let mut vec_res = Vec::new();
        if !vec1.is_empty() {
            let path1 = SkeletonBoundaryPath {
                components: vec1,
                opt_ind_pedge_first: None,
                opt_ind_pedge_last: None,
            };
            if path1.is_valid(skeleton_interface) {
                vec_res.push(path1);
            } else {
                return Err(anyhow::Error::msg("New boundary path 1 is not valid"));
            }
        }
        if !vec2.is_empty() {
            let path2 = SkeletonBoundaryPath {
                components: vec2,
                opt_ind_pedge_first: None,
                opt_ind_pedge_last: None,
            };
            if path2.is_valid(skeleton_interface) {
                vec_res.push(path2);
            } else {
                return Err(anyhow::Error::msg("New boundary path 2 is not valid"));
            }
        }
        Ok(vec_res)
    }

    pub fn create_path_excluding(
        &self,
        skeleton_interface: &mut SkeletonInterface3D,
        ind_in_path: usize,
    ) -> Result<Option<(SkeletonSingularPath, Vec<usize>, Vec<usize>)>> {
        if ind_in_path >= self.components.len() - 1 {
            return Err(anyhow::Error::msg("Index in path out of bounds"));
        }
        println!("");
        let mut ind_first = ind_in_path;
        let mut ind_last = ind_in_path;
        let mut vec_ind_pedges = Vec::new();
        let mut vec_ind_removed_alveola = Vec::new();
        loop {
            let pedge = skeleton_interface.get_partial_edge_uncheck(self.components[ind_first]);

            if pedge.partial_edge_prev().unwrap().is_boundary() {
                if ind_first == 0 {
                    return Ok(None);
                }
                ind_first = ind_first - 1;
            } else {
                if ind_first != 0 {
                    vec_ind_pedges.push(self.components[ind_first - 1]);
                }
                break;
            }
        }
        loop {
            let pedge = skeleton_interface.get_partial_edge_uncheck(self.components[ind_last]);

            if pedge.partial_edge_next().unwrap().is_boundary() {
                if ind_last >= self.components.len() - 1 {
                    return Ok(None);
                }
                ind_last = ind_last + 1;
            } else {
                if ind_last < self.components.len() - 1 {
                    vec_ind_pedges.push(self.components[ind_last + 1]);
                }
                break;
            }
        }
        let ind_alve1 = skeleton_interface
            .get_partial_edge_uncheck(self.components[ind_first])
            .partial_alveola()
            .alveola()
            .ind();
        let ind_alve2 = skeleton_interface
            .get_partial_edge_uncheck(self.components[ind_last])
            .partial_alveola()
            .alveola()
            .ind();

        let alve1 = skeleton_interface
            .get_partial_edge_uncheck(self.components[ind_first])
            .partial_alveola()
            .alveola();
        let alve2 = skeleton_interface
            .get_partial_edge_uncheck(self.components[ind_last])
            .partial_alveola()
            .alveola();
        println!("{} {}", ind_first, ind_last);
        println!(
            "({}, {}) ({}, {})",
            alve1.ind(),
            alve1.label().unwrap_or(0),
            alve2.ind(),
            alve2.label().unwrap_or(0)
        );

        let mut edges_path = HashSet::new();
        if ind_alve1 == ind_alve2 {
            vec_ind_removed_alveola.push(ind_alve1);
            let edges1: HashSet<usize> = skeleton_interface
                .get_alveola_uncheck(ind_alve1)
                .edges()
                .iter()
                .filter_map(|&edge| {
                    if edge.degree() != 1 {
                        Some(edge.ind())
                    } else {
                        None
                    }
                })
                .collect();
            println!("edges1 {:?}", edges1);
            for ind_edge in edges1 {
                edges_path.insert(ind_edge);
            }
        } else {
            vec_ind_removed_alveola.push(ind_alve1);
            vec_ind_removed_alveola.push(ind_alve2);
            let edges1: HashSet<usize> = skeleton_interface
                .get_alveola_uncheck(ind_alve1)
                .edges()
                .iter()
                .filter_map(|&edge| {
                    if edge.degree() != 1 {
                        Some(edge.ind())
                    } else {
                        None
                    }
                })
                .collect();
            let edges2: HashSet<usize> = skeleton_interface
                .get_alveola_uncheck(ind_alve2)
                .edges()
                .iter()
                .filter_map(|&edge| {
                    if edge.degree() != 1 {
                        Some(edge.ind())
                    } else {
                        None
                    }
                })
                .collect();
            println!("edges1 {:?}", edges1);
            println!("edges2 {:?}", edges1);
            for &ind_edge in edges1.iter() {
                if !edges2.contains(&ind_edge) {
                    edges_path.insert(ind_edge);
                }
            }
            for &ind_edge in edges2.iter() {
                if !edges1.contains(&ind_edge) {
                    edges_path.insert(ind_edge);
                }
            }
        }

        let ind_pedge_first = skeleton_interface
            .get_partial_edge_uncheck(self.components[ind_last])
            .partial_edge_next()
            .unwrap()
            .ind();
        let ind_pedge_last = skeleton_interface
            .get_partial_edge_uncheck(self.components[ind_first])
            .partial_edge_prev()
            .unwrap()
            .ind();

        let mut sing_path = SkeletonSingularPath::create(ind_pedge_first);
        sing_path.append_last(skeleton_interface)?;

        println!("{:?}", edges_path);
        for &ind_edge in edges_path.iter() {
            let edge = skeleton_interface.get_edge_uncheck(ind_edge);
            let nods: Vec<usize> = edge.nodes().iter().map(|&nod| nod.ind()).collect();
            println!("{}: {:?}", edge.ind(), nods);
        }
        print!("from ");
        let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge_first);
        println!(
            "({}, {}): {} -> {}",
            pedge.ind(),
            pedge.edge().degree(),
            pedge.partial_node_first().unwrap().node().ind(),
            pedge.partial_node_last().unwrap().node().ind()
        );
        print!("to ");
        let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge_last);
        println!(
            "({}, {}): {} -> {}",
            pedge.ind(),
            pedge.edge().degree(),
            pedge.partial_node_first().unwrap().node().ind(),
            pedge.partial_node_last().unwrap().node().ind()
        );

        loop {
            if let Some(ind_pedge_last) = sing_path.ind_last_partial_edge() {
                let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge_last);
                if edges_path.contains(&pedge.edge().ind()) {
                    sing_path.append_last(skeleton_interface)?;
                } else {
                    sing_path.rotate_last(skeleton_interface)?;
                }
            } else {
                break;
            }
        }

        Ok(Some((sing_path, vec_ind_pedges, vec_ind_removed_alveola)))
    }
}

pub fn next_boundary_edge(
    ind_pedge: usize,
    skeleton_interface: &SkeletonInterface3D,
) -> Option<usize> {
    let pedge_cur = skeleton_interface.get_partial_edge_uncheck(ind_pedge);
    if pedge_cur.edge().degree() != 1 {
        return None;
    }
    if pedge_cur.partial_alveola().alveola().label().is_none() {
        return None;
    }
    let pnode_last = pedge_cur.partial_node_last().unwrap();
    let pedge_next: Vec<usize> = pnode_last
        .partial_edge_next()
        .iter()
        .filter_map(|pedge| {
            if pedge.partial_alveola().alveola().label().is_some()
                && pedge.edge().degree() == 1
                && pedge.edge().ind() != pedge_cur.edge().ind()
            {
                Some(pedge.ind())
            } else {
                None
            }
        })
        .collect();

    if pedge_next.len() != 1 {
        return None;
    }

    Some(pedge_next[0])
}

pub fn compute_saliency(
    ind_pedge: usize,
    skeleton_interface: &SkeletonInterface3D,
) -> Result<Option<f32>> {
    if let Some(ind_aft) = next_boundary_edge(ind_pedge, skeleton_interface) {
        let pedge_cur = skeleton_interface.get_partial_edge_uncheck(ind_pedge);
        let pedge_aft = skeleton_interface.get_partial_edge_uncheck(ind_aft);

        let pedge_nex = pedge_cur.partial_edge_next().unwrap();
        let (vert1, _) = pedge_cur
            .partial_node_first()
            .unwrap()
            .node()
            .center_and_radius()?;
        let (vert2, _) = pedge_cur
            .partial_node_last()
            .unwrap()
            .node()
            .center_and_radius()?;
        let (vert3, _) = pedge_aft
            .partial_node_last()
            .unwrap()
            .node()
            .center_and_radius()?;
        let (vert4, _) = pedge_nex
            .partial_node_last()
            .unwrap()
            .node()
            .center_and_radius()?;

        let vec_cur = (vert2 - vert1).normalize();
        let vec_nex = vert4 - vert2;
        let nor = vec_cur.cross(&vec_nex).normalize();
        let vec_aft = (vert3 - vert2).normalize();
        let vec_cross = vec_cur.cross(&vec_aft);
        if vec_cross.dot(&nor) > 0.0 {
            let sin_ang = vec_cross.norm();
            return Ok(Some(sin_ang));
        }
    }
    Ok(None)
}

pub fn excluded_alveolae(
    ind_pedge: usize,
    skeleton_interface: &SkeletonInterface3D,
) -> HashSet<usize> {
    let pnode = skeleton_interface
        .get_partial_edge_uncheck(ind_pedge)
        .partial_node_last()
        .unwrap();
    pnode
        .partial_edge_next()
        .iter()
        .filter_map(|pedge| {
            let alve = pedge.partial_alveola().alveola();
            if alve.is_full() {
                Some(alve.ind())
            } else {
                None
            }
        })
        .collect()
}

pub fn singular_path_to_exclude_alveolae(
    set_alve_to_exclude: &HashSet<usize>,
    skeleton_interface: &mut SkeletonInterface3D,
) -> Result<Option<SkeletonSingularPath>> {
    let mut opt_ind_first_pedge = None;
    for &ind_alve in set_alve_to_exclude {
        let palve = skeleton_interface
            .get_alveola_uncheck(ind_alve)
            .partial_alveolae()[0];
        for pedge in palve.partial_edges().iter() {
            if pedge.edge().degree() == 1 {
                continue;
            }
            let neigh_alve: HashSet<usize> = pedge
                .edge()
                .alveolae()
                .iter()
                .filter_map(|alve| {
                    if set_alve_to_exclude.contains(&alve.ind()) {
                        Some(alve.ind())
                    } else {
                        None
                    }
                })
                .collect();
            if neigh_alve.len() == 1 {
                opt_ind_first_pedge = Some(pedge.ind());
                break;
            }
        }
        if opt_ind_first_pedge.is_some() {
            break;
        }
    }

    if let Some(ind_first_pedge) = opt_ind_first_pedge {
        let mut sing_path = SkeletonSingularPath::create(ind_first_pedge);
        sing_path.append_last(skeleton_interface)?;

        loop {
            if let Some(ind_pedge_last) = sing_path.ind_last_partial_edge() {
                let pedge = skeleton_interface.get_partial_edge_uncheck(ind_pedge_last);
                let should_append = if pedge.edge().degree() == 1 {
                    false
                } else {
                    let neigh_alve: HashSet<usize> = pedge
                        .edge()
                        .alveolae()
                        .iter()
                        .filter_map(|alve| {
                            if set_alve_to_exclude.contains(&alve.ind()) {
                                Some(alve.ind())
                            } else {
                                None
                            }
                        })
                        .collect();
                    if neigh_alve.len() == 1 {
                        true
                    } else {
                        false
                    }
                };

                if should_append {
                    sing_path.append_last(skeleton_interface)?;
                } else {
                    sing_path.rotate_last(skeleton_interface)?;
                }
            } else {
                break;
            }
        }

        Ok(Some(sing_path))
    } else {
        Ok(None)
    }
}
