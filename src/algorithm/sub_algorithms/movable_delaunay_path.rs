use anyhow::Result;
use nalgebra::Vector3;
use std::collections::HashSet;

use super::SkeletonInterface3D;

pub struct MovableDelaunayPath<'a, 'b> {
    skeleton_interface: &'b SkeletonInterface3D<'a>,
    ind_palves: Vec<usize>,
    normals: Vec<Vector3<f64>>,
    faces_prev: Vec<Option<[usize; 3]>>,
    faces_prev_prior: Vec<Option<i8>>,
    has_face_connected: Vec<bool>,
}

impl<'a, 'b> MovableDelaunayPath<'a, 'b> {
    fn new(skeleton_interface: &'b SkeletonInterface3D<'a>) -> MovableDelaunayPath<'a, 'b> {
        MovableDelaunayPath {
            skeleton_interface,
            ind_palves: Vec::new(),
            normals: Vec::new(),
            faces_prev: Vec::new(),
            faces_prev_prior: Vec::new(),
            has_face_connected: Vec::new(),
        }
    }

    pub fn print(&self) -> () {
        for i in 0..self.ind_palves.len() {
            let ind_palve = self.ind_palves[i];
            let palve = self
                .skeleton_interface
                .get_partial_alveola_uncheck(ind_palve);
            let seg = palve.alveola().delaunay_segment();
            let seg = if seg[0] == palve.corner() {
                seg
            } else {
                [seg[1], seg[0]]
            };

            let strfac = if let Some([v1, v2, v3]) = self.faces_prev[i] {
                format!("({}, {}, {})", v1, v2, v3)
            } else {
                format!("N")
            };
            let hasfac = if self.has_face_connected[i] { "1" } else { "0" };
            print!("[{} -> {}, {}, {}], ", seg[0], seg[1], strfac, hasfac);
        }
        println!("");
    }

    fn compute_face_prev(&mut self, ind: usize) -> Result<()> {
        let ind_prev = (ind - 1 + self.ind_palves.len()) % self.ind_palves.len();
        let ind_palve_prev = self.ind_palves[ind_prev];
        let palve_prev = self
            .skeleton_interface
            .get_partial_alveola(ind_palve_prev)?;
        let ind_vertex = palve_prev.corner();
        let palve = self
            .skeleton_interface
            .get_partial_alveola(self.ind_palves[ind])?;
        let seg = palve.alveola().delaunay_segment();
        let seg = if seg[0] == palve.corner() {
            seg
        } else {
            [seg[1], seg[0]]
        };

        self.faces_prev_prior[ind] = None;
        self.faces_prev[ind] = None;
        if ind_vertex != seg[0] && ind_vertex != seg[1] {
            let mut tri = [ind_vertex, seg[0], seg[1]];
            tri.sort();
            for pedge in palve.partial_edges() {
                let tri = pedge.edge().delaunay_triangle();
                if tri.contains(&ind_vertex) {
                    let vert1 = self
                        .skeleton_interface
                        .get_mesh()
                        .get_vertex(ind_vertex)?
                        .vertex();
                    let vert2 = self
                        .skeleton_interface
                        .get_mesh()
                        .get_vertex(seg[0])?
                        .vertex();
                    let vert3 = self
                        .skeleton_interface
                        .get_mesh()
                        .get_vertex(seg[1])?
                        .vertex();
                    let nor = (vert2 - vert1).cross(&(vert3 - vert2)).normalize();
                    if nor.dot(&self.normals[ind]) < 0.0 {
                        self.faces_prev_prior[ind] = Some(3);
                    } else if !self.has_face_connected[ind] || !self.has_face_connected[ind_prev] {
                        self.faces_prev_prior[ind] = Some(0);
                    } else if self
                        .skeleton_interface
                        .get_mesh()
                        .is_face_in(ind_vertex, seg[0], seg[1])
                        .is_some()
                    {
                        self.faces_prev_prior[ind] = Some(2);
                    } else {
                        self.faces_prev_prior[ind] = Some(1);
                    }
                    let tri = [ind_vertex, seg[0], seg[1]];
                    self.faces_prev[ind] = Some(tri);
                    return Ok(());
                }
            }
        }
        Ok(())
    }

    pub fn create(
        skeleton_interface: &'b SkeletonInterface3D<'a>,
        ind_palves: Vec<usize>,
        unfaced_hedges: &HashSet<[usize; 2]>,
    ) -> Result<MovableDelaunayPath<'a, 'b>> {
        let mut path = MovableDelaunayPath::new(skeleton_interface);
        for ind in 0..ind_palves.len() {
            path.ind_palves.push(ind_palves[ind]);

            let palve = skeleton_interface.get_partial_alveola(ind_palves[ind])?;

            // counting surrouding faces
            let seg = palve.alveola().delaunay_segment();
            let (seg, seg_opp) = if seg[0] == palve.corner() {
                (seg, [seg[1], seg[0]])
            } else {
                ([seg[1], seg[0]], seg)
            };
            let hedge = skeleton_interface
                .get_mesh()
                .is_edge_in(seg[0], seg[1])
                .unwrap();
            let [v1, v2, v3] = hedge.face().unwrap().vertices();
            let vert1 = v1.vertex();
            let vert2 = v2.vertex();
            let vert3 = v3.vertex();
            let normal = (vert2 - vert1).cross(&(vert3 - vert2)).normalize();
            path.normals.push(normal);

            let mut nb_co = 2;
            if unfaced_hedges.contains(&seg) {
                nb_co = nb_co - 1;
            }
            if unfaced_hedges.contains(&seg_opp) {
                nb_co = nb_co - 1;
            }
            path.has_face_connected.push(nb_co != 0);

            path.faces_prev.push(None);
            path.faces_prev_prior.push(None);
        }

        for ind in 0..ind_palves.len() {
            path.compute_face_prev(ind)?;
        }

        Ok(path)
    }

    pub fn get_couple_to_fusion(&self) -> Option<[usize; 2]> {
        for ind in 0..self.ind_palves.len() {
            if !self.has_face_connected[ind] {
                continue;
            }
            let ind_palve = self.ind_palves[ind];

            let palve = self
                .skeleton_interface
                .get_partial_alveola_uncheck(ind_palve);

            let ind_palve_opp = palve.partial_alveola_opposite().ind();
            let opt_position = self
                .ind_palves
                .iter()
                .position(|&ind_alve| ind_alve == ind_palve_opp);
            if let Some(ind_opp) = opt_position {
                return Some([ind, ind_opp]);
            }
        }
        None
    }

    pub fn fusion_couple(&self, couple: [usize; 2]) -> Result<Vec<MovableDelaunayPath<'a, 'b>>> {
        let mut path1 = MovableDelaunayPath::new(self.skeleton_interface);
        let mut path2 = MovableDelaunayPath::new(self.skeleton_interface);
        let (ind_min, ind_max) = if couple[0] < couple[1] {
            (couple[0], couple[1])
        } else {
            (couple[1], couple[0])
        };
        for i in 0..self.ind_palves.len() {
            if i < ind_min || i > ind_max {
                path1.ind_palves.push(self.ind_palves[i]);
                path1.normals.push(self.normals[i]);
                path1.faces_prev.push(self.faces_prev[i]);
                path1.faces_prev_prior.push(self.faces_prev_prior[i]);
                path1.has_face_connected.push(self.has_face_connected[i]);
            }
            if i > ind_min && i < ind_max {
                path2.ind_palves.push(self.ind_palves[i]);
                path2.normals.push(self.normals[i]);
                path2.faces_prev.push(self.faces_prev[i]);
                path2.faces_prev_prior.push(self.faces_prev_prior[i]);
                path2.has_face_connected.push(self.has_face_connected[i]);
            }
        }

        let mut res = Vec::new();
        if !path1.closed() {
            for i in 0..path1.ind_palves.len() {
                path1.compute_face_prev(i)?;
            }
            // path1.compute_face_prev(ind_min % path1.ind_palves.len())?;
            // path1.compute_face_prev((ind_min + 1) % path1.ind_palves.len())?;
            res.push(path1);
        }
        if !path2.closed() {
            for i in 0..path2.ind_palves.len() {
                path2.compute_face_prev(i)?;
            }
            // path2.compute_face_prev(0)?;
            // path2.compute_face_prev(path2.ind_palves.len() - 1)?;
            res.push(path2);
        }
        Ok(res)
    }

    pub fn get_ind_to_expand(&self) -> Option<usize> {
        let mut ind_exp = None;
        let mut prior = 4;
        for ind in 0..self.ind_palves.len() {
            if let Some(prior_cur) = self.faces_prev_prior[ind] {
                if prior_cur < prior {
                    prior = prior_cur;
                    ind_exp = Some(ind);
                }
            }
        }
        ind_exp
    }

    pub fn expand_ind(
        &mut self,
        ind_exp: usize,
        closing_faces: &mut Vec<[usize; 3]>,
    ) -> Result<()> {
        let face_cur = self.faces_prev[ind_exp]
            .ok_or(anyhow::Error::msg("No previous face associated to edge"))?;
        closing_faces.push(face_cur);

        let new_seg = if face_cur[0] < face_cur[2] {
            [face_cur[0], face_cur[2]]
        } else {
            [face_cur[2], face_cur[0]]
        };
        let &ind_alveola =
            self.skeleton_interface
                .del_seg
                .get(&new_seg)
                .ok_or(anyhow::Error::msg(format!(
                    "Partial alveola ({}, {}) not in skeleton",
                    new_seg[0], new_seg[1]
                )))?;
        let alve = self.skeleton_interface.get_alveola_uncheck(ind_alveola);
        let ind_palve = if alve.partial_alveolae()[0].corner() == face_cur[0] {
            alve.partial_alveolae()[0].ind()
        } else {
            alve.partial_alveolae()[1].ind()
        };
        self.ind_palves[ind_exp] = ind_palve;
        self.has_face_connected[ind_exp] = true;

        let vert1 = self
            .skeleton_interface
            .get_mesh()
            .get_vertex(face_cur[0])?
            .vertex();
        let vert2 = self
            .skeleton_interface
            .get_mesh()
            .get_vertex(face_cur[1])?
            .vertex();
        let vert3 = self
            .skeleton_interface
            .get_mesh()
            .get_vertex(face_cur[2])?
            .vertex();
        let normal = (vert2 - vert1).cross(&(vert3 - vert2)).normalize();
        self.normals[ind_exp] = normal;
        if ind_exp > 0 {
            self.ind_palves.remove(ind_exp - 1);
            self.normals.remove(ind_exp - 1);
            self.has_face_connected.remove(ind_exp - 1);
            self.faces_prev.remove(ind_exp - 1);
            self.faces_prev_prior.remove(ind_exp - 1);
            // self.compute_face_prev(ind_exp - 1)?;
            // self.compute_face_prev(ind_exp % self.ind_palves.len())?;
        } else {
            self.ind_palves.pop();
            self.normals.pop();
            self.has_face_connected.pop();
            self.faces_prev.pop();
            self.faces_prev_prior.pop();
            // self.compute_face_prev(self.ind_palves.len() - 1)?;
            // self.compute_face_prev(0)?;
        }
        for i in 0..self.ind_palves.len() {
            self.compute_face_prev(i)?;
        }
        Ok(())
    }

    pub fn closed(&self) -> bool {
        self.ind_palves.is_empty()
    }
}
