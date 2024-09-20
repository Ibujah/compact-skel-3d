use anyhow::Result;
use nalgebra::base::*;
use std::collections::{HashMap, HashSet};

use simple_delaunay_lib::delaunay_3d::delaunay_struct_3d::{
    DelaunayStructure3D, ExtendedTetrahedron,
};
use simple_delaunay_lib::delaunay_3d::simplicial_struct_3d::Node;

use crate::mesh3d::mesh_operations;
use crate::mesh3d::{manifold_mesh3d, ManifoldMesh3D};

use crate::geometry::geometry_operations::center_and_radius;

pub type Edge = [usize; 2];
pub type Triangle = [usize; 3];
pub type Tetrahedron = [usize; 4];

/// Delaunay structure
pub struct DelaunayInterface<'a> {
    mesh: &'a mut ManifoldMesh3D,
    del_struct: DelaunayStructure3D,

    vertex_edges: Vec<Vec<(usize, usize)>>,

    non_del_edges: Vec<usize>,
    non_del_faces: Vec<usize>,

    initial_vertices_number: usize,
}

impl<'a> DelaunayInterface<'a> {
    fn generate_struct(&mut self) -> Result<()> {
        let mut points = Vec::new();
        for v in 0..self.mesh.get_nb_vertices() {
            let vert = self.mesh.get_vertex(v)?.vertex();
            points.push([vert[0], vert[1], vert[2]]);
            self.vertex_edges.push(Vec::new());
        }
        self.del_struct.insert_vertices(&points, true)?;

        for ind_tet in 0..self.del_struct.get_simplicial().get_nb_tetrahedra() {
            let tetra = self.del_struct.get_simplicial().get_tetrahedron(ind_tet)?;

            for tri in tetra.halftriangles() {
                let hes = tri.halfedges();
                for i in 0..3 {
                    if let (Node::Value(i1), Node::Value(_)) =
                        (hes[i].first_node(), hes[i].last_node())
                    {
                        self.vertex_edges[i1].push((tri.ind(), i));
                    }
                }
            }
        }

        Ok(())
    }

    /// compute if each triangle is cocone or not
    pub fn compute_cocone(&self) -> Result<HashMap<[usize; 3], bool>> {
        let mut sphere_centers = Vec::new();
        let mut voro_cell = vec![Vec::new(); self.mesh.get_nb_vertices()];
        let mut cell_axis = Vec::new();

        for ind_tet in 0..self.del_struct.get_simplicial().get_nb_tetrahedra() {
            let tetra = self.del_struct.get_simplicial().get_tetrahedron(ind_tet)?;

            let ext_tet = self.del_struct.get_extended_tetrahedron(ind_tet)?;
            let center = match ext_tet {
                ExtendedTetrahedron::Tetrahedron([p1, p2, p3, p4]) => {
                    let pts = [
                        Vector3::new(p1[0], p1[1], p1[2]),
                        Vector3::new(p2[0], p2[1], p2[2]),
                        Vector3::new(p3[0], p3[1], p3[2]),
                        Vector3::new(p4[0], p4[1], p4[2]),
                    ];
                    let (center, _) = center_and_radius(pts, None).unwrap();
                    [center[0], center[1], center[2], 1.0]
                }
                ExtendedTetrahedron::Triangle([p1, p2, p3]) => {
                    let pt1 = Vector3::new(p1[0], p1[1], p1[2]);
                    let pt2 = Vector3::new(p2[0], p2[1], p2[2]);
                    let pt3 = Vector3::new(p3[0], p3[1], p3[2]);
                    let v12 = pt2 - pt1;
                    let v23 = pt3 - pt2;
                    let nor = v12.cross(&v23).normalize();
                    [nor[0], nor[1], nor[2], 0.]
                }
            };
            sphere_centers.push(center);

            for &nod in tetra.nodes().iter() {
                if let Node::Value(ind_v) = nod {
                    voro_cell[ind_v].push(ind_tet);
                }
            }
        }

        for ind_cell in 0..voro_cell.len() {
            let vert = self.mesh.get_vertex(ind_cell)?.vertex();
            let mut axis_opt = None;
            let mut farthest_opt = None;
            let mut dmin_opt = None;
            for sub_ind_sph in 0..voro_cell[ind_cell].len() {
                let ind_sph = voro_cell[ind_cell][sub_ind_sph];
                let center = sphere_centers[ind_sph];
                let ctr = Vector3::new(center[0], center[1], center[2]);
                if center[3] == 0. {
                    axis_opt = Some(ctr);
                } else {
                    let vec = ctr - vert;
                    let dist = vec.norm();
                    if let Some(dmin) = dmin_opt {
                        if dmin > dist {
                            farthest_opt = Some(vec / dist);
                            dmin_opt = Some(dist);
                        }
                    } else {
                        farthest_opt = Some(vec / dist);
                        dmin_opt = Some(dist);
                    };
                }
            }
            let farthest = farthest_opt.unwrap();
            let axis = if let Some(axis) = axis_opt {
                if axis.dot(&farthest) > 0. {
                    -axis
                } else {
                    axis
                }
            } else {
                farthest
            };
            cell_axis.push(axis);
        }

        let mut tri_cocone = HashMap::new();
        for ind_tri in 0..self.del_struct.get_simplicial().get_nb_tetrahedra() * 4 {
            let tri = self.del_struct.get_simplicial().get_halftriangle(ind_tri)?;

            if let [Node::Value(i1), Node::Value(i2), Node::Value(i3)] = tri.nodes() {
                let mut tri_ind = [i1, i2, i3];
                tri_ind.sort();
                if !tri_cocone.contains_key(&tri_ind) {
                    let ind_tet1 = tri.tetrahedron().ind();
                    let ind_tet2 = tri.opposite().tetrahedron().ind();

                    // let mut respect_cocone = self.mesh.is_face_in(i1, i2, i3).is_some();
                    let mut respect_cocone = false;
                    if respect_cocone {
                        for nod in tri.nodes() {
                            if let Node::Value(ind_v) = nod {
                                let p = self.mesh.get_vertex(ind_v)?.vertex();
                                let vp = cell_axis[ind_v];
                                let ctr1h = sphere_centers[ind_tet1];
                                let ctr2h = sphere_centers[ind_tet2];
                                if !self.respects_cocone_conditions(p, vp, ctr1h, ctr2h)? {
                                    respect_cocone = false;
                                    break;
                                }
                            } else {
                                respect_cocone = false;
                                break;
                            }
                        }
                    }
                    tri_cocone.insert(tri_ind, respect_cocone);
                }
            }
        }

        Ok(tri_cocone)
    }

    fn insert_vertex(&mut self, ind_vertex: usize, near_to: usize) -> Result<()> {
        let vert = self.mesh.get_vertex(ind_vertex)?.vertex();
        self.del_struct.insert_vertex(
            [vert[0] as f64, vert[1] as f64, vert[2] as f64],
            Some(near_to),
        )?;
        self.vertex_edges.push(Vec::new());

        let tet_update = self
            .del_struct
            .get_simplicial()
            .get_tetrahedra_containing(&Node::Value(self.vertex_edges.len() - 1));

        let mut vert_to_check = HashSet::new();
        for tetra in tet_update {
            for tri in tetra.halftriangles() {
                let hes = tri.halfedges();
                for i in 0..3 {
                    if let (Node::Value(i1), Node::Value(i2)) =
                        (hes[i].first_node(), hes[i].last_node())
                    {
                        self.vertex_edges[i1].push((tri.ind(), i));
                        vert_to_check.insert(i1);
                        vert_to_check.insert(i2);
                    }
                }
            }
        }

        for iv in vert_to_check {
            self.vertex_edges[iv] = self.vertex_edges[iv]
                .iter()
                .filter_map(|&(it, i)| {
                    if let Ok(tri) = self.del_struct.get_simplicial().get_halftriangle(it) {
                        if tri.halfedges()[i].first_node().equals(&Node::Value(iv)) {
                            Some((it, i))
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                })
                .collect();
        }

        Ok(())
    }

    fn fill_non_del(&mut self) -> () {
        self.non_del_edges.clear();
        self.non_del_faces.clear();

        for ind_fac in 0..self.mesh.get_nb_faces() {
            let face = self.mesh.get_face(ind_fac).unwrap();
            let face_vert = face.vertices_inds();
            if !self.is_face_in(&face_vert) {
                self.non_del_faces.push(ind_fac);
                for he in face.halfedges() {
                    let hedg_vert = he.halfedge();
                    if !self.is_edge_in(&hedg_vert) {
                        self.non_del_edges.push(he.ind());
                    };
                }
            };
        }
    }

    /// Creates Delaunay structure from mesh
    pub fn from_mesh(mesh: &'a mut ManifoldMesh3D) -> Result<DelaunayInterface<'a>> {
        let initial_vertices_number = mesh.get_nb_vertices();
        let mut deltet = DelaunayInterface {
            mesh,
            del_struct: DelaunayStructure3D::new(),
            vertex_edges: Vec::new(),
            non_del_edges: Vec::new(),
            non_del_faces: Vec::new(),
            initial_vertices_number,
        };

        deltet.generate_struct()?;

        deltet.fill_non_del();

        Ok(deltet)
    }

    /// Mesh getter
    pub fn get_mesh(&self) -> &ManifoldMesh3D {
        self.mesh
    }

    /// Tetrahedra set getter
    pub fn get_faces(&self) -> HashMap<Triangle, Vec<Tetrahedron>> {
        let mut face_set = HashMap::new();
        for ind_tet in 0..self.del_struct.get_simplicial().get_nb_tetrahedra() {
            let tetra = self
                .del_struct
                .get_simplicial()
                .get_tetrahedron(ind_tet)
                .unwrap();
            if let [Node::Value(i1), Node::Value(i2), Node::Value(i3), Node::Value(i4)] =
                tetra.nodes()
            {
                let mut tetra_ind = [i1, i2, i3, i4];
                tetra_ind.sort();

                face_set
                    .entry([tetra_ind[0], tetra_ind[1], tetra_ind[2]])
                    .or_insert(Vec::new())
                    .push(tetra_ind);
                face_set
                    .entry([tetra_ind[0], tetra_ind[1], tetra_ind[3]])
                    .or_insert(Vec::new())
                    .push(tetra_ind);
                face_set
                    .entry([tetra_ind[0], tetra_ind[2], tetra_ind[3]])
                    .or_insert(Vec::new())
                    .push(tetra_ind);
                face_set
                    .entry([tetra_ind[1], tetra_ind[2], tetra_ind[3]])
                    .or_insert(Vec::new())
                    .push(tetra_ind);
            }
        }
        face_set
    }

    /// Checks if vertex was an original mesh vertex
    pub fn is_original_vertex(&self, ind_vertex: usize) -> bool {
        ind_vertex < self.initial_vertices_number
    }

    /// Checks if edge is in Delaunay
    fn is_edge_in(&self, edge: &Edge) -> bool {
        self.vertex_edges[edge[0]]
            .iter()
            .position(|&(it, i)| {
                if let Ok(tri) = self.del_struct.get_simplicial().get_halftriangle(it) {
                    tri.halfedges()[i].last_node().equals(&Node::Value(edge[1]))
                } else {
                    false
                }
            })
            .is_some()
    }

    fn respects_cocone_conditions(
        &self,
        p: Vector3<f64>,
        vp: Vector3<f64>,
        ctr1h: [f64; 4],
        ctr2h: [f64; 4],
    ) -> Result<bool> {
        // let p = self.mesh.get_vertex(ind_cell)?.vertex();
        // let vp = self.cell_axis[ind_cell];
        // let ctr1h = self.sphere_centers[ind_sph1];
        // let ctr2h = self.sphere_centers[ind_sph2];

        let a = if ctr1h[3] == 0. {
            Vector3::new(ctr1h[0], ctr1h[1], ctr1h[2]).normalize()
        } else {
            let ctr1 = Vector3::new(ctr1h[0], ctr1h[1], ctr1h[2]);
            (ctr1 - p).normalize()
        };
        let b = if ctr2h[3] == 0. {
            Vector3::new(ctr2h[0], ctr2h[1], ctr2h[2]).normalize()
        } else {
            let ctr2 = Vector3::new(ctr2h[0], ctr2h[1], ctr2h[2]);
            (ctr2 - p).normalize()
        };

        let vp_a = vp.dot(&a);
        let vp_b = vp.dot(&b);
        let cos_3_pi_on_8 = -0.125;

        Ok(vp_a * vp_b < 0. || vp_a < cos_3_pi_on_8 || vp_b < cos_3_pi_on_8)
    }

    /// Checks if face is in Delaunay
    fn is_face_in(&self, face: &Triangle) -> bool {
        self.vertex_edges[face[0]]
            .iter()
            .position(|&(it, i)| {
                if let Ok(tri) = self.del_struct.get_simplicial().get_halftriangle(it) {
                    let he = tri.halfedges()[i];
                    if he.last_node().equals(&Node::Value(face[1])) {
                        he.next().last_node().equals(&Node::Value(face[2]))
                    } else {
                        false
                    }
                } else {
                    false
                }
            })
            .is_some()
    }

    /// Count number of non Delaunay halfedges
    pub fn count_non_del_halfedges(&mut self) -> usize {
        if self.non_del_edges.len() == 0 {
            self.fill_non_del();
        }
        self.non_del_edges.len()
    }

    /// Count number of non Delaunay faces
    pub fn count_non_del_faces(&mut self) -> usize {
        if self.non_del_faces.len() == 0 {
            self.fill_non_del();
        }
        self.non_del_faces.len()
    }

    /// Gets first globally non Delaunay halfedge, starting from a shift
    pub fn get_non_del_halfedge(&mut self) -> Result<Option<manifold_mesh3d::IterHalfEdge>> {
        let mut length: Vec<(usize, f64)> = self
            .non_del_edges
            .iter()
            .map(|&ind_he| {
                if let Ok(edge) = self.mesh.get_halfedge(ind_he) {
                    (
                        ind_he,
                        (edge.first_vertex().vertex() - edge.last_vertex().vertex()).norm_squared(),
                    )
                } else {
                    (ind_he, 0.0)
                }
            })
            .collect();
        length.sort_by(|(_, l1), (_, l2)| l1.partial_cmp(l2).unwrap());
        self.non_del_edges = length.iter().map(|&(ind_he, _)| ind_he).collect();

        loop {
            if let Some(ind_he) = self.non_del_edges.pop() {
                if let Ok(edge) = self.mesh.get_halfedge(ind_he) {
                    let seg = edge.halfedge();
                    if !self.is_edge_in(&seg) {
                        let he = self.mesh.get_halfedge(ind_he)?;
                        break Ok(Some(he));
                    };
                }
            } else {
                break Ok(None);
            }
        }
    }

    /// Gets first globally non Delaunay face, starting from a shift
    pub fn get_non_del_face(&mut self) -> Result<Option<manifold_mesh3d::IterFace>> {
        loop {
            if let Some(ind_fac) = self.non_del_faces.pop() {
                if let Ok(face) = self.mesh.get_face(ind_fac) {
                    let face_vert = face.vertices_inds();
                    if !self.is_face_in(&face_vert) {
                        break Ok(Some(face));
                    };
                }
            } else {
                break Ok(None);
            }
        }
    }

    /// Flips given halfedge
    pub fn flip_halfedge(&mut self, ind_halfedge: usize) -> Result<bool> {
        mesh_operations::flip_halfedge(self.mesh, ind_halfedge)
    }

    /// Splits given halfedge
    pub fn split_halfedge(
        &mut self,
        vert: &manifold_mesh3d::Vertex,
        ind_halfedge: usize,
    ) -> Result<()> {
        let ind_near_vert = self.mesh.get_halfedge(ind_halfedge)?.first_vertex().ind();
        let (it, _) = self.vertex_edges[ind_near_vert][0];
        let ind_tet = self
            .del_struct
            .get_simplicial()
            .get_halftriangle(it)?
            .tetrahedron()
            .ind();

        log::debug!("he to split {}", ind_halfedge);
        let ind_vertex = mesh_operations::split_halfedge(self.mesh, vert, ind_halfedge)?;
        log::debug!(
            "Added {} , ({}, {}, {})",
            ind_vertex,
            vert[0],
            vert[1],
            vert[2]
        );
        self.insert_vertex(ind_vertex, ind_tet)
    }

    /// Splits given face
    pub fn split_face(&mut self, vert: &manifold_mesh3d::Vertex, ind_face: usize) -> Result<()> {
        let ind_near_vert = self.mesh.get_face(ind_face)?.vertices_inds()[0];
        let (it, _) = self.vertex_edges[ind_near_vert][0];
        let ind_tet = self
            .del_struct
            .get_simplicial()
            .get_halftriangle(it)?
            .tetrahedron()
            .ind();
        let ind_vertex = mesh_operations::split_face(self.mesh, vert, ind_face)?;
        self.insert_vertex(ind_vertex, ind_tet)
    }

    /// Check if tetrahedra are in or out
    pub fn compute_tetras_in_out(&mut self) -> Result<HashMap<[usize; 4], bool>> {
        // create Option<bool> vector for each tetrahedra
        let mut is_tetra_in = vec![None; self.del_struct.get_simplicial().get_nb_tetrahedra()];
        let mut to_check: Vec<usize> = vec![];

        // for all mesh faces, classify neighbor tetrahedra
        for ind_face in 0..self.mesh.get_nb_faces() {
            let face = self.mesh.get_face(ind_face)?;
            let [node1, node2, node3] = face.vertices_inds();
            let half_tri = self
                .del_struct
                .get_simplicial()
                .get_halftriangle_containing(
                    &Node::Value(node1),
                    &Node::Value(node2),
                    &Node::Value(node3),
                )
                .unwrap();
            let ind_tetra_in = half_tri.tetrahedron().ind();
            let ind_tetra_out = half_tri.opposite().tetrahedron().ind();

            // checks if tetrahedron are already in or not
            if let Some(val) = is_tetra_in[ind_tetra_in] {
                if val == false {
                    return Err(anyhow::Error::msg("Tetrahedron is both in and out"));
                }
            }
            if let Some(val) = is_tetra_in[ind_tetra_out] {
                if val == true {
                    return Err(anyhow::Error::msg("Tetrahedron is both in and out"));
                }
            }

            is_tetra_in[ind_tetra_in] = Some(true);
            is_tetra_in[ind_tetra_out] = Some(false);

            to_check.push(ind_tetra_in);
            to_check.push(ind_tetra_out);
        }

        // propagate insideness/ousideness for each neighbor
        while let Some(ind_tet) = to_check.pop() {
            let tet = self.del_struct.get_simplicial().get_tetrahedron(ind_tet)?;

            let val_cur = is_tetra_in[ind_tet].unwrap();

            for halftri in tet.halftriangles().iter() {
                // checks if the triangle is a mesh face
                if let [Node::Value(ind_vertex1), Node::Value(ind_vertex2), Node::Value(ind_vertex3)] =
                    halftri.nodes()
                {
                    if self
                        .mesh
                        .is_face_in(ind_vertex1, ind_vertex2, ind_vertex3)
                        .is_some()
                    {
                        continue;
                    }
                }
                // if not, sets neighbor tetrahedron as in or out
                let ind_neighbor = halftri.opposite().tetrahedron().ind();

                if let Some(val) = is_tetra_in[ind_neighbor] {
                    if val_cur != val {
                        return Err(anyhow::Error::msg("Tetrahedron went outside"));
                    }
                } else {
                    is_tetra_in[ind_neighbor] = Some(val_cur);
                    to_check.push(ind_neighbor);
                }
            }
        }

        let mut tetra_in = HashMap::new();
        for ind_tet in 0..self.del_struct.get_simplicial().get_nb_tetrahedra() {
            let tetra = self.del_struct.get_simplicial().get_tetrahedron(ind_tet)?;
            if let [Node::Value(i1), Node::Value(i2), Node::Value(i3), Node::Value(i4)] =
                tetra.nodes()
            {
                let mut tetra_ind = [i1, i2, i3, i4];
                tetra_ind.sort();
                tetra_in.insert(tetra_ind, is_tetra_in[ind_tet].unwrap());
            }
        }

        Ok(tetra_in)
    }
}
