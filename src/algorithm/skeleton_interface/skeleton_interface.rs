use anyhow::Result;
use std::collections::HashMap;

use crate::algorithm::delaunay_interface::DelaunayInterface;
use crate::mesh3d::Mesh3D;
use crate::skeleton3d::Skeleton3D;

pub struct SkeletonInterface3D<'a, 'b> {
    pub(super) mesh: &'a mut Mesh3D,
    pub(super) skeleton: &'b mut Skeleton3D,

    // existing delaunay
    pub(super) faces: HashMap<[usize; 3], Vec<[usize; 4]>>,

    // delaunay related
    pub(super) del_tet: HashMap<[usize; 4], usize>, // list of delaunay tetrahedra
    pub(super) del_tri: HashMap<[usize; 3], usize>, // list of delaunay triangles
    pub(super) del_seg: HashMap<[usize; 2], usize>, // list of delaunay segments

    // node related
    pub(super) node_tet: Vec<[usize; 4]>, // link to tetrahedron
    pub(super) node_pnode: Vec<[usize; 4]>, // partial nodes associated to each node, ordered with corners
    pub(super) node_edge: Vec<[usize; 4]>, // edges associated to each node, ordered with opposite corners

    // edge related
    pub(super) edge_tri: Vec<[usize; 3]>, // link to delaunay triangles
    pub(super) edge_pedge_dir: Vec<[usize; 3]>, // direct partial edges associated to each edge, ordered with corners
    pub(super) edge_pedge_opp: Vec<[usize; 3]>, // opposite partial edges associated to each edge, ordered with corners
    pub(super) edge_node: Vec<[Option<usize>; 2]>, // links two nodes (ordered)
    pub(super) edge_alve: Vec<[usize; 3]>,      // alveolae indices

    // alveola related
    pub(super) alve_seg: Vec<[usize; 2]>, // link to delaunay segments
    pub(super) alve_palve: Vec<[usize; 2]>, // partial alveolae associated to each face, same direction then opposite orientation
    pub(super) alve_edge: Vec<Vec<usize>>,  // lists surrouding edges

    // partial node related
    pub(super) pnode_corner: Vec<usize>, // refers to associated mesh point
    pub(super) pnode_node: Vec<usize>,   // points to a node
    pub(super) pnode_pedge_next: Vec<HashMap<usize, usize>>, // starting partial edge, indexed by partial alveolae
    pub(super) pnode_pedge_prev: Vec<HashMap<usize, usize>>, // ending partial edge, indexed by partial alveolae

    // partial edge related
    pub(super) pedge_corner: Vec<usize>, // refers to associated mesh point
    pub(super) pedge_edge: Vec<usize>,   // points to an edge
    pub(super) pedge_pnode: Vec<[Option<usize>; 2]>, // links two partial nodes (ordered)
    pub(super) pedge_palve: Vec<usize>,  // partial alveola containing partial edge
    pub(super) pedge_neigh: Vec<usize>,  // on neighbor alveola
    pub(super) pedge_opp: Vec<usize>,    // on same alveola, opposite side

    // partial alveola related
    pub(super) palve_corner: Vec<usize>, // refers to associated mesh point
    pub(super) palve_alve: Vec<usize>,   // points to an alveola
    pub(super) palve_pedge: Vec<Vec<usize>>, // pedges surrouding alveola
    pub(super) palve_opp: Vec<usize>,    // opposite partial alveola
}

#[derive(Copy, Clone)]
pub struct IterNode<'a, 'b, 'c> {
    voronoi: &'c SkeletonInterface3D<'a, 'b>,
    ind_node: usize,
}

#[derive(Copy, Clone)]
pub struct IterEdge<'a, 'b, 'c> {
    voronoi: &'c SkeletonInterface3D<'a, 'b>,
    ind_edge: usize,
}

#[derive(Copy, Clone)]
pub struct IterAlveola<'a, 'b, 'c> {
    voronoi: &'c SkeletonInterface3D<'a, 'b>,
    ind_alveola: usize,
}

#[derive(Copy, Clone)]
pub struct IterPartialNode<'a, 'b, 'c> {
    voronoi: &'c SkeletonInterface3D<'a, 'b>,
    ind_pnode: usize,
}

#[derive(Copy, Clone)]
pub struct IterPartialEdge<'a, 'b, 'c> {
    voronoi: &'c SkeletonInterface3D<'a, 'b>,
    ind_pedge: usize,
}

#[derive(Copy, Clone)]
pub struct IterPartialAlveola<'a, 'b, 'c> {
    voronoi: &'c SkeletonInterface3D<'a, 'b>,
    ind_palveola: usize,
}

impl<'a, 'b, 'c> SkeletonInterface3D<'a, 'b> {
    pub fn init(
        mesh: &'a mut Mesh3D,
        skeleton: &'b mut Skeleton3D,
    ) -> Result<SkeletonInterface3D<'a, 'b>> {
        let deltet = DelaunayInterface::from_mesh(mesh)?;
        let nb_non_del_hedges = deltet.count_non_del_halfedges()?;
        let nb_non_del_faces = deltet.count_non_del_faces()?;

        let faces = deltet
            .get_faces()
            .iter()
            .map(|(&tri, tetras)| {
                let tetras_cpy = tetras.iter().map(|&tetra| tetra).collect();
                (tri, tetras_cpy)
            })
            .collect();

        if nb_non_del_hedges != 0 || nb_non_del_faces != 0 {
            Err(anyhow::Error::msg("Mesh is not Delaunay"))
        } else {
            Ok(SkeletonInterface3D {
                mesh,
                skeleton,
                faces,
                del_tet: HashMap::new(),
                del_tri: HashMap::new(),
                del_seg: HashMap::new(),
                node_tet: Vec::new(),
                node_pnode: Vec::new(),
                node_edge: Vec::new(),
                edge_tri: Vec::new(),
                edge_pedge_dir: Vec::new(),
                edge_pedge_opp: Vec::new(),
                edge_node: Vec::new(),
                edge_alve: Vec::new(),
                alve_seg: Vec::new(),
                alve_palve: Vec::new(),
                alve_edge: Vec::new(),
                pnode_corner: Vec::new(),
                pnode_node: Vec::new(),
                pnode_pedge_next: Vec::new(),
                pnode_pedge_prev: Vec::new(),
                pedge_corner: Vec::new(),
                pedge_edge: Vec::new(),
                pedge_pnode: Vec::new(),
                pedge_palve: Vec::new(),
                pedge_neigh: Vec::new(),
                pedge_opp: Vec::new(),
                palve_corner: Vec::new(),
                palve_alve: Vec::new(),
                palve_pedge: Vec::new(),
                palve_opp: Vec::new(),
            })
        }
    }

    pub fn add_node(&'c mut self, del_tet: &[usize; 4]) -> Result<IterNode<'a, 'b, 'c>> {
        if let Some(&ind_node) = self.del_tet.get(del_tet) {
            return Ok(IterNode {
                voronoi: self,
                ind_node,
            });
        }

        let ind_node = self.del_tet.len();

        self.del_tet.insert(*del_tet, ind_node);

        // node
        self.node_tet.push(*del_tet);

        // partial nodes
        let pnodes = self.add_partial_nodes(ind_node, del_tet);
        self.node_pnode.push(pnodes);

        // edges
        let ind_edg0 = self.add_edge(&[del_tet[1], del_tet[2], del_tet[3]]);
        let ind_edg1 = self.add_edge(&[del_tet[0], del_tet[2], del_tet[3]]);
        let ind_edg2 = self.add_edge(&[del_tet[0], del_tet[1], del_tet[3]]);
        let ind_edg3 = self.add_edge(&[del_tet[0], del_tet[1], del_tet[2]]);

        self.link_node_edges(ind_node, [ind_edg0, ind_edg1, ind_edg2, ind_edg3])?;

        Ok(IterNode {
            voronoi: self,
            ind_node,
        })
    }

    fn add_partial_nodes(&mut self, ind_node: usize, del_tet: &[usize; 4]) -> [usize; 4] {
        let mut pnodes = [0; 4];
        for i in 0..4 {
            let ind_pnode = self.pnode_node.len();
            self.pnode_corner.push(del_tet[i]);
            self.pnode_node.push(ind_node);
            self.pnode_pedge_next.push(HashMap::new());
            self.pnode_pedge_prev.push(HashMap::new());
            pnodes[i] = ind_pnode;
        }
        pnodes
    }

    fn add_edge(&mut self, del_tri: &[usize; 3]) -> usize {
        let ind_edge = match self.del_tri.get(del_tri) {
            Some(&ind_edge) => ind_edge,
            None => {
                let ind_edge = self.del_tri.len();
                self.del_tri.insert(*del_tri, ind_edge);

                self.edge_tri.push(*del_tri);
                self.edge_node.push([None, None]);
                let (pedges_dir, pedges_opp) = self.add_partial_edges(ind_edge, del_tri);
                self.edge_pedge_dir.push(pedges_dir);
                self.edge_pedge_opp.push(pedges_opp);

                let ind_alv0 = self.add_alveola(&[del_tri[1], del_tri[2]]);
                let ind_alv1 = self.add_alveola(&[del_tri[0], del_tri[2]]);
                let ind_alv2 = self.add_alveola(&[del_tri[0], del_tri[1]]);

                self.link_edge_alves(ind_edge, [ind_alv0, ind_alv1, ind_alv2]);

                ind_edge
            }
        };
        ind_edge
    }

    fn add_partial_edges(
        &mut self,
        ind_edge: usize,
        del_tri: &[usize; 3],
    ) -> ([usize; 3], [usize; 3]) {
        let mut pedges_dir = [0; 3];
        let mut pedges_opp = [0; 3];
        for i in 0..3 {
            let ind_pedge_dir = self.pedge_edge.len();
            self.pedge_corner.push(del_tri[i]);
            self.pedge_edge.push(ind_edge);
            self.pedge_pnode.push([None, None]);
            self.pedge_palve.push(0); // temp value
            pedges_dir[i] = ind_pedge_dir;
        }
        for i in 0..3 {
            let ind_pedge_opp = self.pedge_edge.len();
            self.pedge_corner.push(del_tri[i]);
            self.pedge_edge.push(ind_edge);
            self.pedge_pnode.push([None, None]);
            self.pedge_palve.push(0); // temp value
            pedges_opp[i] = ind_pedge_opp;
        }
        // neigh: same corner
        self.pedge_neigh.push(pedges_opp[0]);
        self.pedge_neigh.push(pedges_opp[1]);
        self.pedge_neigh.push(pedges_opp[2]);
        self.pedge_neigh.push(pedges_dir[0]);
        self.pedge_neigh.push(pedges_dir[1]);
        self.pedge_neigh.push(pedges_dir[2]);
        // opp: same alveola
        self.pedge_opp.push(pedges_opp[1]);
        self.pedge_opp.push(pedges_opp[2]);
        self.pedge_opp.push(pedges_opp[0]);
        self.pedge_opp.push(pedges_dir[2]);
        self.pedge_opp.push(pedges_dir[0]);
        self.pedge_opp.push(pedges_dir[1]);
        (pedges_dir, pedges_opp)
    }

    fn add_alveola(&mut self, del_seg: &[usize; 2]) -> usize {
        let ind_alve = match self.del_seg.get(del_seg) {
            Some(&ind_alve) => ind_alve,
            None => {
                let ind_alve = self.del_seg.len();
                self.del_seg.insert(*del_seg, ind_alve);
                self.alve_seg.push(*del_seg);
                self.alve_edge.push(Vec::new());
                self.add_partial_alveolae(ind_alve, del_seg);
                ind_alve
            }
        };
        ind_alve
    }

    fn add_partial_alveolae(&mut self, ind_alve: usize, del_seg: &[usize; 2]) -> () {
        let mut array_palveolae = [0; 2];
        for i in 0..2 {
            let ind_palve = self.palve_alve.len();
            self.palve_alve.push(ind_alve);
            self.palve_corner.push(del_seg[i]);
            self.palve_pedge.push(Vec::new());
            array_palveolae[i] = ind_palve;
        }
        self.palve_opp.push(self.palve_alve.len() - 1);
        self.palve_opp.push(self.palve_alve.len() - 2);
        self.alve_palve.push(array_palveolae);
    }

    fn link_node_edges(&mut self, ind_node: usize, ind_edges: [usize; 4]) -> Result<()> {
        self.node_edge.push(ind_edges);
        let side = if self.edge_node[ind_edges[0]][0] == None
            && self.edge_node[ind_edges[1]][1] == None
            && self.edge_node[ind_edges[2]][0] == None
            && self.edge_node[ind_edges[3]][1] == None
        {
            Ok([0, 1, 0, 1])
        } else if self.edge_node[ind_edges[0]][1] == None
            && self.edge_node[ind_edges[1]][0] == None
            && self.edge_node[ind_edges[2]][1] == None
            && self.edge_node[ind_edges[3]][0] == None
        {
            Ok([1, 0, 1, 0])
        } else {
            Err(anyhow::Error::msg("Contradiction in skeletons nodes"))
        }?;

        for i in 0..4 {
            let ind_edg = ind_edges[i];
            self.edge_node[ind_edg][side[i]] = Some(ind_node);

            let (ind_pedges_beg, ind_pedges_end) = if side[i] == 0 {
                (self.edge_pedge_dir[ind_edg], self.edge_pedge_opp[ind_edg])
            } else {
                (self.edge_pedge_opp[ind_edg], self.edge_pedge_dir[ind_edg])
            };

            for j in 0..3 {
                let assoc_nod = if j < i { j } else { j + 1 };
                let ind_pedge_beg = ind_pedges_beg[j];
                let ind_pedge_end = ind_pedges_end[j];
                let ind_pnode = self.node_pnode[ind_node][assoc_nod];

                let ind_palve_beg = self.pedge_palve[ind_pedge_beg];
                self.pnode_pedge_next[ind_pnode].insert(ind_palve_beg, ind_pedge_beg);

                let ind_palve_end = self.pedge_palve[ind_pedge_end];
                self.pnode_pedge_prev[ind_pnode].insert(ind_palve_end, ind_pedge_end);

                self.pedge_pnode[ind_pedge_beg][0] = Some(ind_pnode);
                self.pedge_pnode[ind_pedge_end][1] = Some(ind_pnode);
            }
        }
        Ok(())
    }

    fn link_edge_alves(&mut self, ind_edge: usize, ind_alve: [usize; 3]) -> () {
        self.edge_alve.push(ind_alve);

        for alv in ind_alve {
            self.alve_edge[alv].push(ind_edge);
        }

        self.pedge_palve[self.edge_pedge_dir[ind_edge][0]] = self.alve_palve[ind_alve[2]][0];
        self.pedge_palve[self.edge_pedge_dir[ind_edge][1]] = self.alve_palve[ind_alve[0]][0];
        self.pedge_palve[self.edge_pedge_dir[ind_edge][2]] = self.alve_palve[ind_alve[1]][1];
        self.pedge_palve[self.edge_pedge_opp[ind_edge][0]] = self.alve_palve[ind_alve[1]][0];
        self.pedge_palve[self.edge_pedge_opp[ind_edge][1]] = self.alve_palve[ind_alve[2]][1];
        self.pedge_palve[self.edge_pedge_opp[ind_edge][2]] = self.alve_palve[ind_alve[0]][1];

        self.palve_pedge[self.alve_palve[ind_alve[2]][0]].push(self.edge_pedge_dir[ind_edge][0]);
        self.palve_pedge[self.alve_palve[ind_alve[0]][0]].push(self.edge_pedge_dir[ind_edge][1]);
        self.palve_pedge[self.alve_palve[ind_alve[1]][1]].push(self.edge_pedge_dir[ind_edge][2]);
        self.palve_pedge[self.alve_palve[ind_alve[1]][0]].push(self.edge_pedge_opp[ind_edge][0]);
        self.palve_pedge[self.alve_palve[ind_alve[2]][1]].push(self.edge_pedge_opp[ind_edge][1]);
        self.palve_pedge[self.alve_palve[ind_alve[0]][1]].push(self.edge_pedge_opp[ind_edge][2]);
    }

    pub fn get_node(&'c self, ind_node: usize) -> Result<IterNode<'a, 'b, 'c>> {
        if ind_node >= self.node_tet.len() {
            return Err(anyhow::Error::msg("Node index out of bounds"));
        }
        Ok(IterNode {
            voronoi: self,
            ind_node,
        })
    }

    pub fn get_partial_node(&'c self, ind_pnode: usize) -> Result<IterPartialNode<'a, 'b, 'c>> {
        if ind_pnode >= self.pnode_node.len() {
            return Err(anyhow::Error::msg("Partial node index out of bounds"));
        }
        Ok(IterPartialNode {
            voronoi: self,
            ind_pnode,
        })
    }

    pub fn get_edge(&'c self, ind_edge: usize) -> Result<IterEdge<'a, 'b, 'c>> {
        if ind_edge >= self.edge_tri.len() {
            return Err(anyhow::Error::msg("Edge index out of bounds"));
        }
        Ok(IterEdge {
            voronoi: self,
            ind_edge,
        })
    }

    pub fn get_partial_edge(&'c self, ind_pedge: usize) -> Result<IterPartialEdge<'a, 'b, 'c>> {
        if ind_pedge >= self.pedge_edge.len() {
            return Err(anyhow::Error::msg("Partial edge index out of bounds"));
        }
        Ok(IterPartialEdge {
            voronoi: self,
            ind_pedge,
        })
    }

    pub fn get_alveola(&'c self, ind_alveola: usize) -> Result<IterAlveola<'a, 'b, 'c>> {
        if ind_alveola >= self.alve_seg.len() {
            return Err(anyhow::Error::msg("Alveola index out of bounds"));
        }
        Ok(IterAlveola {
            voronoi: self,
            ind_alveola,
        })
    }

    pub fn get_partial_alveola(
        &'c self,
        ind_palveola: usize,
    ) -> Result<IterPartialAlveola<'a, 'b, 'c>> {
        if ind_palveola >= self.palve_alve.len() {
            return Err(anyhow::Error::msg("Partial alveola index out of bounds"));
        }
        Ok(IterPartialAlveola {
            voronoi: self,
            ind_palveola,
        })
    }

    pub fn get_skeleton(&self) -> &Skeleton3D {
        self.skeleton
    }

    pub fn get_mesh(&self) -> &Mesh3D {
        self.mesh
    }

    pub fn get_tetrahedra_from_triangle(&self, del_tri: [usize; 3]) -> Result<Vec<[usize; 4]>> {
        let vec = self
            .faces
            .get(&del_tri)
            .ok_or(anyhow::Error::msg("Triangle does not exist"))?
            .iter()
            .map(|&x| x)
            .collect();
        Ok(vec)
    }

    fn check_node(&self, ind_node: usize) -> Result<()> {
        let node = self.get_node(ind_node)?;

        // checking link between edges and node
        for edge in node.edges() {
            let in_edge = edge
                .nodes()
                .iter()
                .fold(false, |b, nod| b || nod.ind() == node.ind());
            if !in_edge {
                let msg = format!("Node {} not contained in edge {}", node.ind(), edge.ind());
                return Err(anyhow::Error::msg(msg));
            }
        }

        for pnode in node.partial_nodes() {
            if pnode.node().ind() != node.ind() {
                let msg = format!(
                    "Node {} not contained in partial node {}",
                    node.ind(),
                    pnode.ind()
                );
                return Err(anyhow::Error::msg(msg));
            }

            for pedge in pnode.partial_edge_next() {
                match pedge.partial_node_first() {
                    None => {
                        let msg = format!("No first node for partial edge {}", pedge.ind());
                        return Err(anyhow::Error::msg(msg));
                    }
                    Some(pn) => {
                        if pn.ind() != pnode.ind() {
                            let msg = format!(
                                "Partial node {} not first node of partial edge {}",
                                pnode.ind(),
                                pedge.ind()
                            );
                            return Err(anyhow::Error::msg(msg));
                        }
                    }
                }
            }

            for pedge in pnode.partial_edge_prev() {
                match pedge.partial_node_last() {
                    None => {
                        let msg = format!("No last node for partial edge {}", pedge.ind());
                        return Err(anyhow::Error::msg(msg));
                    }
                    Some(pn) => {
                        if pn.ind() != pnode.ind() {
                            let msg = format!(
                                "Partial node {} not last node of partial edge {}",
                                pnode.ind(),
                                pedge.ind()
                            );
                            return Err(anyhow::Error::msg(msg));
                        }
                    }
                }
            }
        }

        Ok(())
    }

    fn check_edge(&self, ind_edge: usize) -> Result<()> {
        let edge = self.get_edge(ind_edge)?;

        for pedge in edge.partial_edges() {
            if pedge.edge().ind() != edge.ind() {
                let msg = format!(
                    "Partial edge {} not connected to edge {}",
                    pedge.ind(),
                    edge.ind()
                );
                return Err(anyhow::Error::msg(msg));
            }
        }

        for alve in edge.alveolae() {
            let is_in = alve
                .edges()
                .iter()
                .fold(false, |b, edg| b || edg.ind() == edge.ind());
            if !is_in {
                let msg = format!("Edge {} not in alveola {}", edge.ind(), alve.ind());
                return Err(anyhow::Error::msg(msg));
            }
        }
        Ok(())
    }

    fn check_partial_edge(&self, ind_pedge: usize) -> Result<()> {
        let pedge = self.get_partial_edge(ind_pedge)?;

        let edge = pedge.edge();
        let is_in_edg = edge
            .partial_edges()
            .iter()
            .fold(false, |b, ped| b || ped.ind() == pedge.ind());
        if !is_in_edg {
            let msg = format!("Partial edge {} not in edge {}", pedge.ind(), edge.ind());
            return Err(anyhow::Error::msg(msg));
        }

        if let Some(pnode) = pedge.partial_node_first() {
            let is_in_pnd = pnode
                .partial_edge_next()
                .iter()
                .fold(false, |b, ped| b || ped.ind() == pedge.ind());
            if !is_in_pnd {
                let msg = format!(
                    "Partial edge {} not a next of partial node {}",
                    pedge.ind(),
                    pnode.ind()
                );
                return Err(anyhow::Error::msg(msg));
            }
        }

        if let Some(pnode) = pedge.partial_node_last() {
            let is_in_pnd = pnode
                .partial_edge_prev()
                .iter()
                .fold(false, |b, ped| b || ped.ind() == pedge.ind());
            if !is_in_pnd {
                let msg = format!(
                    "Partial edge {} not a prev of partial node {}",
                    pedge.ind(),
                    pnode.ind()
                );
                return Err(anyhow::Error::msg(msg));
            }
        }

        if let Some(pedge_next) = pedge.partial_edge_next()? {
            if pedge_next.partial_alveolae().ind() != pedge.partial_alveolae().ind() {
                let msg = format!(
                    "Partial edge {} and next partial edge {} not on same alveola",
                    pedge.ind(),
                    pedge_next.ind()
                );
                return Err(anyhow::Error::msg(msg));
            }

            let nod_last = pedge.partial_node_last().ok_or(anyhow::Error::msg(
                "Last node of current partial edge should exist",
            ))?;

            let nod_first = pedge_next.partial_node_first().ok_or(anyhow::Error::msg(
                "First node of next partial edge should exist",
            ))?;

            if nod_first.ind() != nod_last.ind() {
                let msg = format!(
                    "Last node of partial edge {} different to first node of partial edge {}",
                    pedge.ind(),
                    pedge_next.ind()
                );
                return Err(anyhow::Error::msg(msg));
            }
        }

        if let Some(pedge_prev) = pedge.partial_edge_prev()? {
            if pedge_prev.partial_alveolae().ind() != pedge.partial_alveolae().ind() {
                let msg = format!(
                    "Partial edge {} and prev partial edge {} not on same alveola",
                    pedge.ind(),
                    pedge_prev.ind()
                );
                return Err(anyhow::Error::msg(msg));
            }

            let nod_last = pedge_prev.partial_node_last().ok_or(anyhow::Error::msg(
                "Last node of previous partial edge should exist",
            ))?;

            let nod_first = pedge.partial_node_first().ok_or(anyhow::Error::msg(
                "First node of current partial edge should exist",
            ))?;

            if nod_first.ind() != nod_last.ind() {
                let msg = format!(
                    "First node of partial edge {} different to last node of partial edge {}",
                    pedge.ind(),
                    pedge_prev.ind()
                );
                return Err(anyhow::Error::msg(msg));
            }
        }

        Ok(())
    }

    fn check_alveola(&self, ind_alve: usize) -> Result<()> {
        let alve = self.get_alveola(ind_alve)?;

        for edge in alve.edges() {
            let is_in = edge
                .alveolae()
                .iter()
                .fold(false, |b, alv| b || alv.ind() == alve.ind());
            if !is_in {
                let msg = format!(
                    "Edge {} not connected to alveola {}",
                    edge.ind(),
                    alve.ind()
                );
                return Err(anyhow::Error::msg(msg));
            }
        }

        for palve in alve.partial_alveolae() {
            if palve.alveola().ind() != alve.ind() {
                let msg = format!(
                    "Partial alveola {} not connected to alveola {}",
                    palve.ind(),
                    alve.ind()
                );
                return Err(anyhow::Error::msg(msg));
            }
        }

        Ok(())
    }

    pub fn check(&self) -> Result<()> {
        for ind_node in 0..self.node_tet.len() {
            self.check_node(ind_node)?;
        }
        for ind_edge in 0..self.edge_tri.len() {
            self.check_edge(ind_edge)?;
        }
        for ind_alve in 0..self.alve_seg.len() {
            self.check_alveola(ind_alve)?;
        }
        for ind_pedge in 0..self.pedge_edge.len() {
            self.check_partial_edge(ind_pedge)?;
        }
        Ok(())
    }

    pub fn print_node(&self, ind_node: usize) -> () {
        let tet = self.node_tet[ind_node];
        let pnods = self.node_pnode[ind_node];
        let edges = self.node_edge[ind_node];
        println!("{} node", ind_node);
        println!(
            "  Tetrahedron: ({}, {}, {}, {})",
            tet[0], tet[1], tet[2], tet[3]
        );
        println!(
            "  Partial nodes: {}, {}, {}, {}",
            pnods[0], pnods[1], pnods[2], pnods[3]
        );
        println!(
            "  Edges: {}, {}, {}, {}",
            edges[0], edges[1], edges[2], edges[3]
        );
    }

    pub fn print_edge(&self, ind_edge: usize) -> () {
        let tri = self.edge_tri[ind_edge];
        let nods = self.edge_node[ind_edge];
        let pedg_dir = self.edge_pedge_dir[ind_edge];
        let pedg_opp = self.edge_pedge_opp[ind_edge];
        let alv = self.edge_alve[ind_edge];
        println!("{} edge", ind_edge);
        println!("  Triangle: ({}, {}, {})", tri[0], tri[1], tri[2]);
        print!("  Nodes, (");
        if let Some(nod) = nods[0] {
            print!("{}, ", nod);
        } else {
            print!("[], ");
        }
        if let Some(nod) = nods[1] {
            println!("{})", nod);
        } else {
            println!("[])");
        }
        println!("  Alveolae: {}, {}, {}", alv[0], alv[1], alv[2]);
        println!(
            "  Partial edges dir: ({}, {}, {})",
            pedg_dir[0], pedg_dir[1], pedg_dir[2]
        );
        println!(
            "  Partial edges opp: ({}, {}, {})",
            pedg_opp[0], pedg_opp[1], pedg_opp[2]
        );
    }

    pub fn print_alveola(&self, ind_alve: usize) -> () {
        let seg = self.alve_seg[ind_alve];
        let edgs = self.alve_edge[ind_alve].iter();
        let palv = self.alve_palve[ind_alve];
        println!("{} alveola", ind_alve);
        println!("  Segment: ({}, {})", seg[0], seg[1]);
        print!("  Edges: ");
        for edg in edgs {
            print!("{}, ", edg);
        }
        println!("");
        println!("  Partial alveolae: ({}, {})", palv[0], palv[1]);
    }

    pub fn print_partial_node(&self, ind_pnode: usize) -> () {
        let nod = self.pnode_node[ind_pnode];
        let corner = self.pnode_corner[ind_pnode];
        let pedg_next = self.pnode_pedge_next[ind_pnode].iter();
        let pedg_prev = self.pnode_pedge_prev[ind_pnode].iter();
        println!("{} partial node; node {}", ind_pnode, nod);
        println!("  Corner: {}", corner);
        print!("  Next: ");
        for (a, e) in pedg_next {
            print!(" ({}, {}), ", a, e);
        }
        println!("");
        print!("  Prev: ");
        for (a, e) in pedg_prev {
            print!(" ({}, {}), ", a, e);
        }
        println!("");
    }

    pub fn print_partial_edge(&self, ind_pedge: usize) -> () {
        let edg = self.pedge_edge[ind_pedge];
        let corner = self.pedge_corner[ind_pedge];
        let pnods = self.pedge_pnode[ind_pedge];
        let palve = self.pedge_palve[ind_pedge];
        let pneigh = self.pedge_neigh[ind_pedge];
        let popp = self.pedge_opp[ind_pedge];
        println!("{} partial edge; edge {}", ind_pedge, edg);
        println!("  Corner: {}", corner);
        print!("  Partial nodes, (");
        if let Some(pnod) = pnods[0] {
            print!("{}, ", pnod);
        } else {
            print!("[], ");
        }
        if let Some(pnod) = pnods[1] {
            println!("{})", pnod);
        } else {
            println!("[])");
        }
        println!("  Partial alveola: {}", palve);
        println!("  Neighbor partial edge: {}", pneigh);
        println!("  Opposite partial edge: {}", popp);
    }

    pub fn print_all(&self) -> () {
        println!("Nodes");
        for ind_node in 0..self.node_tet.len() {
            self.print_node(ind_node);
        }
        println!("");
        println!("Edges");
        for ind_edge in 0..self.edge_tri.len() {
            self.print_edge(ind_edge);
        }
        println!("");
        println!("Alveolae");
        for ind_alve in 0..self.alve_seg.len() {
            self.print_alveola(ind_alve);
        }
        println!("");
        println!("Partial nodes");
        for ind_pnode in 0..self.pnode_node.len() {
            self.print_partial_node(ind_pnode);
        }
        println!("");
        println!("Partial edges");
        for ind_pedge in 0..self.pedge_edge.len() {
            self.print_partial_edge(ind_pedge);
        }
        println!("");
    }

    pub fn fully_computed(&self) -> Result<bool> {
        for i in 0..self.alve_seg.len() {
            let alve = self.get_alveola(i)?;

            if !alve.is_computed() && alve.is_in()? {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

impl<'a, 'b, 'c> IterNode<'a, 'b, 'c> {
    pub fn ind(&self) -> usize {
        self.ind_node
    }

    pub fn delaunay_tetrahedron(&self) -> [usize; 4] {
        self.voronoi.node_tet[self.ind_node]
    }

    pub fn partial_nodes(&self) -> [IterPartialNode<'a, 'b, 'c>; 4] {
        [
            IterPartialNode {
                voronoi: self.voronoi,
                ind_pnode: self.voronoi.node_pnode[self.ind_node][0],
            },
            IterPartialNode {
                voronoi: self.voronoi,
                ind_pnode: self.voronoi.node_pnode[self.ind_node][1],
            },
            IterPartialNode {
                voronoi: self.voronoi,
                ind_pnode: self.voronoi.node_pnode[self.ind_node][2],
            },
            IterPartialNode {
                voronoi: self.voronoi,
                ind_pnode: self.voronoi.node_pnode[self.ind_node][3],
            },
        ]
    }

    pub fn edges(&self) -> [IterEdge<'a, 'b, 'c>; 4] {
        [
            IterEdge {
                voronoi: self.voronoi,
                ind_edge: self.voronoi.node_edge[self.ind_node][0],
            },
            IterEdge {
                voronoi: self.voronoi,
                ind_edge: self.voronoi.node_edge[self.ind_node][1],
            },
            IterEdge {
                voronoi: self.voronoi,
                ind_edge: self.voronoi.node_edge[self.ind_node][2],
            },
            IterEdge {
                voronoi: self.voronoi,
                ind_edge: self.voronoi.node_edge[self.ind_node][3],
            },
        ]
    }
}

impl<'a, 'b, 'c> IterEdge<'a, 'b, 'c> {
    pub fn ind(&self) -> usize {
        self.ind_edge
    }

    pub fn delaunay_triangle(&self) -> [usize; 3] {
        self.voronoi.edge_tri[self.ind_edge]
    }

    pub fn nodes(&self) -> Vec<IterNode<'a, 'b, 'c>> {
        let mut nods: Vec<IterNode> = Vec::new();
        self.voronoi.edge_node[self.ind_edge]
            .iter()
            .fold(&mut nods, |v, &nod| {
                if let Some(ind_node) = nod {
                    v.push(IterNode {
                        voronoi: self.voronoi,
                        ind_node,
                    })
                };
                v
            });
        nods
    }

    pub fn is_computed(&self) -> bool {
        let nods = self.voronoi.edge_node[self.ind_edge];
        nods[0].is_some() && nods[1].is_some()
    }

    pub fn is_in(&self) -> Result<bool> {
        let tri = self.delaunay_triangle();
        Ok(self
            .voronoi
            .mesh
            .is_face_in(tri[0], tri[1], tri[2])?
            .is_none())
    }

    pub fn alveolae(&self) -> [IterAlveola<'a, 'b, 'c>; 3] {
        [
            IterAlveola {
                voronoi: self.voronoi,
                ind_alveola: self.voronoi.edge_alve[self.ind_edge][0],
            },
            IterAlveola {
                voronoi: self.voronoi,
                ind_alveola: self.voronoi.edge_alve[self.ind_edge][1],
            },
            IterAlveola {
                voronoi: self.voronoi,
                ind_alveola: self.voronoi.edge_alve[self.ind_edge][2],
            },
        ]
    }

    pub fn partial_edges(&self) -> [IterPartialEdge<'a, 'b, 'c>; 6] {
        [
            IterPartialEdge {
                voronoi: self.voronoi,
                ind_pedge: self.voronoi.edge_pedge_dir[self.ind_edge][0],
            },
            IterPartialEdge {
                voronoi: self.voronoi,
                ind_pedge: self.voronoi.edge_pedge_dir[self.ind_edge][1],
            },
            IterPartialEdge {
                voronoi: self.voronoi,
                ind_pedge: self.voronoi.edge_pedge_dir[self.ind_edge][2],
            },
            IterPartialEdge {
                voronoi: self.voronoi,
                ind_pedge: self.voronoi.edge_pedge_opp[self.ind_edge][0],
            },
            IterPartialEdge {
                voronoi: self.voronoi,
                ind_pedge: self.voronoi.edge_pedge_opp[self.ind_edge][1],
            },
            IterPartialEdge {
                voronoi: self.voronoi,
                ind_pedge: self.voronoi.edge_pedge_opp[self.ind_edge][2],
            },
        ]
    }
}

impl<'a, 'b, 'c> IterAlveola<'a, 'b, 'c> {
    pub fn ind(&self) -> usize {
        self.ind_alveola
    }

    pub fn delaunay_segment(&self) -> [usize; 2] {
        self.voronoi.alve_seg[self.ind_alveola]
    }

    pub fn edges(&self) -> Vec<IterEdge<'a, 'b, 'c>> {
        self.voronoi.alve_edge[self.ind_alveola]
            .iter()
            .map(|&ind_edge| IterEdge {
                voronoi: self.voronoi,
                ind_edge,
            })
            .collect()
    }

    pub fn is_computed(&self) -> bool {
        self.voronoi.alve_edge[self.ind_alveola]
            .iter()
            .map(|&ind_edge| IterEdge {
                voronoi: self.voronoi,
                ind_edge,
            })
            .fold(true, |b, ed| b && ed.is_computed())
    }

    pub fn is_in(&self) -> Result<bool> {
        let seg = self.delaunay_segment();
        Ok(self.voronoi.mesh.is_edge_in(seg[0], seg[1])?.is_none())
    }

    pub fn partial_alveolae(&self) -> [IterPartialAlveola<'a, 'b, 'c>; 2] {
        [
            IterPartialAlveola {
                voronoi: self.voronoi,
                ind_palveola: self.voronoi.alve_palve[self.ind_alveola][0],
            },
            IterPartialAlveola {
                voronoi: self.voronoi,
                ind_palveola: self.voronoi.alve_palve[self.ind_alveola][1],
            },
        ]
    }
}

impl<'a, 'b, 'c> IterPartialNode<'a, 'b, 'c> {
    pub fn ind(&self) -> usize {
        self.ind_pnode
    }

    pub fn corner(&self) -> usize {
        self.voronoi.pnode_corner[self.ind_pnode]
    }

    pub fn node(&self) -> IterNode<'a, 'b, 'c> {
        IterNode {
            voronoi: self.voronoi,
            ind_node: self.voronoi.pnode_node[self.ind_pnode],
        }
    }

    pub fn partial_edge_prev(&self) -> Vec<IterPartialEdge<'a, 'b, 'c>> {
        self.voronoi.pnode_pedge_prev[self.ind_pnode]
            .iter()
            .map(|(_ind_palve, &ind_pedge)| IterPartialEdge {
                voronoi: self.voronoi,
                ind_pedge,
            })
            .collect()
    }

    fn partial_edge_prev_on_alve(&self, ind_palve: usize) -> Result<IterPartialEdge<'a, 'b, 'c>> {
        let &ind_pedge = self.voronoi.pnode_pedge_prev[self.ind_pnode]
            .get(&ind_palve)
            .ok_or(anyhow::Error::msg(
                "Could not find previous partial edge from partial alveola",
            ))?;
        Ok(IterPartialEdge {
            voronoi: self.voronoi,
            ind_pedge,
        })
    }

    pub fn partial_edge_next(&self) -> Vec<IterPartialEdge<'a, 'b, 'c>> {
        self.voronoi.pnode_pedge_next[self.ind_pnode]
            .iter()
            .map(|(_ind_palve, &ind_pedge)| IterPartialEdge {
                voronoi: self.voronoi,
                ind_pedge,
            })
            .collect()
    }

    fn partial_edge_next_on_alve(&self, ind_palve: usize) -> Result<IterPartialEdge<'a, 'b, 'c>> {
        let &ind_pedge = self.voronoi.pnode_pedge_next[self.ind_pnode]
            .get(&ind_palve)
            .ok_or(anyhow::Error::msg(
                "Could not find next partial edge from partial alveola",
            ))?;
        Ok(IterPartialEdge {
            voronoi: self.voronoi,
            ind_pedge,
        })
    }
}

impl<'a, 'b, 'c> IterPartialEdge<'a, 'b, 'c> {
    pub fn ind(&self) -> usize {
        self.ind_pedge
    }

    pub fn corner(&self) -> usize {
        self.voronoi.pedge_corner[self.ind_pedge]
    }

    pub fn edge(&self) -> IterEdge<'a, 'b, 'c> {
        IterEdge {
            voronoi: self.voronoi,
            ind_edge: self.voronoi.pedge_edge[self.ind_pedge],
        }
    }

    pub fn partial_node_first(&self) -> Option<IterPartialNode<'a, 'b, 'c>> {
        match self.voronoi.pedge_pnode[self.ind_pedge][0] {
            None => None,
            Some(ind_pnode) => Some(IterPartialNode {
                voronoi: self.voronoi,
                ind_pnode,
            }),
        }
    }

    pub fn partial_node_last(&self) -> Option<IterPartialNode<'a, 'b, 'c>> {
        match self.voronoi.pedge_pnode[self.ind_pedge][1] {
            None => None,
            Some(ind_pnode) => Some(IterPartialNode {
                voronoi: self.voronoi,
                ind_pnode,
            }),
        }
    }

    pub fn partial_edge_opposite(&self) -> IterPartialEdge<'a, 'b, 'c> {
        IterPartialEdge {
            voronoi: self.voronoi,
            ind_pedge: self.voronoi.pedge_opp[self.ind_pedge],
        }
    }

    pub fn partial_edge_neighbor(&self) -> IterPartialEdge<'a, 'b, 'c> {
        IterPartialEdge {
            voronoi: self.voronoi,
            ind_pedge: self.voronoi.pedge_neigh[self.ind_pedge],
        }
    }

    pub fn partial_edge_prev(&self) -> Result<Option<IterPartialEdge<'a, 'b, 'c>>> {
        match self.partial_node_first() {
            None => Ok(None),
            Some(pnode) => {
                let palve = self.partial_alveolae().ind();
                let pedge = pnode.partial_edge_prev_on_alve(palve)?;
                Ok(Some(pedge))
            }
        }
    }

    pub fn partial_edge_next(&self) -> Result<Option<IterPartialEdge<'a, 'b, 'c>>> {
        match self.partial_node_last() {
            None => Ok(None),
            Some(pnode) => {
                let palve = self.partial_alveolae().ind();
                let pedge = pnode.partial_edge_next_on_alve(palve)?;
                Ok(Some(pedge))
            }
        }
    }

    pub fn partial_alveolae(&self) -> IterPartialAlveola<'a, 'b, 'c> {
        IterPartialAlveola {
            voronoi: self.voronoi,
            ind_palveola: self.voronoi.pedge_palve[self.ind_pedge],
        }
    }
}

impl<'a, 'b, 'c> IterPartialAlveola<'a, 'b, 'c> {
    pub fn ind(&self) -> usize {
        self.ind_palveola
    }

    pub fn corner(&self) -> usize {
        self.voronoi.palve_corner[self.ind_palveola]
    }

    pub fn alveola(&self) -> IterAlveola<'a, 'b, 'c> {
        IterAlveola {
            voronoi: self.voronoi,
            ind_alveola: self.voronoi.palve_alve[self.ind_palveola],
        }
    }

    pub fn partial_alveola_opposite(&self) -> IterPartialAlveola<'a, 'b, 'c> {
        IterPartialAlveola {
            voronoi: self.voronoi,
            ind_palveola: self.voronoi.palve_opp[self.ind_palveola],
        }
    }

    pub fn partial_edges(&self) -> Vec<IterPartialEdge<'a, 'b, 'c>> {
        self.voronoi.palve_pedge[self.ind_palveola]
            .iter()
            .map(|&ind_pedge| IterPartialEdge {
                voronoi: self.voronoi,
                ind_pedge,
            })
            .collect()
    }
}
