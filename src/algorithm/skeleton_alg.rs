use crate::algorithm::skeleton_interface::{skeleton_operations, SkeletonInterface3D};
use crate::mesh3d::Mesh3D;
use crate::skeleton3d::Skeleton3D;
use anyhow::Result;

pub fn full_skeletonization(mesh: &mut Mesh3D, skeleton: &mut Skeleton3D) -> Result<()> {
    println!("Init skeleton interface");
    let mut skeleton_interface = SkeletonInterface3D::init(mesh, skeleton)?;

    println!("Finding some first node");
    let ind_first_node = skeleton_operations::first_node_in(&mut skeleton_interface)?;

    println!("Finding some first alveolae");
    let cur_node = skeleton_interface.get_node(ind_first_node)?;

    let edges = cur_node.edges();

    let mut vec_alveola = Vec::new();
    for edge in edges {
        let tri = edge.delaunay_triangle();

        if skeleton_interface
            .get_mesh()
            .is_face_in(tri[0], tri[1], tri[2])?
            .is_none()
        {
            for alve in edge.alveolae() {
                let seg = alve.delaunay_segment();
                if skeleton_interface
                    .get_mesh()
                    .is_edge_in(seg[0], seg[1])?
                    .is_none()
                {
                    vec_alveola.push(alve.ind());
                }
            }
        }
    }

    println!("Propagating skeleton");
    loop {
        if let Some(ind_alveola) = vec_alveola.pop() {
            let alveola = skeleton_interface.get_alveola(ind_alveola)?;
            let alveola_in = alveola.is_in()?;
            if !alveola.is_computed() && alveola_in {
                skeleton_operations::compute_alveola(&mut skeleton_interface, ind_alveola)?;
                let mut vec_neigh =
                    skeleton_operations::neighbor_alveolae(&mut skeleton_interface, ind_alveola)?;
                vec_alveola.append(&mut vec_neigh);
            }
            if alveola_in {
                skeleton_operations::include_alveola_in_skel(&mut skeleton_interface, ind_alveola)?;
            }
            print!("\r{} alveolae remaining     ", vec_alveola.len());
        } else {
            break;
        }
    }
    println!("");

    println!("Checking skeleton");
    skeleton_interface.check()?;

    Ok(())
}
