use anyhow::Result;
use compact_skel_3d::algorithm::skeleton_interface::{skeleton_operations, SkeletonInterface3D};
use compact_skel_3d::mesh3d::{self, Mesh3D};
use compact_skel_3d::skeleton3d::{self, Skeleton3D};
use nalgebra::base::*;

fn generate_test_mesh() -> Result<Mesh3D> {
    let mut mesh = Mesh3D::new();
    let up_vert = mesh.add_vertex(&Vector3::new(0.0, 0.0, 0.5));
    let down_vert = mesh.add_vertex(&Vector3::new(0.0, 0.0, -0.5));
    let mut surr_vert = Vec::new();
    for i in 0..9 {
        let ang = 2.0 * std::f32::consts::PI * (i as f32) / 9.0;
        let cur_vert = mesh.add_vertex(&Vector3::new(ang.cos(), ang.sin(), 0.0));
        surr_vert.push(cur_vert);
    }

    for i in 0..9 {
        let j = (i + 1) % 9;
        mesh.add_face(up_vert, surr_vert[i], surr_vert[j])?;
        mesh.add_face(down_vert, surr_vert[j], surr_vert[i])?;
    }

    Ok(mesh)
}

fn main() -> Result<()> {
    // let mut mesh = generate_test_mesh()?;

    let mut mesh = mesh3d::io::load_obj("./ressources/hand_del.obj")?;
    mesh.check_mesh()?;

    let mut skeleton = Skeleton3D::new();

    println!("Init skeleton interface");
    let mut skeleton_interface = SkeletonInterface3D::init(&mut mesh, &mut skeleton)?;

    println!("Finding first node");
    let ind_first_node = skeleton_operations::first_node_in(&mut skeleton_interface)?;

    println!("Checking skeleton");
    skeleton_interface.check()?;

    println!("Finding alveola");
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
            print!("\r{} alveolae remaining     ", vec_alveola.len());
            if !skeleton_interface.get_alveola(ind_alveola)?.is_computed()
                && skeleton_interface.get_alveola(ind_alveola)?.is_in()?
            {
                skeleton_operations::compute_alveola(&mut skeleton_interface, ind_alveola)?;
                skeleton_operations::include_alveola_in_skel(&mut skeleton_interface, ind_alveola)?;
                let mut vec_neigh =
                    skeleton_operations::neighbor_alveolae(&mut skeleton_interface, ind_alveola)?;
                vec_alveola.append(&mut vec_neigh);
            }
        } else {
            break;
        }
    }
    println!("");

    println!("Checking skeleton");
    skeleton_interface.check()?;

    println!("Saving skeleton and mesh");
    mesh3d::io::save_obj("./ressources/mesh.obj", skeleton_interface.get_mesh())?;

    skeleton3d::io::save_obj(
        "./ressources/skeleton.obj",
        skeleton_interface.get_skeleton(),
    )?;

    Ok(())
}
