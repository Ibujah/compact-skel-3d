use nalgebra::Vector3;

use crate::mesh3d::ManifoldMesh3D;
use crate::skeleton3d::Skeleton3D;

/// Computes hausdorff distance between a skeleton and a mesh
pub fn hausdorff_distance(mesh: &ManifoldMesh3D, skel: &Skeleton3D) -> f64 {
    let mut dmax = 0.0;
    let mut min_vert: Option<Vector3<f64>> = None;
    let mut max_vert: Option<Vector3<f64>> = None;
    for vert in mesh.vertices().iter() {
        min_vert = if let Some(v) = min_vert {
            let x = if v[0] < vert[0] { v[0] } else { vert[0] };
            let y = if v[1] < vert[1] { v[1] } else { vert[1] };
            let z = if v[2] < vert[2] { v[2] } else { vert[2] };
            Some(Vector3::new(x, y, z))
        } else {
            Some(*vert)
        };

        max_vert = if let Some(v) = max_vert {
            let x = if v[0] > vert[0] { v[0] } else { vert[0] };
            let y = if v[1] > vert[1] { v[1] } else { vert[1] };
            let z = if v[2] > vert[2] { v[2] } else { vert[2] };
            Some(Vector3::new(x, y, z))
        } else {
            Some(*vert)
        };

        let dmin = skel
            .get_nodes()
            .iter()
            .map(|(_, sph)| {
                let dist = (vert - sph.center).norm() - sph.radius;
                if dist > 0.0 {
                    dist
                } else {
                    0.0
                }
            })
            .fold(None, |vmin_opt, vcur| {
                if let Some(vmin) = vmin_opt {
                    if vmin < vcur {
                        Some(vmin)
                    } else {
                        Some(vcur)
                    }
                } else {
                    Some(vcur)
                }
            })
            .unwrap();
        if dmin > dmax {
            dmax = dmin
        }
    }

    let diag_length = (min_vert.unwrap() - max_vert.unwrap()).norm();

    dmax / diag_length
}
