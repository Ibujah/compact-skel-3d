use nalgebra::base::*;

pub fn is_flat(pts: [Vector3<f32>; 4], eps: Option<f32>) -> bool {
    let eps_val = eps.unwrap_or(0.00001);
    let vec_3_0n = (pts[0] - pts[3]).normalize();
    let vec_3_1n = (pts[1] - pts[3]).normalize();
    let vec_3_2n = (pts[2] - pts[3]).normalize();

    #[rustfmt::skip]
    let mat_eval = Matrix3::new(
        vec_3_0n[0], vec_3_0n[1], vec_3_0n[2], 
        vec_3_1n[0], vec_3_1n[1], vec_3_1n[2], 
        vec_3_2n[0], vec_3_2n[1], vec_3_2n[2], 
    );

    let det = mat_eval.determinant().abs();

    det < eps_val
}

pub fn sphere_center(pts: [Vector3<f32>; 4]) -> Option<Vector3<f32>> {
    let vec_0_1 = pts[1] - pts[0];
    let vec_1_2 = pts[2] - pts[1];
    let vec_2_0 = pts[0] - pts[2];
    let vec_3_0 = pts[0] - pts[3];
    let vec_3_1 = pts[1] - pts[3];
    let vec_3_2 = pts[2] - pts[3];

    #[rustfmt::skip]
    let mat_slv = Matrix6x3::new(
        vec_0_1[0], vec_0_1[1], vec_0_1[2], 
        vec_1_2[0], vec_1_2[1], vec_1_2[2], 
        vec_2_0[0], vec_2_0[1], vec_2_0[2], 
        vec_3_0[0], vec_3_0[1], vec_3_0[2], 
        vec_3_1[0], vec_3_1[1], vec_3_1[2], 
        vec_3_2[0], vec_3_2[1], vec_3_2[2], 
    );

    let sqn0 = pts[0].norm_squared();
    let sqn1 = pts[1].norm_squared();
    let sqn2 = pts[2].norm_squared();
    let sqn3 = pts[3].norm_squared();

    let vec_slv = Matrix6x1::new(
        0.5 * (sqn1 - sqn0),
        0.5 * (sqn2 - sqn1),
        0.5 * (sqn0 - sqn2),
        0.5 * (sqn0 - sqn3),
        0.5 * (sqn1 - sqn3),
        0.5 * (sqn2 - sqn3),
    );
    let mat_slv_mod = mat_slv.transpose() * mat_slv;
    let vec_slv_mod = mat_slv.transpose() * vec_slv;

    mat_slv_mod.lu().solve(&vec_slv_mod)
}

pub fn circle_center(pts: [Vector3<f32>; 3]) -> Option<Vector3<f32>> {
    let vec_0_1 = pts[1] - pts[0];
    let vec_1_2 = pts[2] - pts[1];
    let vec_2_0 = pts[0] - pts[2];
    let vec_c = (pts[2] - pts[0]).cross(&(pts[1] - pts[0])).normalize();

    #[rustfmt::skip]
    let mat_slv = Matrix4x3::new(
        vec_0_1[0], vec_0_1[1], vec_0_1[2], 
        vec_1_2[0], vec_1_2[1], vec_1_2[2], 
        vec_2_0[0], vec_2_0[1], vec_2_0[2], 
        vec_c[0], vec_c[1], vec_c[2], 
    );

    #[rustfmt::skip]
    let vec_slv = Matrix4x1::new(
        0.5 * vec_0_1.norm_squared() + vec_0_1.dot(&pts[0]),
        0.5 * vec_1_2.norm_squared() + vec_1_2.dot(&pts[1]),
        0.5 * vec_2_0.norm_squared() + vec_2_0.dot(&pts[2]),
        vec_c.dot(&pts[0]),
    );

    let mat_slv_mod = mat_slv.transpose() * mat_slv;
    let vec_slv_mod = mat_slv.transpose() * vec_slv;

    mat_slv_mod.lu().solve(&vec_slv_mod)
}

pub fn center_and_radius(pts: [Vector3<f32>; 4], eps: Option<f32>) -> Option<(Vector3<f32>, f32)> {
    let center = if is_flat(pts, eps) {
        circle_center([pts[0], pts[1], pts[2]])
    } else {
        sphere_center(pts)
    };
    match center {
        Some(center) => {
            let radius = (center - pts[0]).norm();
            Some((center, radius))
        }
        None => None,
    }
}
