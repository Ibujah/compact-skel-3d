use nalgebra::base::*;

pub fn center_and_radius(pts: [Vector3<f32>; 4]) -> Option<(Vector3<f32>, f32)> {
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

    match mat_slv_mod.lu().solve(&vec_slv_mod) {
        Some(center) => {
            let radius = (center - pts[0]).norm();
            Some((center, radius))
        }
        None => None,
    }
}
