use super::manifold_mesh3d::{IterHalfEdge, ManifoldMesh3D};
use anyhow::Result;
use nalgebra::base::*;

/// Checks if a halfedge can be flipped
fn can_flip_halfedge(halfedge: &IterHalfEdge) -> bool {
    let opp_vert1 = halfedge.next_halfedge().last_vertex();
    let opp_vert2 = halfedge
        .opposite_halfedge()
        .unwrap()
        .next_halfedge()
        .last_vertex();

    let edg_found = opp_vert1.halfedges().iter().fold(false, |res, &he| {
        res || he.last_vertex().ind() == opp_vert2.ind()
    });

    !edg_found
}

/// Flips a halfedge
///
/// Given halfedge (1->2):
/// ```text
///     1             1
///   / | \         /   \
///  4  |  3  -->  4 --- 3
///   \ | /         \   /
///     2             2
/// ```
pub fn flip_halfedge(mesh: &mut ManifoldMesh3D, ind_halfedge: usize) -> Result<bool> {
    // Get the halfedge to flip
    let ind_he_12 = ind_halfedge;
    let he12 = mesh.get_halfedge(ind_he_12)?;

    // Check if we can flip this halfedge
    if !can_flip_halfedge(&he12) {
        return Ok(false);
    }

    // Get the next and previous halfedges of he12
    let he23 = he12.next_halfedge();
    let he31 = he12.prev_halfedge();

    // Get the opposite halfedge of he12
    let he21 = he12.opposite_halfedge().unwrap();

    // Get the next and previous halfedges of he21
    let he14 = he21.next_halfedge();
    let he42 = he21.prev_halfedge();

    // Get the halfedges indices
    let ind_he_23 = he23.ind();
    let ind_he_31 = he31.ind();
    let ind_he_21 = he21.ind();
    let ind_he_14 = he14.ind();
    let ind_he_42 = he42.ind();

    // Get the indices of the faces attached to the halfedges
    let ind_fa_123 = ind_he_12 - ind_he_12 % 3;
    let ind_fa_214 = ind_he_21 - ind_he_21 % 3;

    // Get the indices of non modified halfedges
    let ind_he_32 = he23.opposite_halfedge().unwrap().ind();
    let ind_he_13 = he31.opposite_halfedge().unwrap().ind();
    let ind_he_41 = he14.opposite_halfedge().unwrap().ind();
    let ind_he_24 = he42.opposite_halfedge().unwrap().ind();

    // Get the indices of the vertices attached to the halfedges
    let ind_v1 = mesh.hedg_vert_inds[ind_he_12];
    let ind_v2 = mesh.hedg_vert_inds[ind_he_21];
    let ind_v3 = mesh.hedg_vert_inds[ind_he_32];
    let ind_v4 = mesh.hedg_vert_inds[ind_he_41];

    // Delete removed halfedges from vertex-halfedge mappings
    mesh.vert_hedg[ind_v1].retain(|&ind_he| ind_he != ind_he_14 && ind_he != ind_he_12);
    mesh.vert_hedg[ind_v2].retain(|&ind_he| ind_he != ind_he_23 && ind_he != ind_he_21);
    mesh.vert_hedg[ind_v3].retain(|&ind_he| ind_he != ind_he_31);
    mesh.vert_hedg[ind_v4].retain(|&ind_he| ind_he != ind_he_42);

    // Face renaming
    let ind_fa_143 = ind_fa_123;
    let ind_fa_234 = ind_fa_214;

    // Added halfedges
    let ind_he_14 = ind_fa_143;
    let ind_he_43 = ind_fa_143 + 1;
    let ind_he_31 = ind_fa_143 + 2;
    let ind_he_23 = ind_fa_234;
    let ind_he_34 = ind_fa_234 + 1;
    let ind_he_42 = ind_fa_234 + 2;

    // Update the halfedge-vertex mappings
    mesh.hedg_vert_inds[ind_he_14] = ind_v1;
    mesh.hedg_vert_inds[ind_he_43] = ind_v4;
    mesh.hedg_vert_inds[ind_he_31] = ind_v3;
    mesh.hedg_vert_inds[ind_he_23] = ind_v2;
    mesh.hedg_vert_inds[ind_he_34] = ind_v3;
    mesh.hedg_vert_inds[ind_he_42] = ind_v4;

    // Update opposite halfedges mappings
    mesh.hedg_opp[ind_he_32] = Some(ind_he_23);
    mesh.hedg_opp[ind_he_13] = Some(ind_he_31);
    mesh.hedg_opp[ind_he_41] = Some(ind_he_14);
    mesh.hedg_opp[ind_he_24] = Some(ind_he_42);
    mesh.hedg_opp[ind_he_14] = Some(ind_he_41);
    mesh.hedg_opp[ind_he_43] = Some(ind_he_34);
    mesh.hedg_opp[ind_he_31] = Some(ind_he_13);
    mesh.hedg_opp[ind_he_23] = Some(ind_he_32);
    mesh.hedg_opp[ind_he_34] = Some(ind_he_43);
    mesh.hedg_opp[ind_he_42] = Some(ind_he_24);

    // Insert added halfedges to vertex-halfedge mappings
    mesh.vert_hedg[ind_v1].push(ind_he_14);
    mesh.vert_hedg[ind_v2].push(ind_he_23);
    mesh.vert_hedg[ind_v3].push(ind_he_31);
    mesh.vert_hedg[ind_v3].push(ind_he_34);
    mesh.vert_hedg[ind_v4].push(ind_he_42);
    mesh.vert_hedg[ind_v4].push(ind_he_43);

    Ok(true)
}

/// Splits an halfedge
///
/// Given halfedge (1->2):
/// ```text
///     1             1
///   / | \         / | \
///  4  |  3  -->  4 -5- 3
///   \ | /         \ | /
///     2             2
/// ```
pub fn split_halfedge(
    mesh: &mut ManifoldMesh3D,
    vert: &Vector3<f64>,
    ind_halfedge: usize,
) -> Result<usize> {
    // Get the halfedge to remove
    let ind_he_12 = ind_halfedge;
    let he12 = mesh.get_halfedge(ind_he_12)?;

    // Get the next and previous halfedges of he12
    let he23 = he12.next_halfedge();
    let he31 = he12.prev_halfedge();

    // Get the opposite halfedge of he12
    let he21 = he12.opposite_halfedge().unwrap();

    // Get the next and previous halfedges of he21
    let he14 = he21.next_halfedge();
    let he42 = he21.prev_halfedge();

    // Get the indices of halfedges to remove
    let ind_he_23 = he23.ind();
    let ind_he_31 = he31.ind();
    let ind_he_21 = he21.ind();
    let ind_he_14 = he14.ind();
    let ind_he_42 = he42.ind();

    // Get the indices of the faces attached to the halfedges
    let ind_fa_123 = he12.face().ind() * 3;
    let ind_fa_214 = he21.face().ind() * 3;

    // Get the indices of non modified halfedges
    let ind_he_32 = he23.opposite_halfedge().unwrap().ind();
    let ind_he_13 = he31.opposite_halfedge().unwrap().ind();
    let ind_he_41 = he14.opposite_halfedge().unwrap().ind();
    let ind_he_24 = he42.opposite_halfedge().unwrap().ind();

    // Get the indices of the vertices attached to the halfedges
    let ind_v1 = mesh.hedg_vert_inds[ind_he_12];
    let ind_v2 = mesh.hedg_vert_inds[ind_he_21];
    let ind_v3 = mesh.hedg_vert_inds[ind_he_32];
    let ind_v4 = mesh.hedg_vert_inds[ind_he_41];
    let ind_v5 = mesh.add_vertex(vert);

    // Delete removed halfedges from vertex-halfedge mappings
    mesh.vert_hedg[ind_v1].retain(|&ind_he| ind_he != ind_he_14 && ind_he != ind_he_12);
    mesh.vert_hedg[ind_v2].retain(|&ind_he| ind_he != ind_he_23 && ind_he != ind_he_21);
    mesh.vert_hedg[ind_v3].retain(|&ind_he| ind_he != ind_he_31);
    mesh.vert_hedg[ind_v4].retain(|&ind_he| ind_he != ind_he_42);

    // Added faces
    let ind_fa_153 = ind_fa_123;
    let ind_fa_254 = ind_fa_214;
    let ind_fa_352 = mesh.hedg_vert_inds.len();
    let ind_fa_451 = ind_fa_352 + 3;

    // Added halfedges
    let ind_he_15 = ind_fa_153;
    let ind_he_53 = ind_fa_153 + 1;
    let ind_he_31 = ind_fa_153 + 2;

    let ind_he_25 = ind_fa_254;
    let ind_he_54 = ind_fa_254 + 1;
    let ind_he_42 = ind_fa_254 + 2;

    let ind_he_35 = ind_fa_352;
    let ind_he_52 = ind_fa_352 + 1;
    let ind_he_23 = ind_fa_352 + 2;

    let ind_he_45 = ind_fa_451;
    let ind_he_51 = ind_fa_451 + 1;
    let ind_he_14 = ind_fa_451 + 2;

    // Update the halfedge-vertex mappings
    mesh.hedg_vert_inds[ind_he_15] = ind_v1;
    mesh.hedg_vert_inds[ind_he_53] = ind_v5;
    mesh.hedg_vert_inds[ind_he_31] = ind_v3;

    mesh.hedg_vert_inds[ind_he_25] = ind_v2;
    mesh.hedg_vert_inds[ind_he_54] = ind_v5;
    mesh.hedg_vert_inds[ind_he_42] = ind_v4;

    mesh.hedg_vert_inds.push(ind_v3);
    mesh.hedg_vert_inds.push(ind_v5);
    mesh.hedg_vert_inds.push(ind_v2);

    mesh.hedg_vert_inds.push(ind_v4);
    mesh.hedg_vert_inds.push(ind_v5);
    mesh.hedg_vert_inds.push(ind_v1);

    // Update opposite halfedges mappings
    mesh.hedg_opp[ind_he_32] = Some(ind_he_23);
    mesh.hedg_opp[ind_he_13] = Some(ind_he_31);
    mesh.hedg_opp[ind_he_41] = Some(ind_he_14);
    mesh.hedg_opp[ind_he_24] = Some(ind_he_42);

    mesh.hedg_opp[ind_he_15] = Some(ind_he_51);
    mesh.hedg_opp[ind_he_53] = Some(ind_he_35);
    mesh.hedg_opp[ind_he_31] = Some(ind_he_13);

    mesh.hedg_opp[ind_he_25] = Some(ind_he_52);
    mesh.hedg_opp[ind_he_54] = Some(ind_he_45);
    mesh.hedg_opp[ind_he_42] = Some(ind_he_24);

    mesh.hedg_opp.push(Some(ind_he_53));
    mesh.hedg_opp.push(Some(ind_he_25));
    mesh.hedg_opp.push(Some(ind_he_32));

    mesh.hedg_opp.push(Some(ind_he_54));
    mesh.hedg_opp.push(Some(ind_he_15));
    mesh.hedg_opp.push(Some(ind_he_41));

    // Insert added halfedges to vertex-halfedge mappings
    mesh.vert_hedg[ind_v1].push(ind_he_14);
    mesh.vert_hedg[ind_v1].push(ind_he_15);
    mesh.vert_hedg[ind_v2].push(ind_he_23);
    mesh.vert_hedg[ind_v2].push(ind_he_25);
    mesh.vert_hedg[ind_v3].push(ind_he_31);
    mesh.vert_hedg[ind_v3].push(ind_he_35);
    mesh.vert_hedg[ind_v4].push(ind_he_42);
    mesh.vert_hedg[ind_v4].push(ind_he_45);
    mesh.vert_hedg[ind_v5].push(ind_he_51);
    mesh.vert_hedg[ind_v5].push(ind_he_52);
    mesh.vert_hedg[ind_v5].push(ind_he_53);
    mesh.vert_hedg[ind_v5].push(ind_he_54);

    mesh.face_groups.push(None);
    mesh.face_groups.push(None);

    Ok(ind_v5)
}

/// Splits a face
///
/// Given face (1, 2, 3):
/// ```text
///      1                1
///    /   \            / | \
///   /     \   -->    /  4  \
///  /       \        / /   \ \
/// 2 ------- 3      2 ------- 3
/// ```
pub fn split_face(
    mesh: &mut ManifoldMesh3D,
    vert: &Vector3<f64>,
    ind_face: usize,
) -> Result<usize> {
    // Get face to remove
    let fac_123 = mesh.get_face(ind_face)?;
    let ind_fa_123 = ind_face * 3;

    // Get halfedges surrouding face
    let [he12, he23, he31] = fac_123.halfedges();

    // Get the indices of the vertices attached to the halfedges
    let ind_v1 = he12.first_vertex().ind();
    let ind_v2 = he23.first_vertex().ind();
    let ind_v3 = he31.first_vertex().ind();

    // Get the indices of halfedges to remove
    let ind_he_12 = he12.ind();
    let ind_he_23 = he23.ind();
    let ind_he_31 = he31.ind();

    // Get the indices of non modified halfedges
    let ind_he_21 = he12.opposite_halfedge().unwrap().ind();
    let ind_he_32 = he23.opposite_halfedge().unwrap().ind();
    let ind_he_13 = he31.opposite_halfedge().unwrap().ind();

    // Add new vertex
    let ind_v4 = mesh.add_vertex(vert);

    // Delete removed halfedges from vertex-halfedge mappings
    mesh.vert_hedg[ind_v1].retain(|&ind_he| ind_he != ind_he_12);
    mesh.vert_hedg[ind_v2].retain(|&ind_he| ind_he != ind_he_23);
    mesh.vert_hedg[ind_v3].retain(|&ind_he| ind_he != ind_he_31);

    // Added faces
    let ind_fa_124 = ind_fa_123;
    let ind_fa_234 = mesh.hedg_vert_inds.len();
    let ind_fa_314 = ind_fa_234 + 3;

    // Added halfedges
    let ind_he_12 = ind_fa_124;
    let ind_he_24 = ind_fa_124 + 1;
    let ind_he_41 = ind_fa_124 + 2;

    let ind_he_23 = ind_fa_234;
    let ind_he_34 = ind_fa_234 + 1;
    let ind_he_42 = ind_fa_234 + 2;

    let ind_he_31 = ind_fa_314;
    let ind_he_14 = ind_fa_314 + 1;
    let ind_he_43 = ind_fa_314 + 2;

    // Update the halfedge-vertex mappings
    mesh.hedg_vert_inds[ind_he_12] = ind_v1;
    mesh.hedg_vert_inds[ind_he_24] = ind_v2;
    mesh.hedg_vert_inds[ind_he_41] = ind_v4;

    mesh.hedg_vert_inds.push(ind_v2);
    mesh.hedg_vert_inds.push(ind_v3);
    mesh.hedg_vert_inds.push(ind_v4);

    mesh.hedg_vert_inds.push(ind_v3);
    mesh.hedg_vert_inds.push(ind_v1);
    mesh.hedg_vert_inds.push(ind_v4);

    // Update opposite halfedges mappings
    mesh.hedg_opp[ind_he_21] = Some(ind_he_12);
    mesh.hedg_opp[ind_he_32] = Some(ind_he_23);
    mesh.hedg_opp[ind_he_13] = Some(ind_he_31);

    mesh.hedg_opp[ind_he_12] = Some(ind_he_21);
    mesh.hedg_opp[ind_he_24] = Some(ind_he_42);
    mesh.hedg_opp[ind_he_41] = Some(ind_he_14);

    mesh.hedg_opp.push(Some(ind_he_32));
    mesh.hedg_opp.push(Some(ind_he_43));
    mesh.hedg_opp.push(Some(ind_he_24));

    mesh.hedg_opp.push(Some(ind_he_13));
    mesh.hedg_opp.push(Some(ind_he_41));
    mesh.hedg_opp.push(Some(ind_he_34));

    // Insert added halfedges to vertex-halfedge mappings
    mesh.vert_hedg[ind_v1].push(ind_he_12);
    mesh.vert_hedg[ind_v1].push(ind_he_14);
    mesh.vert_hedg[ind_v2].push(ind_he_23);
    mesh.vert_hedg[ind_v2].push(ind_he_24);
    mesh.vert_hedg[ind_v3].push(ind_he_31);
    mesh.vert_hedg[ind_v3].push(ind_he_34);
    mesh.vert_hedg[ind_v4].push(ind_he_41);
    mesh.vert_hedg[ind_v4].push(ind_he_42);
    mesh.vert_hedg[ind_v4].push(ind_he_43);

    mesh.face_groups.push(None);
    mesh.face_groups.push(None);

    Ok(ind_v4)
}
