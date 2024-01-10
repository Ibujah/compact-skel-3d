use crate::mesh3d::ManifoldMesh3D;
use anyhow::Result;
use nalgebra::base::*;

/// Checks if a halfedge can be flipped
pub fn can_flip_halfedge(mesh: &ManifoldMesh3D, ind_halfedge: usize) -> Result<bool> {
    let halfedge = mesh.get_halfedge(ind_halfedge)?;

    let opp_vert1 = halfedge
        .next_halfedge()
        .ok_or(anyhow::Error::msg(
            "can_flip_halfedge(): Halfedge should have next",
        ))?
        .last_vertex();
    let opp_vert2 = halfedge
        .opposite_halfedge()
        .ok_or(anyhow::Error::msg(
            "can_flip_halfedge(): Halfedge should have opposite",
        ))?
        .next_halfedge()
        .ok_or(anyhow::Error::msg(
            "can_flip_halfedge(): Opposite halfedge should have next",
        ))?
        .last_vertex();

    let edg_found = opp_vert1.halfedges().iter().fold(false, |res, &he| {
        res || he.last_vertex().ind() == opp_vert2.ind()
    });

    Ok(!edg_found)
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
    if !mesh.halfedges.contains_key(&ind_halfedge) {
        return Err(anyhow::Error::msg("flip_halfedge(): Index out of bounds"));
    }
    let flippable = can_flip_halfedge(mesh, ind_halfedge)?;

    if !flippable {
        return Ok(false);
    }

    let ind_he_12 = ind_halfedge;
    let &ind_fac_123 = mesh
        .map_hedg_face
        .get(&ind_he_12)
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty face"))?;
    let &ind_he_23 = mesh
        .map_hedg_next
        .get(&ind_he_12)
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;

    let &ind_he_21 = mesh
        .map_hedg_opp
        .get(&ind_he_12)
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
    let &ind_fac_214 = mesh
        .map_hedg_face
        .get(&ind_he_21)
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty face"))?;
    let &ind_he_14 = mesh
        .map_hedg_next
        .get(&ind_he_21)
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;

    let ind_v1 = mesh.halfedges.get(&ind_he_12).unwrap()[0];
    let ind_v2 = mesh.halfedges.get(&ind_he_12).unwrap()[1];
    let ind_v3 = mesh.halfedges.get(&ind_he_23).unwrap()[1];
    let ind_v4 = mesh.halfedges.get(&ind_he_14).unwrap()[1];

    mesh.remove_face(ind_fac_214)?;
    mesh.remove_face(ind_fac_123)?;
    mesh.add_face(ind_v1, ind_v4, ind_v3)?;
    mesh.add_face(ind_v2, ind_v3, ind_v4)?;

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
    vert: &Vector3<f32>,
    ind_halfedge: usize,
) -> Result<usize> {
    if !mesh.halfedges.contains_key(&ind_halfedge) {
        return Err(anyhow::Error::msg("split_halfedge(): Index out of bounds"));
    }

    let ind_he_12 = ind_halfedge;
    let &ind_fac_123 = mesh
        .map_hedg_face
        .get(&ind_he_12)
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty face"))?;
    let &ind_he_23 = mesh
        .map_hedg_next
        .get(&ind_he_12)
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;

    let &ind_he_21 = mesh
        .map_hedg_opp
        .get(&ind_he_12)
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
    let &ind_fac_214 = mesh
        .map_hedg_face
        .get(&ind_he_21)
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty face"))?;
    let &ind_he_14 = mesh
        .map_hedg_next
        .get(&ind_he_21)
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;

    let ind_v1 = mesh.halfedges.get(&ind_he_12).unwrap()[0];
    let ind_v2 = mesh.halfedges.get(&ind_he_12).unwrap()[1];
    let ind_v3 = mesh.halfedges.get(&ind_he_23).unwrap()[1];
    let ind_v4 = mesh.halfedges.get(&ind_he_14).unwrap()[1];
    let ind_v5 = mesh.add_vertex(vert);

    mesh.remove_face(ind_fac_123)?;
    mesh.remove_face(ind_fac_214)?;
    mesh.add_face(ind_v1, ind_v5, ind_v3)?;
    mesh.add_face(ind_v2, ind_v3, ind_v5)?;
    mesh.add_face(ind_v2, ind_v5, ind_v4)?;
    mesh.add_face(ind_v1, ind_v4, ind_v5)?;

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
    vert: &Vector3<f32>,
    ind_face: usize,
) -> Result<usize> {
    if !mesh.faces.contains_key(&ind_face) {
        return Err(anyhow::Error::msg("split_face(): Index out of bounds"));
    }

    let ind_fac_123 = ind_face;
    let [ind_he_12, ind_he_23, ind_he_31] = mesh.get_face(ind_face)?.face_halfedges();

    let ind_v1 = mesh.halfedges.get(&ind_he_12).unwrap()[0];
    let ind_v2 = mesh.halfedges.get(&ind_he_23).unwrap()[0];
    let ind_v3 = mesh.halfedges.get(&ind_he_31).unwrap()[0];
    let ind_v4 = mesh.add_vertex(vert);

    mesh.remove_face(ind_fac_123)?;
    mesh.add_face(ind_v1, ind_v2, ind_v4)?;
    mesh.add_face(ind_v2, ind_v3, ind_v4)?;
    mesh.add_face(ind_v1, ind_v4, ind_v3)?;

    Ok(ind_v4)
}
