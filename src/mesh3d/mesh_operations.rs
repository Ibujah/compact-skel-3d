use crate::mesh3d::Mesh3D;
use anyhow::Result;
use nalgebra::base::*;

pub fn can_flip_halfedge(mesh: &Mesh3D, ind_halfedge: usize) -> Result<bool> {
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

pub fn flip_halfedge(mesh: &mut Mesh3D, ind_halfedge: usize) -> Result<bool> {
    if ind_halfedge >= mesh.halfedges.len() {
        return Err(anyhow::Error::msg("flip_halfedge(): Index out of bounds"));
    }
    let flippable = can_flip_halfedge(mesh, ind_halfedge)?;

    if !flippable {
        return Ok(false);
    }

    // direct order
    //     1             1
    //   / | \         /   \
    //  4  |  3  -->  4 --- 3
    //   \ | /         \   /
    //     2             2

    let ind_he_12 = ind_halfedge;
    let ind_fac_123 =
        mesh.map_hedg_face[ind_he_12].ok_or(anyhow::Error::msg("flip_halfedge(): Empty face"))?;
    let ind_he_23 = mesh.map_hedg_next[ind_he_12]
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
    let ind_he_31 = mesh.map_hedg_prev[ind_he_12]
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;

    let ind_he_21 = mesh.map_hedg_opp[ind_he_12]
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
    let ind_fac_214 =
        mesh.map_hedg_face[ind_he_21].ok_or(anyhow::Error::msg("flip_halfedge(): Empty face"))?;
    let ind_he_14 = mesh.map_hedg_next[ind_he_21]
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
    let ind_he_42 = mesh.map_hedg_prev[ind_he_21]
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;

    let ind_v1 = mesh.halfedges[ind_he_12][0];
    let ind_v2 = mesh.halfedges[ind_he_12][1];
    let ind_v3 = mesh.halfedges[ind_he_23][1];
    let ind_v4 = mesh.halfedges[ind_he_14][1];

    // rename edges and faces
    let ind_fac_143 = ind_fac_123;
    let ind_fac_234 = ind_fac_214;

    let ind_he_43 = ind_he_12;
    let ind_he_34 = ind_he_21;

    // update faces
    mesh.faces[ind_fac_143] = [ind_he_14, ind_he_43, ind_he_31];
    mesh.faces[ind_fac_234] = [ind_he_23, ind_he_34, ind_he_42];

    // update edges
    mesh.halfedges[ind_he_43] = [ind_v4, ind_v3];
    mesh.halfedges[ind_he_34] = [ind_v3, ind_v4];

    // update connectivity
    mesh.map_hedg_face[ind_he_14] = Some(ind_fac_143);
    mesh.map_hedg_face[ind_he_23] = Some(ind_fac_234);

    mesh.map_hedg_next[ind_he_43] = Some(ind_he_31);
    mesh.map_hedg_next[ind_he_31] = Some(ind_he_14);
    mesh.map_hedg_next[ind_he_14] = Some(ind_he_43);

    mesh.map_hedg_prev[ind_he_43] = Some(ind_he_14);
    mesh.map_hedg_prev[ind_he_14] = Some(ind_he_31);
    mesh.map_hedg_prev[ind_he_31] = Some(ind_he_43);

    mesh.map_hedg_next[ind_he_34] = Some(ind_he_42);
    mesh.map_hedg_next[ind_he_42] = Some(ind_he_23);
    mesh.map_hedg_next[ind_he_23] = Some(ind_he_34);

    mesh.map_hedg_prev[ind_he_34] = Some(ind_he_23);
    mesh.map_hedg_prev[ind_he_23] = Some(ind_he_42);
    mesh.map_hedg_prev[ind_he_42] = Some(ind_he_34);

    mesh.map_vert_hedg[ind_v1].retain(|&ind| ind != ind_he_12);
    mesh.map_vert_hedg[ind_v2].retain(|&ind| ind != ind_he_21);
    mesh.map_vert_hedg[ind_v3].push(ind_he_34);
    mesh.map_vert_hedg[ind_v4].push(ind_he_43);

    Ok(true)
}

pub fn split_halfedge(mesh: &mut Mesh3D, vert: &Vector3<f32>, ind_halfedge: usize) -> Result<()> {
    if ind_halfedge >= mesh.halfedges.len() {
        return Err(anyhow::Error::msg("split_halfedge(): Index out of bounds"));
    }

    // direct order
    //     1             1
    //   / | \         / | \
    //  4  |  3  -->  4 -5- 3
    //   \ | /         \ | /
    //     2             2

    let ind_he_12 = ind_halfedge;
    let ind_fac_123 =
        mesh.map_hedg_face[ind_he_12].ok_or(anyhow::Error::msg("flip_halfedge(): Empty face"))?;
    let ind_he_23 = mesh.map_hedg_next[ind_he_12]
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
    let ind_he_31 = mesh.map_hedg_prev[ind_he_12]
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;

    let ind_he_21 = mesh.map_hedg_opp[ind_he_12]
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
    let ind_fac_214 =
        mesh.map_hedg_face[ind_he_21].ok_or(anyhow::Error::msg("flip_halfedge(): Empty face"))?;
    let ind_he_14 = mesh.map_hedg_next[ind_he_21]
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
    let ind_he_42 = mesh.map_hedg_prev[ind_he_21]
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;

    let ind_he_13 = mesh.map_hedg_opp[ind_he_31]
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
    let ind_he_32 = mesh.map_hedg_opp[ind_he_23]
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
    let ind_he_24 = mesh.map_hedg_opp[ind_he_42]
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;
    let ind_he_41 = mesh.map_hedg_opp[ind_he_14]
        .ok_or(anyhow::Error::msg("flip_halfedge(): Empty halfedge"))?;

    let ind_v1 = mesh.halfedges[ind_he_12][0];
    let ind_v2 = mesh.halfedges[ind_he_12][1];
    let ind_v3 = mesh.halfedges[ind_he_23][1];
    let ind_v4 = mesh.halfedges[ind_he_14][1];
    let ind_v5 = mesh.add_vertex(vert);

    mesh.map_vert_hedg[ind_v1].retain(|&ind| ind != ind_he_12);
    mesh.map_vert_hedg[ind_v2].retain(|&ind| ind != ind_he_21);

    // rename edges
    let ind_he_52 = ind_he_12;
    let ind_he_51 = ind_he_21;

    // reinit edges
    mesh.halfedges[ind_he_51] = [ind_v5, ind_v1];
    mesh.halfedges[ind_he_52] = [ind_v5, ind_v2];

    mesh.map_hedg_face[ind_he_51] = None;
    mesh.map_hedg_face[ind_he_52] = None;

    mesh.map_hedg_prev[ind_he_51] = None;
    mesh.map_hedg_prev[ind_he_52] = None;

    mesh.map_hedg_next[ind_he_51] = None;
    mesh.map_hedg_next[ind_he_52] = None;

    mesh.map_hedg_opp[ind_he_51] = None;
    mesh.map_hedg_opp[ind_he_52] = None;

    mesh.map_vert_hedg[ind_v5].push(ind_he_51);
    mesh.map_vert_hedg[ind_v5].push(ind_he_52);

    // add new halfedges
    let ind_he_53 = mesh.add_halfedge(ind_v5, ind_v3)?;
    let ind_he_54 = mesh.add_halfedge(ind_v5, ind_v4)?;
    let ind_he_15 = mesh.add_halfedge(ind_v1, ind_v5)?;
    let ind_he_25 = mesh.add_halfedge(ind_v2, ind_v5)?;
    let ind_he_35 = mesh.add_halfedge(ind_v3, ind_v5)?;
    let ind_he_45 = mesh.add_halfedge(ind_v4, ind_v5)?;

    // rename faces
    let ind_fac_235 = ind_fac_123;
    let ind_fac_145 = ind_fac_214;

    // reinit faces
    mesh.faces[ind_fac_235] = [ind_he_23, ind_he_35, ind_he_52];
    mesh.faces[ind_fac_145] = [ind_he_14, ind_he_45, ind_he_51];

    // add new faces
    mesh.faces.push([ind_he_15, ind_he_53, ind_he_31]);
    let ind_fac_153 = mesh.faces.len() - 1;

    mesh.faces.push([ind_he_25, ind_he_54, ind_he_42]);
    let ind_fac_254 = mesh.faces.len() - 1;

    mesh.fill_face(
        ind_fac_235,
        ind_he_23,
        ind_he_35,
        ind_he_52,
        ind_he_32,
        ind_he_53,
        ind_he_25,
    );
    mesh.fill_face(
        ind_fac_145,
        ind_he_14,
        ind_he_45,
        ind_he_51,
        ind_he_41,
        ind_he_54,
        ind_he_15,
    );
    mesh.fill_face(
        ind_fac_153,
        ind_he_15,
        ind_he_53,
        ind_he_31,
        ind_he_51,
        ind_he_35,
        ind_he_13,
    );
    mesh.fill_face(
        ind_fac_254,
        ind_he_25,
        ind_he_54,
        ind_he_42,
        ind_he_52,
        ind_he_45,
        ind_he_24,
    );

    Ok(())
}

pub fn split_face(mesh: &mut Mesh3D, vert: &Vector3<f32>, ind_face: usize) -> Result<()> {
    if ind_face >= mesh.faces.len() {
        return Err(anyhow::Error::msg("split_face(): Index out of bounds"));
    }

    // direct order
    //      1                1
    //    /   \            / | \
    //   /     \   -->    /  4  \
    //  /       \        / /   \ \
    // 2 ------- 3      2 ------- 3

    let ind_fac_123 = ind_face;
    let [ind_he_12, ind_he_23, ind_he_31] = mesh.get_face(ind_face)?.face_halfedges();
    let ind_he_21 =
        mesh.map_hedg_opp[ind_he_12].ok_or(anyhow::Error::msg("split_face(): Empty halfedge"))?;
    let ind_he_32 =
        mesh.map_hedg_opp[ind_he_23].ok_or(anyhow::Error::msg("split_face(): Empty halfedge"))?;
    let ind_he_13 =
        mesh.map_hedg_opp[ind_he_31].ok_or(anyhow::Error::msg("split_face(): Empty halfedge"))?;

    let ind_v1 = mesh.halfedges[ind_he_12][0];
    let ind_v2 = mesh.halfedges[ind_he_23][0];
    let ind_v3 = mesh.halfedges[ind_he_31][0];
    let ind_v4 = mesh.add_vertex(vert);

    // rename face
    let ind_fac_124 = ind_fac_123;

    // add halfedges
    let ind_he_14 = mesh.add_halfedge(ind_v1, ind_v4)?;
    let ind_he_41 = mesh.add_halfedge(ind_v4, ind_v1)?;
    let ind_he_24 = mesh.add_halfedge(ind_v2, ind_v4)?;
    let ind_he_42 = mesh.add_halfedge(ind_v4, ind_v2)?;
    let ind_he_34 = mesh.add_halfedge(ind_v3, ind_v4)?;
    let ind_he_43 = mesh.add_halfedge(ind_v4, ind_v3)?;

    // reinit faces
    mesh.faces[ind_fac_124] = [ind_he_12, ind_he_24, ind_he_41];

    // add faces
    mesh.faces.push([ind_he_23, ind_he_34, ind_he_42]);
    let ind_fac_234 = mesh.faces.len() - 1;
    mesh.faces.push([ind_he_14, ind_he_43, ind_he_31]);
    let ind_fac_143 = mesh.faces.len() - 1;

    // fill faces
    mesh.fill_face(
        ind_fac_124,
        ind_he_12,
        ind_he_24,
        ind_he_41,
        ind_he_21,
        ind_he_42,
        ind_he_14,
    );

    mesh.fill_face(
        ind_fac_234,
        ind_he_23,
        ind_he_34,
        ind_he_42,
        ind_he_32,
        ind_he_43,
        ind_he_24,
    );

    mesh.fill_face(
        ind_fac_143,
        ind_he_14,
        ind_he_43,
        ind_he_31,
        ind_he_41,
        ind_he_34,
        ind_he_13,
    );

    Ok(())
}
