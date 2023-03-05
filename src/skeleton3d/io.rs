use anyhow::Result;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;

use crate::skeleton3d::Skeleton3D;

pub fn save_obj(filename: &str, skeleton: &Skeleton3D) -> Result<()> {
    let mut file = File::create(filename)?;

    let mut skel_ind_to_ind = HashMap::new();
    let mut ind = 1;
    for (skel_ind, sph) in skeleton.nodes.iter() {
        let vert = sph.center;
        writeln!(file, "v {} {} {}", vert[0], vert[1], vert[2])?;
        skel_ind_to_ind.insert(skel_ind, ind);
        ind = ind + 1;
    }

    for (_, alv) in skeleton.alveolae.iter() {
        for i in 1..(alv.len() >> 1) {
            writeln!(
                file,
                "f {}// {}// {}//",
                skel_ind_to_ind[&alv[alv.len() - i]],
                skel_ind_to_ind[&alv[i - 1]],
                skel_ind_to_ind[&alv[i]],
            )?;
            writeln!(
                file,
                "f {}// {}// {}//",
                skel_ind_to_ind[&alv[alv.len() - i - 1]],
                skel_ind_to_ind[&alv[alv.len() - i]],
                skel_ind_to_ind[&alv[i]],
            )?;
        }
        if alv.len() % 2 == 1 {
            let ind = alv.len() >> 1;
            writeln!(
                file,
                "f {}// {}// {}//",
                skel_ind_to_ind[&alv[ind - 1]],
                skel_ind_to_ind[&alv[ind]],
                skel_ind_to_ind[&alv[ind + 1]],
            )?;
        }
    }

    Ok(())
}
