use std::collections::HashMap;

use crate::skeleton3d::Skeleton3D;

/// Counts number of connected components in a skeleton
pub fn connected_components_count(skel: &Skeleton3D) -> usize {
    let mut node_cc: HashMap<usize, usize> = HashMap::new();

    for (&i, _) in skel.get_nodes().iter() {
        node_cc.insert(i, i);
    }

    for (_, [indv1, indv2]) in skel.get_lone_edges().iter() {
        let &cc1 = node_cc.get(indv1).unwrap();
        let &cc2 = node_cc.get(indv2).unwrap();

        let (cc1, cc2) = if cc1 < cc2 { (cc1, cc2) } else { (cc2, cc1) };
        node_cc.iter_mut().for_each(|(_, cc)| {
            if *cc == cc2 {
                *cc = cc1
            }
        });
    }

    for (_, [indv1, indv2]) in skel.get_edges_on_alv().iter() {
        let &cc1 = node_cc.get(indv1).unwrap();
        let &cc2 = node_cc.get(indv2).unwrap();

        if cc1 == cc2 {
            continue;
        }

        let (cc1, cc2) = if cc1 < cc2 { (cc1, cc2) } else { (cc2, cc1) };
        node_cc.iter_mut().for_each(|(_, cc)| {
            if *cc == cc2 {
                *cc = cc1
            }
        });
    }

    let mut cc_id: Vec<usize> = node_cc.iter().map(|(_, &cc)| cc).collect();
    cc_id.sort();
    cc_id.dedup();

    cc_id.len()
}
