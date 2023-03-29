use anyhow::Result;

use super::SkeletonInterface3D;
use super::SkeletonSingularPath;

/// Sepration on skeleton
pub struct SkeletonSeparation<'a, 'b> {
    skeleton_interface: &'b mut SkeletonInterface3D<'a>,
    external_path: SkeletonSingularPath,
    internal_paths: Vec<SkeletonSingularPath>,
}

impl<'a, 'b> SkeletonSeparation<'a, 'b> {
    /// Separation constructor
    pub fn create(
        skeleton_interface: &'b mut SkeletonInterface3D<'a>,
        ind_pedge: usize,
    ) -> Result<SkeletonSeparation<'a, 'b>> {
        let external_path = SkeletonSingularPath::create(ind_pedge, skeleton_interface);
        Ok(SkeletonSeparation {
            skeleton_interface,
            external_path,
            internal_paths: Vec::new(),
        })
    }

    pub fn from_singular_path(
        skeleton_interface: &'b mut SkeletonInterface3D<'a>,
        sing_path: SkeletonSingularPath,
    ) -> SkeletonSeparation<'a, 'b> {
        SkeletonSeparation {
            skeleton_interface,
            external_path: sing_path,
            internal_paths: Vec::new(),
        }
    }

    /// Skeleton interface getter
    pub fn skeleton_interface(&self) -> &SkeletonInterface3D {
        self.skeleton_interface
    }

    /// External path getter
    pub fn external_path(&self) -> &SkeletonSingularPath {
        &self.external_path
    }

    /// Internal path getter
    pub fn internal_paths(&self) -> &Vec<SkeletonSingularPath> {
        &self.internal_paths
    }

    fn follow_external_path(&mut self) -> Result<()> {
        self.external_path
            .follow_singular_path(&mut self.skeleton_interface)
    }

    fn internal_partial_edges(&self) -> Vec<usize> {
        let vec_pedges_ext = self.external_path.ind_partial_edges();
        let mut vec_pedges_int = Vec::new();
        for &ind_pedge in vec_pedges_ext.iter() {
            let ind_pedge_opp = self
                .skeleton_interface
                .get_partial_edge_uncheck(ind_pedge)
                .partial_edge_opposite()
                .ind();
            if vec_pedges_ext
                .iter()
                .find(|&&ind| ind == ind_pedge_opp)
                .is_none()
            {
                vec_pedges_int.push(ind_pedge_opp);
            }
        }
        vec_pedges_int.sort();
        vec_pedges_int.dedup();
        vec_pedges_int
    }

    fn follow_internal_paths(&mut self) -> Result<()> {
        let mut vec_internal_pedges = self.internal_partial_edges();
        loop {
            if let Some(ind_pedge) = vec_internal_pedges.pop() {
                let mut skeleton_path_int =
                    SkeletonSingularPath::create(ind_pedge, &mut self.skeleton_interface);
                skeleton_path_int.follow_singular_path(&mut self.skeleton_interface)?;
                for &ind_pedge_new in skeleton_path_int.ind_partial_edges().iter() {
                    if let Some(pos) = vec_internal_pedges
                        .iter()
                        .position(|&ind| ind == ind_pedge_new)
                    {
                        vec_internal_pedges.remove(pos);
                    }
                }
                self.internal_paths.push(skeleton_path_int);
            } else {
                break;
            }
        }
        Ok(())
    }

    /// Computes external path and internal paths associate to separation
    pub fn follow_separation(&mut self) -> Result<()> {
        self.follow_external_path()?;
        self.follow_internal_paths()?;
        Ok(())
    }

    /// Checks if path is closable
    pub fn closable_path(&self) -> Result<bool> {
        self.external_path.closable_path(&self.skeleton_interface)
        // Ok(true)
    }
}
