use geom::Vec3;

/// Performs cosine weighted sampling, e.g. for monte carlo estimators.
pub trait CosineWeighted {
    fn cosine_weighted(&self) -> Vec3;
}
