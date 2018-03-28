use geom::Vec3;

/// Samples any point on the surface of the implementer.
/// All possible points have approximately equal probability.
pub trait Uniform {
    fn uniform(&self) -> Vec3;
}
