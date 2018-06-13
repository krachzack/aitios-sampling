use super::uniform::Uniform;
use geom::Vec3;
use geom::prelude::*;
use rand;
use rand::distributions::{IndependentSample, Range};
use std::f32::consts::PI;

pub struct UnitSphere;

impl Uniform for UnitSphere {
    fn uniform(&self) -> Vec3 {
        let mut rng = rand::thread_rng();
        let between = Range::new(-1.0f32, 1.0);

        let dir = Vec3::new(
            between.ind_sample(&mut rng),
            between.ind_sample(&mut rng),
            between.ind_sample(&mut rng)
        );

        if !dir.is_zero() {
            dir.normalize()
        } else {
            // If, by pure chance got a zero vector, try again so we can normalize it
            self.uniform()
        }
    }
}

/// A hemisphere of radius 1, with the bottom disk aligned to a plane
pub enum UnitHemisphere {
    PosX, NegX,
    PosY, NegY,
    PosZ, NegZ
}

impl Uniform for UnitHemisphere {

    /// Uniformly samples the given hemisphere.
    ///
    /// Sampling method from this
    /// [blog article](http://www.rorydriscoll.com/2009/01/07/better-sampling/).
    fn uniform(&self) -> Vec3 {
        let u1 : f32 = rand::random();
        let u2 : f32 = rand::random();

        let r = u1.sqrt();
        let theta = 2.0 * PI * u2;

        let x = r * theta.cos();
        let y = r * theta.sin();
        let z = ((1.0 - u1).max(0.0)).sqrt();

        match self {
            &UnitHemisphere::PosZ => Vec3::new(x, y, z),
            &UnitHemisphere::NegZ => Vec3::new(x, y, -z),
            &UnitHemisphere::PosY => Vec3::new(x, z, y),
            &UnitHemisphere::NegY => Vec3::new(x, -z, y),
            &UnitHemisphere::PosX => Vec3::new(z, y, x),
            &UnitHemisphere::NegX => Vec3::new(-z, y, x)
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_normalized_unit_sphere() {
        (1..100).for_each(|_| assert_ulps_eq!(1.0, UnitSphere.uniform().magnitude()));
    }

    #[test]
    fn test_normalized_unit_hemisphere() {
        (1..100).for_each(|_| assert_ulps_eq!(1.0, UnitHemisphere::PosX.uniform().magnitude()));
        (1..100).for_each(|_| assert_ulps_eq!(1.0, UnitHemisphere::NegX.uniform().magnitude()));
        (1..100).for_each(|_| assert_ulps_eq!(1.0, UnitHemisphere::PosY.uniform().magnitude()));
        (1..100).for_each(|_| assert_ulps_eq!(1.0, UnitHemisphere::NegY.uniform().magnitude()));
        (1..100).for_each(|_| assert_ulps_eq!(1.0, UnitHemisphere::PosZ.uniform().magnitude()));
        (1..100).for_each(|_| assert_ulps_eq!(1.0, UnitHemisphere::NegZ.uniform().magnitude()));
    }

    fn uniform_naive(hemisphere: &UnitHemisphere) -> Vec3 {
        let (range_x, range_y, range_z) = match hemisphere {
            &UnitHemisphere::PosX => (Range::new(0.0f32, 1.0), Range::new(-1.0f32, 1.0), Range::new(-1.0f32, 1.0)),
            &UnitHemisphere::NegX => (Range::new(-1.0f32, 0.0), Range::new(-1.0f32, 1.0), Range::new(-1.0f32, 1.0)),
            &UnitHemisphere::PosY => (Range::new(-1.0f32, 1.0), Range::new(0.0f32, 1.0), Range::new(-1.0f32, 1.0)),
            &UnitHemisphere::NegY => (Range::new(-1.0f32, 1.0), Range::new(-1.0f32, 0.0), Range::new(-1.0f32, 1.0)),
            &UnitHemisphere::PosZ => (Range::new(-1.0f32, 1.0), Range::new(-1.0f32, 1.0), Range::new(0.0f32, 1.0)),
            &UnitHemisphere::NegZ => (Range::new(-1.0f32, 1.0), Range::new(-1.0f32, 1.0), Range::new(-1.0f32, 0.0)),
        };

        let mut rng = rand::thread_rng();
        let dir = Vec3::new(
            range_x.ind_sample(&mut rng),
            range_y.ind_sample(&mut rng),
            range_z.ind_sample(&mut rng)
        );

        if !dir.is_zero() {
            dir.normalize()
        } else {
            // If, by pure chance got a zero vector, try again so we can normalize it
            uniform_naive(hemisphere)
        }
    }
}
