use super::uniform::Uniform;
use geom::Vec3;
use geom::prelude::*;
use rand;
use rand::distributions::{IndependentSample, Range};

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
    fn uniform(&self) -> Vec3 {
        let (range_x, range_y, range_z) = match self {
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
            self.uniform()
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
}
