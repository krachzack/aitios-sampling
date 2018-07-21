use super::Uniform;
use geom::{Position, Triangle, Vec3};
use rand;

impl<T: Triangle> Uniform for T {
    fn uniform(&self) -> Vec3 {
        let (a, b, c) = self.positions();
        let bary = sample_bary();

        bary[0] * a.position() + bary[1] * b.position() + bary[2] * c.position()
    }
}

pub fn sample_bary() -> [f32; 3] {
    let u = rand::random::<f32>();
    let v = rand::random::<f32>();

    let sqrt_u = u.sqrt();

    [1.0 - sqrt_u, (sqrt_u * (1.0 - v)), (sqrt_u * v)]
}

#[cfg(test)]
mod test {
    use super::*;
    use geom::{FromVertices, TupleTriangle};

    #[test]
    fn test_sample_tri_point() {
        let tri = TupleTriangle::new(
            Vec3::new(100.0, 100.0, 100.0),
            Vec3::new(200.0, 100.0, 100.0),
            Vec3::new(100.0, 200.0, 100.0),
        );

        for _ in 0..100 {
            let on_there = tri.uniform();
            assert_ulps_eq!(on_there.z, 100.0);
            assert!(on_there.x >= 100.0 && on_there.x < 200.0);
            assert!(on_there.y >= 100.0 && on_there.y < 200.0);
        }
    }
}
