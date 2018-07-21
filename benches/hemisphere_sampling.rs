#![feature(test)]

extern crate aitios_geom;
extern crate aitios_sampling;
extern crate rand;
extern crate test;

use aitios_geom::{InnerSpace, Vec3, Zero};
use aitios_sampling::{CosineWeighted, Uniform, UnitHemisphere};
use rand::distributions::{IndependentSample, Range};

#[bench]
fn bench_uniform_hemisphere_sampling(bencher: &mut test::Bencher) {
    bencher.iter(|| UnitHemisphere::PosZ.uniform());
}

#[bench]
fn bench_cosine_weighted_hemisphere_sampling(bencher: &mut test::Bencher) {
    bencher.iter(|| UnitHemisphere::PosZ.cosine_weighted());
}

#[bench]
fn bench_naive_hemisphere_sampling(bencher: &mut test::Bencher) {
    bencher.iter(|| uniform_naive(&UnitHemisphere::PosZ));
}

/// Naive algorithm for performance comparison
fn uniform_naive(hemisphere: &UnitHemisphere) -> Vec3 {
    let (range_x, range_y, range_z) = match hemisphere {
        &UnitHemisphere::PosX => (
            Range::new(0.0f32, 1.0),
            Range::new(-1.0f32, 1.0),
            Range::new(-1.0f32, 1.0),
        ),
        &UnitHemisphere::NegX => (
            Range::new(-1.0f32, 0.0),
            Range::new(-1.0f32, 1.0),
            Range::new(-1.0f32, 1.0),
        ),
        &UnitHemisphere::PosY => (
            Range::new(-1.0f32, 1.0),
            Range::new(0.0f32, 1.0),
            Range::new(-1.0f32, 1.0),
        ),
        &UnitHemisphere::NegY => (
            Range::new(-1.0f32, 1.0),
            Range::new(-1.0f32, 0.0),
            Range::new(-1.0f32, 1.0),
        ),
        &UnitHemisphere::PosZ => (
            Range::new(-1.0f32, 1.0),
            Range::new(-1.0f32, 1.0),
            Range::new(0.0f32, 1.0),
        ),
        &UnitHemisphere::NegZ => (
            Range::new(-1.0f32, 1.0),
            Range::new(-1.0f32, 1.0),
            Range::new(-1.0f32, 0.0),
        ),
    };

    let mut rng = rand::thread_rng();
    let dir = Vec3::new(
        range_x.ind_sample(&mut rng),
        range_y.ind_sample(&mut rng),
        range_z.ind_sample(&mut rng),
    );

    if !dir.is_zero() {
        dir.normalize()
    } else {
        // If, by pure chance got a zero vector, try again so we can normalize it
        uniform_naive(hemisphere)
    }
}
