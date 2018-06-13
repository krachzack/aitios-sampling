//! Provides sampling functionality, for instance:
//! * Uniformly sampling a point on a `Triangle`, the [`UnitSphere`](struct.UnitSphere.html) or a [`UnitHemisphere`](enum.UnitHemisphere.html) with the [`Uniform`](trait.Uniform.html) trait,
//! * efficiently selecting a triangle out of a large set by area with  [`TriangleBins`](struct.TriangleBins.html),
//! * sequences of uniformly sampled points, e.g. [`Poisson`](sequence/struct.Poisson.html) for a poisson disk set obtained from triangles.

#![feature(test)]
#[cfg(test)]
extern crate test;

#[cfg_attr(test, macro_use)]
extern crate aitios_geom as geom;
extern crate rand;
extern crate float_extras;
extern crate kdtree;

mod cosine_weighted;
mod uniform;
mod unit;
mod tri;
mod triangle_bins;
pub mod sequence;

pub use self::cosine_weighted::CosineWeighted;
pub use self::uniform::Uniform;
pub use self::unit::*;
pub use self::tri::sample_bary;
pub use self::triangle_bins::TriangleBins;
pub use self::sequence::{Poisson, poisson_disk_set, into_poisson_disk_set};
