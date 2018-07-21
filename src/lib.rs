//! Provides sampling functionality, for instance:
//! * Uniformly sampling a point on a `Triangle`, the [`UnitSphere`](struct.UnitSphere.html) or a [`UnitHemisphere`](enum.UnitHemisphere.html) with the [`Uniform`](trait.Uniform.html) trait,
//! * efficiently selecting a triangle out of a large set by area with  [`TriangleBins`](struct.TriangleBins.html),
//! * sequences of uniformly sampled points, e.g. [`Poisson`](sequence/struct.Poisson.html) for a poisson disk set obtained from triangles.

#[cfg_attr(test, macro_use)]
extern crate aitios_geom as geom;
extern crate float_extras;
extern crate kdtree;
extern crate rand;

mod cosine_weighted;
pub mod sequence;
mod tri;
mod triangle_bins;
mod uniform;
mod unit;

pub use self::cosine_weighted::CosineWeighted;
pub use self::sequence::{into_poisson_disk_set, poisson_disk_set, Poisson};
pub use self::tri::sample_bary;
pub use self::triangle_bins::TriangleBins;
pub use self::uniform::Uniform;
pub use self::unit::*;
