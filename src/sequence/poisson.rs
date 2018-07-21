use tri::sample_bary;
use triangle_bins::TriangleBins;

use geom::prelude::*;
use geom::tri::split_at_edge_midpoints;
use geom::{InterpolateVertex, Position, TangentSpace, Triangle, Vec3};

use kdtree::distance::squared_euclidean;
use kdtree::ErrorKind;
use kdtree::KdTree;

use std::f32::consts::PI;
use std::f32::EPSILON;

pub fn poisson_disk_set<'a, T, V, I>(triangles: I, min_point_distance: f32) -> Poisson<T>
where
    T: 'a + Clone + InterpolateVertex<Vertex = V> + FromVertices<Vertex = V>,
    V: Clone + Position,
    I: IntoIterator<Item = &'a T>,
{
    Poisson::new(triangles.into_iter().cloned(), min_point_distance)
}

pub fn into_poisson_disk_set<T, V, I>(triangles: I, min_point_distance: f32) -> Poisson<T>
where
    T: Clone + InterpolateVertex<Vertex = V> + FromVertices<Vertex = V>,
    V: Clone + Position,
    I: IntoIterator<Item = T>,
{
    Poisson::new(triangles.into_iter(), min_point_distance)
}

/// Represents a poisson disk set over triangles.
pub struct Poisson<T> {
    active_triangles: TriangleBins<T>,
    min_point_distance: f32,
    /// Do not split a triangle if the resulting subfragments would have a smaller area than this
    disregard_area: f32,
    /// Remembers samples already generated as the f64 array and also remembers the associated normal
    previous_samples: KdTree<Sample, [f64; 3]>,
}

#[derive(Debug)]
struct Sample {
    position: Vec3,
    normal: Vec3,
}

impl<T, V> Poisson<T>
where
    T: InterpolateVertex<Vertex = V> + FromVertices<Vertex = V>,
    V: Clone + Position,
{
    fn new<I>(triangles: I, min_point_distance: f32) -> Self
    where
        I: IntoIterator<Item = T>,
    {
        let active_triangles: TriangleBins<T> = triangles.into_iter().collect();

        // disregardiness is a factor for running time vs amount of generated points
        // A value closer to zero will result in a more even spacing of points by
        // generating more points and making the poisson disk set more maximal.
        // A higher value will decrease the number of performed triangle splits
        // and thus improve running time. The set will be less evenly spaced though.
        let disregardiness = 0.3;
        let disregard_area = disregardiness * min_point_distance * min_point_distance * PI;

        let previous_samples = KdTree::new(3);

        Poisson {
            active_triangles,
            min_point_distance,
            disregard_area,
            previous_samples,
        }
    }

    fn vtx_to_arr(vtx: &T::Vertex) -> [f64; 3] {
        let Vec3 { x, y, z } = vtx.position();
        [x as f64, y as f64, z as f64]
    }

    fn meets_minimum_distance_requirement(&self, vtx: &T::Vertex, vtx_normal: Vec3) -> bool {
        let position = Self::vtx_to_arr(vtx);

        // TODO do I even need to square this?
        //      API is undocumented, could come up with own kdtree
        let dist_sqr = (self.min_point_distance * self.min_point_distance) as f64;

        let within = self
            .previous_samples
            .within(&position, dist_sqr, &squared_euclidean);

        match within {
            // If no points found or all of the points within the distance have a normal rotated by
            // more than 90° with respect to the proposed new point, they are considered to be in proximity
            // and the test fails.
            // Otherwise the found points could be a close but unrelated surface on the other
            // side of a thin object
            Ok(within_points) => within_points.iter().all(|&(_, &Sample { normal, .. })| {
                !Self::holds_angle_requirement(normal, vtx_normal)
            }),
            Err(ErrorKind::ZeroCapacity) => true,
            Err(ErrorKind::NonFiniteCoordinate) => panic!("Vertex has infinite or NaN coordinate"),
            Err(other_error) => panic!("{:?}", other_error),
        }
    }

    fn add_sample(&mut self, vtx: &T::Vertex, normal: Vec3) {
        let position = Self::vtx_to_arr(vtx);
        self.previous_samples
            .add(
                position,
                Sample {
                    position: vtx.position(),
                    normal,
                },
            )
            .expect("Failed to remember a sample for later in kdtree");
    }

    /// Samples a vertex by interpolating a new vertex at random barycentric coordinates
    fn generate_sample(fragment: &T) -> T::Vertex {
        fragment.interpolate_vertex_bary(sample_bary())
    }

    fn is_covered(&self, fragment: &T) -> bool {
        let (proposed_covering_point, min_cover_radius_sqr) = minimum_bounding_sphere_sqr(fragment);
        let r = 0.5 * self.min_point_distance;

        if min_cover_radius_sqr > (r * r) {
            false
        } else {
            let proposed_covering_point = [
                proposed_covering_point.x as f64,
                proposed_covering_point.y as f64,
                proposed_covering_point.z as f64,
            ];

            let within = self.previous_samples.within(
                &proposed_covering_point,
                min_cover_radius_sqr as f64,
                &squared_euclidean,
            );

            let fragment_normal = fragment.normal();

            match within {
                Ok(within_points) => within_points.iter().any(|&(_, sample)| {
                    let &Sample { position, normal } = sample;

                    // If the normal of an existing sample differs by more than 90° compared
                    // to the normal of the fragment, do not take it into account in the covered
                    // check. This ensures surfel generation will not be influenced by generation
                    // of surfels on the other side of a thin surface with respect to the min
                    // point distance.
                    Self::holds_angle_requirement(normal, fragment_normal)
                        && fragment.is_inside_sphere(
                            Vec3::new(position[0] as f32, position[1] as f32, position[2] as f32),
                            r,
                        )
                }),
                Err(ErrorKind::ZeroCapacity) => false,
                Err(other_error) => panic!("{:?}", other_error),
            }
        }
    }

    /// Checks if two normals have an interior angle of less than 90°
    fn holds_angle_requirement(normal1: Vec3, normal2: Vec3) -> bool {
        const DISREGARD_ANGLE_COS: f32 = 0.0;
        normal1.dot(normal2) > DISREGARD_ANGLE_COS
    }
}

impl<T, V> Iterator for Poisson<T>
where
    T: InterpolateVertex<Vertex = V> + FromVertices<Vertex = V>,
    V: Clone + Position,
{
    type Item = T::Vertex;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.active_triangles.triangle_count() == 0 {
                return None;
            }

            let fragment = self.active_triangles.pop();

            let sample = {
                let sample_candidate = Self::generate_sample(&fragment);
                let fragment_normal = fragment.normal();

                if self.meets_minimum_distance_requirement(&sample_candidate, fragment_normal) {
                    self.add_sample(&sample_candidate, fragment_normal);
                    Some(sample_candidate)
                } else {
                    None
                }
            };

            let split_area = 0.25 * fragment.area();
            if split_area > self.disregard_area {
                if !self.is_covered(&fragment) {
                    for sub_fragment in split_at_edge_midpoints(&fragment) {
                        if !self.is_covered(&sub_fragment) {
                            self.active_triangles.push(sub_fragment)
                        }
                    }
                }
            }

            if let Some(_) = sample {
                return sample;
            }
        }
    }
}

fn minimum_bounding_sphere_sqr<T: Triangle>(tri: &T) -> (Vec3, f32) {
    let (a, b, c) = tri.positions();
    let dot_abab = (b - a).dot(b - a);
    let dot_abac = (b - a).dot(c - a);
    let dot_acac = (c - a).dot(c - a);
    let d = 2.0 * (dot_abab * dot_acac - dot_abac * dot_abac);
    let mut reference_point = a;

    let center = if d.abs() <= EPSILON {
        // a, b, and c lie on a line. Circle center is center of AABB of the
        // points, and radius is distance from circle center to AABB corner
        let bbox = tri.bounds();
        reference_point = bbox.min;
        0.5 * (bbox.min + bbox.max)
    } else {
        let s = (dot_abab * dot_acac - dot_acac * dot_abac) / d;
        let t = (dot_acac * dot_abab - dot_abab * dot_abac) / d;
        // s controls height over AC, t over AB, (1-s-t) over BC
        if s <= 0.0 {
            0.5 * (a + c)
        } else if t <= 0.0 {
            0.5 * (a + b)
        } else if (s + t) >= 1.0 {
            reference_point = b;
            0.5 * (b + c)
        } else {
            a + s * (b - a) + t * (c - a)
        }
    };

    let radius_sqr = center.distance2(reference_point);

    (center, radius_sqr)
}

#[cfg(test)]
mod test {
    use super::*;
    use geom::{FromVertices, TupleTriangle};

    #[test]
    fn test_sample_quad_poisson_disk_set() {
        let quad = [
            TupleTriangle::new(
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(1.0, 1.0, 0.0),
            ),
            TupleTriangle::new(
                Vec3::new(1.0, 1.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(0.0, 0.0, 0.0),
            ),
        ];

        let minimum_distance = 0.05;
        let vertices: Vec<_> = poisson_disk_set(quad.iter(), minimum_distance).collect();

        assert!(
            vertices.len() > 50,
            "Unexpectedly low number of points yielded"
        );

        for vtx0 in &vertices {
            let mut below_min_distance_cnt = 0;
            let mut above_min_distance_cnt = 0;

            for vtx1 in &vertices {
                if vtx1.position().distance2(*vtx0) >= (minimum_distance * minimum_distance) {
                    above_min_distance_cnt += 1;
                } else {
                    below_min_distance_cnt += 1;
                }
            }

            assert_eq!(
                1, below_min_distance_cnt,
                "Only the same vertex should have distance < min distance"
            );
            assert_eq!(
                vertices.len() - 1,
                above_min_distance_cnt,
                "All vertices except the same should have distance > min distance"
            );
        }
    }
}
