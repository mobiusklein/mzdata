pub mod coordinate;
pub mod peak;
pub mod peak_set;

pub use crate::peaks::peak::{CentroidPeak, DeconvolutedPeak};
pub use crate::peaks::coordinate::{CoordinateDimension, CoordinateLike, IndexedCoordinate};
pub use crate::peaks::peak_set::{PeakCollection, PeakSet};