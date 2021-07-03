pub struct MZ {}
pub struct Mass {}
pub struct Time {}

pub enum CoordinateDimension {
    MZ(MZ),
    Mass(Mass),
    Time(Time),
}

pub trait CoordinateLike<T> {
    fn get_coordinate(&self) -> f64;
}