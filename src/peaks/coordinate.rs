
#[derive(Default, Debug, Clone, Copy)]
pub struct MZ {}
#[derive(Default, Debug, Clone, Copy)]
pub struct Mass {}
#[derive(Default, Debug, Clone, Copy)]
pub struct Time {}
#[derive(Default, Debug, Clone, Copy)]
pub struct IonMobility {}

pub enum CoordinateDimension {
    MZ(MZ),
    Mass(Mass),
    Time(Time),
    IonMobility(IonMobility),
}

pub trait CoordinateLike<T> : PartialOrd {
    fn get_coordinate(&self) -> f64;
}

pub trait IndexedCoordinate<T> : CoordinateLike<T> {
    fn get_index(&self) -> u32;
    fn set_index(&mut self, index: u32);
}
