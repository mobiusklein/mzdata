const PROTON: f64 = 1.00727646677;

#[inline]
pub fn mass_charge_ratio(mass: f64, z: i32) -> f64 {
    (mass + z as f64 * PROTON) / (z.abs() as f64)
}

#[inline]
pub fn neutral_mass(mz: f64, z: i32) -> f64 {
    (mz * z.abs() as f64) - z as f64 * PROTON
}
