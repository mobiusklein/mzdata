
pub enum MassErrorType {
    Exact,
    PPM
}

impl MassErrorType {
    pub fn lower_bound(&self, query: f64, alt: f64) -> f64 {
        match self {
            Self::Exact => query - alt,
            Self::PPM => (query - alt) / alt
        }
    }

    pub fn call(&self, query: f64, alt: f64) -> f64 {
        match self {
            Self::Exact => query - alt,
            Self::PPM => (query - alt) / alt
        }
    }

    pub fn upper_bound(&self, query: f64, alt: f64) -> f64 {
        return self.lower_bound(query, -alt)
    }
}