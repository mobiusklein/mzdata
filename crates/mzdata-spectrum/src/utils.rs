use std::cmp::Ordering;

/// The different kinds of orientations ion mobility data may be present in.
#[derive(Debug, Default, Clone, Copy, Hash, PartialEq, Eq)]
pub enum HasIonMobility {
    /// No ion mobility measurement found
    #[default]
    None = 0,
    /// A single ion mobility point measurement
    Point = 1,
    /// Multiple ion mobility point measurements along an axis
    Dimension = 2
}

impl PartialOrd for HasIonMobility {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for HasIonMobility {
    fn cmp(&self, other: &Self) -> Ordering {
        match self {
            HasIonMobility::None => {
                match other {
                    HasIonMobility::None => Ordering::Equal,
                    _ => Ordering::Less
                }
            },
            HasIonMobility::Point => {
                match other {
                    HasIonMobility::None => Ordering::Greater,
                    HasIonMobility::Point => Ordering::Equal,
                    HasIonMobility::Dimension => Ordering::Less,
                }
            },
            HasIonMobility::Dimension => {
                match other {
                    HasIonMobility::Dimension => Ordering::Equal,
                    _ => Ordering::Greater
                }
            },
        }
    }
}
