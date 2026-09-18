use std::{collections::HashSet, convert::TryInto, mem};

use mzpeaks::{
    CentroidLike, CentroidPeak, CoordinateLike, DeconvolutedPeak, DeconvolutedPeakSet, IonMobility,
    MZ, MZPeakSetType, Mass, PeakSet,
    feature::{ChargedFeature, Feature},
    peak::{IonMobilityAwareCentroidPeak, IonMobilityAwareDeconvolutedPeak},
    prelude::*,
};

use mzdata_param::Unit;

use crate::utils::{mass_charge_ratio, neutral_mass};

use super::{
    BinaryArrayMap3D, ByteArrayView,
    array::DataArray,
    encodings::{ArrayRetrievalError, ArrayType, BinaryDataArrayType},
    map::BinaryArrayMap,
};

impl From<&PeakSet> for BinaryArrayMap {
    fn from(peaks: &PeakSet) -> BinaryArrayMap {
        let mut arrays = BinaryArrayMap::new();

        let mut mz_array = DataArray::from_name_type_size(
            &ArrayType::MZArray,
            BinaryDataArrayType::Float64,
            peaks.len() * BinaryDataArrayType::Float64.size_of(),
        );

        let mut intensity_array = DataArray::from_name_type_size(
            &ArrayType::IntensityArray,
            BinaryDataArrayType::Float32,
            peaks.len() * BinaryDataArrayType::Float32.size_of(),
        );

        for p in peaks.iter() {
            let mz: f64 = p.coordinate();
            let inten: f32 = p.intensity();

            let raw_bytes: [u8; mem::size_of::<f64>()] = mz.to_le_bytes();
            mz_array.data.extend(raw_bytes);

            let raw_bytes: [u8; mem::size_of::<f32>()] = inten.to_le_bytes();
            intensity_array.data.extend(raw_bytes);
        }

        arrays.add(mz_array);
        arrays.add(intensity_array);
        arrays
    }
}

impl<C: CentroidLike + From<CentroidPeak>> From<BinaryArrayMap> for MZPeakSetType<C> {
    fn from(arrays: BinaryArrayMap) -> MZPeakSetType<C> {
        (&arrays).into()
    }
}

impl<C: CentroidLike + From<CentroidPeak>> From<&BinaryArrayMap> for MZPeakSetType<C> {
    fn from(arrays: &BinaryArrayMap) -> MZPeakSetType<C> {
        let mz_array = arrays.mzs().unwrap();
        let intensity_array = arrays.intensities().unwrap();
        let mut peaks = Vec::with_capacity(mz_array.len());

        for (i, (mz, intensity)) in mz_array.iter().zip(intensity_array.iter()).enumerate() {
            peaks.push(
                CentroidPeak {
                    mz: *mz,
                    intensity: *intensity,
                    index: i as u32,
                }
                .into(),
            )
        }
        MZPeakSetType::<C>::new(peaks)
    }
}

impl From<&BinaryArrayMap> for DeconvolutedPeakSet {
    fn from(arrays: &BinaryArrayMap) -> DeconvolutedPeakSet {
        let mz_array = arrays.mzs().unwrap();
        let intensity_array = arrays.intensities().unwrap();
        let charge_array = arrays
            .charges()
            .expect("Charge state array is required for deconvoluted peaks");
        let mut peaks = Vec::with_capacity(mz_array.len());
        for (i, ((mz, intensity), charge)) in mz_array
            .iter()
            .zip(intensity_array.iter())
            .zip(charge_array.iter())
            .enumerate()
        {
            peaks.push(DeconvolutedPeak {
                neutral_mass: neutral_mass(*mz, *charge),
                intensity: *intensity,
                charge: *charge,
                index: i as u32,
            })
        }

        DeconvolutedPeakSet::new(peaks)
    }
}

impl From<&DeconvolutedPeakSet> for BinaryArrayMap {
    fn from(peaks: &DeconvolutedPeakSet) -> BinaryArrayMap {
        let mut arrays = BinaryArrayMap::new();

        let mut mz_array = DataArray::from_name_type_size(
            &ArrayType::MZArray,
            BinaryDataArrayType::Float64,
            peaks.len() * BinaryDataArrayType::Float64.size_of(),
        );

        let mut intensity_array = DataArray::from_name_type_size(
            &ArrayType::IntensityArray,
            BinaryDataArrayType::Float32,
            peaks.len() * BinaryDataArrayType::Float32.size_of(),
        );

        let mut charge_array = DataArray::from_name_type_size(
            &ArrayType::ChargeArray,
            BinaryDataArrayType::Int32,
            peaks.len() * BinaryDataArrayType::Int32.size_of(),
        );

        for p in peaks.iter() {
            let mz: f64 = p.mz();
            let inten: f32 = p.intensity();
            let charge = p.charge();

            let raw_bytes: [u8; mem::size_of::<f64>()] = mz.to_le_bytes();
            mz_array.data.extend_from_slice(&raw_bytes);

            let raw_bytes: [u8; mem::size_of::<f32>()] = inten.to_le_bytes();
            intensity_array.data.extend_from_slice(&raw_bytes);

            let raw_bytes: [u8; mem::size_of::<i32>()] = charge.to_le_bytes();
            charge_array.data.extend_from_slice(&raw_bytes);
        }

        arrays.add(mz_array);
        arrays.add(intensity_array);
        arrays.add(charge_array);
        arrays
    }
}

/// Whether or not sufficient arrays are available in a [`BinaryArrayMap`] to
/// construct a peak list using [`BuildFromArrayMap`] or [`BuildFromArrayMap3D`]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ArraysAvailable {
    /// The required arrays were not listed, so cannot say yes or no
    Unknown,
    /// All the required arrays are present
    Ok,
    /// The list off arrays that are missing
    MissingArrays(Vec<ArrayType>),
}

/// A trait for reconstructing a `Vec` of peaks or features from a [`BinaryArrayMap`].
///
/// This is the inverse of [`BuildArrayMapFrom`].
///
/// # Examples
/// ```
/// use mzpeaks::CentroidPeak;
/// use mzdata_bindata::{ArraysAvailable, BinaryArrayMap, BuildFromArrayMap};
///
/// let empty = BinaryArrayMap::new();
/// assert!(matches!(
///     CentroidPeak::has_arrays_for(&empty),
///     ArraysAvailable::MissingArrays(_)
/// ));
/// ```
pub trait BuildFromArrayMap: Sized {
    /// The arrays that are *required* to reconstruct this type
    ///
    /// The default implementation returns `None`.
    fn arrays_required() -> Option<Vec<ArrayType>> {
        None
    }

    /// Try to build a [`Vec`] of [`Self`] from a [`BinaryArrayMap`].
    ///
    /// If arrays are absent or malformed, returns [`ArrayRetrievalError`]
    fn try_from_arrays(arrays: &BinaryArrayMap) -> Result<Vec<Self>, ArrayRetrievalError>;

    /// A shortcut form of [`Self::try_from_arrays`] that panics if it fails
    fn from_arrays(arrays: &BinaryArrayMap) -> Vec<Self> {
        Self::try_from_arrays(arrays).unwrap()
    }

    /// A pre-emptive check for the presence of the required arrays.
    fn has_arrays_for(arrays: &BinaryArrayMap) -> ArraysAvailable {
        if let Some(arrays_required) = Self::arrays_required() {
            let missing: Vec<_> = arrays_required
                .into_iter()
                .filter(|array_type| !arrays.has_array(array_type))
                .collect();
            if !missing.is_empty() {
                ArraysAvailable::MissingArrays(missing)
            } else {
                ArraysAvailable::Ok
            }
        } else {
            ArraysAvailable::Unknown
        }
    }
}

/// A trait for converting a slice of peaks or features into a [`BinaryArrayMap`].
///
/// This is the inverse of [`BuildFromArrayMap`].
///
/// # Examples
/// ```
/// use mzpeaks::CentroidPeak;
/// use mzdata_bindata::{ArrayType, BuildArrayMapFrom, BuildFromArrayMap};
///
/// let peaks = vec![
///     CentroidPeak::new(500.0, 1200.0, 0),
///     CentroidPeak::new(750.5, 800.0, 1),
/// ];
///
/// let arrays = CentroidPeak::as_arrays(&peaks);
/// assert!(arrays.has_array(&ArrayType::MZArray));
/// assert!(arrays.has_array(&ArrayType::IntensityArray));
///
/// let roundtripped = CentroidPeak::try_from_arrays(&arrays).unwrap();
/// assert_eq!(roundtripped.len(), 2);
/// assert_eq!(roundtripped[0].mz, 500.0);
/// ```
pub trait BuildArrayMapFrom: Sized {
    /// The [`ArrayType`] produced by this peak type when converted.
    ///
    /// The default implementation returns `None`.
    fn arrays_included(&self) -> Option<Vec<ArrayType>> {
        None
    }

    /// Serialize a slice of [`Self`] into a [`BinaryArrayMap`].
    ///
    /// This assumes that the operation cannot fail as all the data
    /// is available in memory already.
    fn as_arrays(source: &[Self]) -> BinaryArrayMap;
}

/// An extension of [`BuildArrayMapFrom`] that produces [`BinaryArrayMap3D`]
/// instead.
///
/// It is intended for use with ion mobility frames like [`MultiLayerIonMobilityFrame`](mzdata::MultiLayerIonMobilityFrame)
/// and feature types like [`Feature`](mzpeaks::feature::Feature).
pub trait BuildArrayMap3DFrom: BuildArrayMapFrom {
    /// Serialize a slice of [`Self`] into a [`BinaryArrayMap3D`].
    ///
    /// The default implementation just calls [`BuildArrayMapFrom::as_arrays`] and
    /// then reshapes them into a [`BinaryArrayMap3D`]. This leverages the requirement
    /// that a flat [`BinaryArrayMap`] can already be built from the same type, but may
    /// be more expensive than expected because it then makes copies of all of the data
    /// to re-arrange in a 3D layout.
    ///
    /// # Panics
    /// The default implementation panics if the intermediate [`BinaryArrayMap`] cannot
    /// actually be converted into a [`BinaryArrayMap3D`]
    fn as_arrays_3d(source: &[Self]) -> BinaryArrayMap3D {
        BuildArrayMapFrom::as_arrays(source).try_into().unwrap()
    }
}

/// An extension of [`BuildFromArrayMap`] that consumes [`BinaryArray3D`]
/// to build a [`Vec`] of [`Self`].
///
/// It is intended for use with ion mobility frames like [`MultiLayerIonMobilityFrame`](mzdata::MultiLayerIonMobilityFrame)
/// and feature types like [`Feature`](mzpeaks::feature::Feature).
pub trait BuildFromArrayMap3D: BuildFromArrayMap {
    /// Try to consume a [`BinaryArrayMap3D`] to build a `Vec` of [`Self`].
    ///
    /// The default implementation calls [`BinaryArrayMap3D::unstack`] and tries
    /// to use [`BuildFromArrayMap::try_from_arrays`]. A direct implementation
    /// would be much more efficient.
    fn try_from_arrays_3d(arrays: &BinaryArrayMap3D) -> Result<Vec<Self>, ArrayRetrievalError> {
        BuildFromArrayMap::try_from_arrays(&arrays.unstack()?)
    }

    /// A shortcut form of [`Self::try_from_arrays_3d`] that panics if it fails
    fn from_arrays_3d(arrays: &BinaryArrayMap3D) -> Vec<Self> {
        Self::try_from_arrays_3d(arrays).unwrap()
    }

    /// A pre-emptive check for the presence of the required arrays, particularly
    /// across all ion mobility points.
    fn has_arrays_3d_for(arrays: &BinaryArrayMap3D) -> ArraysAvailable {
        if let Some(arrays_required) = Self::arrays_required() {
            let arrays_required: Vec<_> = arrays_required
                .into_iter()
                .filter(|a| !a.is_ion_mobility())
                .collect();
            let mut arrays_not_seen: HashSet<_> = arrays_required.iter().cloned().collect();
            let mut missing = Vec::new();
            for (_, arr) in arrays.iter() {
                if arr.is_empty() {
                    continue;
                }
                missing.clear();
                for array_type in arrays_required.iter() {
                    if arr.has_array(array_type) {
                        arrays_not_seen.remove(array_type);
                    } else {
                        missing.push(array_type);
                    }
                }
                if missing.is_empty() {
                    return ArraysAvailable::Ok;
                }
            }
            if arrays_not_seen.is_empty() {
                ArraysAvailable::Unknown
            } else {
                ArraysAvailable::MissingArrays(arrays_not_seen.into_iter().collect())
            }
        } else {
            ArraysAvailable::Unknown
        }
    }
}

// Basic peaks

/// The basic [`CentroidPeak`] implements [`BuildArrayMapFrom`]
impl BuildArrayMapFrom for CentroidPeak {
    /// [`CentroidPeak`] produces [`ArrayType::MZArray`] and [`ArrayType::IntensityArray`]
    fn arrays_included(&self) -> Option<Vec<ArrayType>> {
        Some(vec![ArrayType::MZArray, ArrayType::IntensityArray])
    }

    fn as_arrays(source: &[Self]) -> BinaryArrayMap {
        let mut arrays = BinaryArrayMap::new();

        let mut mz_array = DataArray::from_name_type_size(
            &ArrayType::MZArray,
            BinaryDataArrayType::Float64,
            source.len() * BinaryDataArrayType::Float64.size_of(),
        );
        mz_array.unit = Unit::MZ;

        let mut intensity_array = DataArray::from_name_type_size(
            &ArrayType::IntensityArray,
            BinaryDataArrayType::Float32,
            source.len() * BinaryDataArrayType::Float32.size_of(),
        );
        intensity_array.unit = Unit::DetectorCounts;

        for p in source.iter() {
            let mz: f64 = p.coordinate();
            let inten: f32 = p.intensity();

            let raw_bytes: [u8; mem::size_of::<f64>()] = mz.to_le_bytes();
            mz_array.data.extend_from_slice(&raw_bytes);

            let raw_bytes: [u8; mem::size_of::<f32>()] = inten.to_le_bytes();
            intensity_array.data.extend_from_slice(&raw_bytes);
        }

        arrays.add(mz_array);
        arrays.add(intensity_array);
        arrays
    }
}

impl BuildFromArrayMap for CentroidPeak {
    fn try_from_arrays(arrays: &BinaryArrayMap) -> Result<Vec<Self>, ArrayRetrievalError> {
        let mz_array = arrays.mzs()?;
        let intensity_array = arrays.intensities()?;
        let mut peaks = Vec::with_capacity(mz_array.len());

        for (i, (mz, intensity)) in mz_array.iter().zip(intensity_array.iter()).enumerate() {
            peaks.push(CentroidPeak {
                mz: *mz,
                intensity: *intensity,
                index: i as u32,
            })
        }
        Ok(peaks)
    }

    /// [`CentroidPeak`] requires [`ArrayType::MZArray`] and [`ArrayType::IntensityArray`]
    fn arrays_required() -> Option<Vec<ArrayType>> {
        Some(vec![ArrayType::MZArray, ArrayType::IntensityArray])
    }
}

/// The basic [`DeconvolutedPeak`] implements [`BuildArrayMapFrom`].
///
/// It produces the much more common [`ArrayType::MZArray`] instead of
/// [`ArrayType::MassArray`] for compatibility with other tools, and assumes
/// that the charge carrier is a proton.
impl BuildArrayMapFrom for DeconvolutedPeak {
    fn as_arrays(source: &[Self]) -> BinaryArrayMap {
        let mut arrays = BinaryArrayMap::new();

        let mut mz_array = DataArray::from_name_type_size(
            &ArrayType::MZArray,
            BinaryDataArrayType::Float64,
            source.len() * BinaryDataArrayType::Float64.size_of(),
        );
        mz_array.unit = Unit::MZ;

        let mut intensity_array = DataArray::from_name_type_size(
            &ArrayType::IntensityArray,
            BinaryDataArrayType::Float32,
            source.len() * BinaryDataArrayType::Float32.size_of(),
        );
        intensity_array.unit = Unit::DetectorCounts;

        let mut charge_array = DataArray::from_name_type_size(
            &ArrayType::ChargeArray,
            BinaryDataArrayType::Int32,
            source.len() * BinaryDataArrayType::Int32.size_of(),
        );

        for p in source.iter() {
            let mz: f64 = p.mz();
            let inten: f32 = p.intensity();
            let charge = p.charge();

            let raw_bytes: [u8; mem::size_of::<f64>()] = mz.to_le_bytes();
            mz_array.data.extend_from_slice(&raw_bytes);

            let raw_bytes: [u8; mem::size_of::<f32>()] = inten.to_le_bytes();
            intensity_array.data.extend_from_slice(&raw_bytes);

            let raw_bytes: [u8; mem::size_of::<i32>()] = charge.to_le_bytes();
            charge_array.data.extend_from_slice(&raw_bytes);
        }

        arrays.add(mz_array);
        arrays.add(intensity_array);
        arrays.add(charge_array);
        arrays
    }

    /// [`CentroidPeak`] produces [`ArrayType::MZArray`], [`ArrayType::IntensityArray`]
    /// and [`ArrayType::ChargeArray`]
    fn arrays_included(&self) -> Option<Vec<ArrayType>> {
        Some(vec![
            ArrayType::MZArray,
            ArrayType::IntensityArray,
            ArrayType::ChargeArray,
        ])
    }
}

/// The basic [`DeconvolutedPeak`] implements [`BuildFromArrayMap`].
///
/// It assumes that the charge carrier was a proton when converting to neutral mass
/// from [`ArrayType::MZArray`].
impl BuildFromArrayMap for DeconvolutedPeak {
    fn try_from_arrays(arrays: &BinaryArrayMap) -> Result<Vec<Self>, ArrayRetrievalError> {
        let mz_array = arrays.mzs()?;
        let intensity_array = arrays.intensities()?;
        let charge_array = arrays.charges()?;
        let mut peaks = Vec::with_capacity(mz_array.len());
        for (i, ((mz, intensity), charge)) in mz_array
            .iter()
            .zip(intensity_array.iter())
            .zip(charge_array.iter())
            .enumerate()
        {
            peaks.push(DeconvolutedPeak {
                neutral_mass: neutral_mass(*mz, *charge),
                intensity: *intensity,
                charge: *charge,
                index: i as u32,
            })
        }

        Ok(peaks)
    }
}

// Ion mobility features

/// [`Feature`] with an [`IonMobility`] dimension implement
/// [`BuildArrayMapFrom`].
///
/// It uses a non-standard "feature identifier array" which assigns
/// a unique integer identifier to all points belonging to the same
/// [`Feature`].
///
/// It also produces [`ArrayType::RawIonMobilityArray`] which does not have
/// an associated [`Unit`]. The caller *should* update this if context is available.
impl BuildArrayMapFrom for Feature<MZ, IonMobility> {
    /// [`Feature`] produces [`ArrayType::MZArray`] and [`ArrayType::IntensityArray`],
    /// [`ArrayType::RawIonMobilityArray`], and a non-standard "feature identifier array"
    fn arrays_included(&self) -> Option<Vec<ArrayType>> {
        Some(vec![
            ArrayType::MZArray,
            ArrayType::IntensityArray,
            ArrayType::RawIonMobilityArray,
            ArrayType::nonstandard("feature identifier array"),
        ])
    }

    /// [`Feature`] uses a non-standard "feature identifier array" which assigns
    /// a unique integer identifier to all points belonging to the same
    /// [`Feature`].
    ///
    /// It also produces [`ArrayType::RawIonMobilityArray`] which does not have
    /// an associated [`Unit`]. The caller *should* update these for consistency
    fn as_arrays(source: &[Self]) -> BinaryArrayMap {
        let mut arrays = BinaryArrayMap::new();
        let n: usize = source.iter().map(|f| f.len()).sum();

        let mut mz_array = DataArray::from_name_type_size(
            &ArrayType::MZArray,
            BinaryDataArrayType::Float64,
            n * BinaryDataArrayType::Float64.size_of(),
        );
        mz_array.unit = Unit::MZ;

        let mut intensity_array = DataArray::from_name_type_size(
            &ArrayType::IntensityArray,
            BinaryDataArrayType::Float32,
            n * BinaryDataArrayType::Float32.size_of(),
        );
        intensity_array.unit = Unit::DetectorCounts;

        let mut ion_mobility_array = DataArray::from_name_type_size(
            &ArrayType::RawIonMobilityArray,
            BinaryDataArrayType::Float64,
            n * BinaryDataArrayType::Float64.size_of(),
        );

        let mut marker_array = DataArray::from_name_type_size(
            &ArrayType::nonstandard("feature identifier array"),
            BinaryDataArrayType::Int32,
            n * BinaryDataArrayType::Int32.size_of(),
        );

        let mut acc = Vec::with_capacity(n);
        source.iter().enumerate().for_each(|(i, f)| {
            f.iter()
                .for_each(|(mz, im, inten)| acc.push((mz, im, inten, i)))
        });
        acc.sort_by(|(mz_a, im_a, _, key_a), (mz_b, im_b, _, key_b)| {
            mz_a.total_cmp(mz_b)
                .then(im_a.total_cmp(im_b))
                .then(key_a.cmp(key_b))
        });

        for (mz, im, inten, key) in acc.iter() {
            mz_array.data.extend_from_slice(&mz.to_le_bytes());
            intensity_array.data.extend_from_slice(&inten.to_le_bytes());
            ion_mobility_array.data.extend_from_slice(&im.to_le_bytes());
            marker_array
                .data
                .extend_from_slice(&(*key as i32).to_le_bytes());
        }

        arrays.add(mz_array);
        arrays.add(intensity_array);
        arrays.add(ion_mobility_array);
        arrays.add(marker_array);
        arrays
    }
}

/// [`Feature`] with an [`IonMobility`] dimension implement
/// [`BuildFromArrayMap`].
///
/// It requires a non-standard "feature identifier array" which assigns
/// a unique integer identifier to all points belonging to the same
/// [`Feature`].
///
/// While it explicitly requires [`ArrayType::RawIonMobilityArray`], it will accept
/// any raw ion mobility array type.
impl BuildFromArrayMap for Feature<MZ, IonMobility> {
    /// [`Feature`] requires [`ArrayType::MZArray`] and [`ArrayType::IntensityArray`],
    /// [`ArrayType::RawIonMobilityArray`], and a non-standard "feature identifier array"
    fn arrays_required() -> Option<Vec<ArrayType>> {
        Self::default().arrays_included()
    }

    fn try_from_arrays(arrays: &BinaryArrayMap) -> Result<Vec<Self>, ArrayRetrievalError> {
        let mz_array = arrays.mzs()?;
        let intensity_array = arrays.intensities()?;
        let im_array = arrays
            .get(&ArrayType::RawIonMobilityArray)
            .or_else(|| arrays.get(&ArrayType::RawDriftTimeArray))
            .or_else(|| arrays.get(&ArrayType::RawInverseReducedIonMobilityArray))
            .ok_or(ArrayRetrievalError::NotFound(
                ArrayType::RawIonMobilityArray,
            ))?
            .to_f64()?;

        let array_key = ArrayType::nonstandard("feature identifier array");
        let marker_array = arrays
            .get(&array_key)
            .ok_or(ArrayRetrievalError::NotFound(array_key))?
            .to_i32()?;

        let n = marker_array.iter().map(|i| *i as usize).max();

        let mut features = if let Some(n) = n {
            let mut features = Vec::with_capacity(n);
            features.resize(n, Feature::default());
            features
        } else {
            return Ok(Vec::new());
        };

        mz_array
            .iter()
            .zip(intensity_array.iter())
            .zip(im_array.iter())
            .zip(marker_array.iter())
            .for_each(|(((mz, inten), im), key)| {
                features[(*key) as usize].push_raw(*mz, *im, *inten);
            });

        Ok(features)
    }
}

/// [`ChargedFeature`] with an [`IonMobility`] dimension implement
/// [`BuildArrayMapFrom`].
///
/// It uses a non-standard "feature identifier array" which assigns
/// a unique integer identifier to all points belonging to the same
/// [`Feature`].
///
/// It also produces [`ArrayType::RawIonMobilityArray`] which does not have
/// an associated [`Unit`]. The caller *should* update this if context is available.
impl BuildArrayMapFrom for ChargedFeature<Mass, IonMobility> {
    /// [`ChargedFeature`] produces [`ArrayType::MZArray`], [`ArrayType::ChargeArray`],
    /// [`ArrayType::IntensityArray`], [`ArrayType::RawIonMobilityArray`], and a
    /// non-standard "feature identifier array"
    fn arrays_included(&self) -> Option<Vec<ArrayType>> {
        Some(vec![
            ArrayType::MZArray,
            ArrayType::IntensityArray,
            ArrayType::RawIonMobilityArray,
            ArrayType::ChargeArray,
            ArrayType::nonstandard("feature identifier array"),
        ])
    }

    /// [`ChargedFeature`] uses a non-standard "feature identifier array" which assigns
    /// a unique integer identifier to all points belonging to the same
    /// [`ChargedFeature`].
    ///
    /// Like [`DeconvolutedPeak`], [`ChargedFeature`] will compute its m/z value from
    /// neutral mass and charge, and assumes that the charge carrier is a proton.
    ///
    /// It also produces [`ArrayType::RawIonMobilityArray`] which does not have
    /// an associated [`Unit`]. The caller *should* update these for consistency
    fn as_arrays(source: &[Self]) -> BinaryArrayMap {
        let mut arrays = BinaryArrayMap::new();
        let n: usize = source.iter().map(|f| f.len()).sum();

        let mut mz_array = DataArray::from_name_type_size(
            &ArrayType::MZArray,
            BinaryDataArrayType::Float64,
            n * BinaryDataArrayType::Float64.size_of(),
        );
        mz_array.unit = Unit::MZ;

        let mut intensity_array = DataArray::from_name_type_size(
            &ArrayType::IntensityArray,
            BinaryDataArrayType::Float32,
            n * BinaryDataArrayType::Float32.size_of(),
        );
        intensity_array.unit = Unit::DetectorCounts;

        let mut charge_array = DataArray::from_name_type_size(
            &ArrayType::ChargeArray,
            BinaryDataArrayType::Int32,
            n * BinaryDataArrayType::Int32.size_of(),
        );

        let mut ion_mobility_array = DataArray::from_name_type_size(
            &ArrayType::RawIonMobilityArray,
            BinaryDataArrayType::Float64,
            n * BinaryDataArrayType::Float64.size_of(),
        );

        let mut marker_array = DataArray::from_name_type_size(
            &ArrayType::nonstandard("feature identifier array"),
            BinaryDataArrayType::Int32,
            n * BinaryDataArrayType::Int32.size_of(),
        );

        let mut acc = Vec::with_capacity(n);
        source.iter().enumerate().for_each(|(i, f)| {
            f.iter().for_each(|(mass, im, inten)| {
                acc.push((mass_charge_ratio(mass, f.charge), im, inten, f.charge, i))
            })
        });
        acc.sort_by(|(mz_a, im_a, _, _, key_a), (mz_b, im_b, _, _, key_b)| {
            mz_a.total_cmp(mz_b)
                .then(im_a.total_cmp(im_b))
                .then(key_a.cmp(key_b))
        });

        for (mz, im, inten, charge, key) in acc.iter() {
            mz_array.data.extend(mz.to_le_bytes());
            intensity_array.data.extend(inten.to_le_bytes());
            ion_mobility_array.data.extend(im.to_le_bytes());
            charge_array.data.extend(charge.to_le_bytes());
            marker_array.data.extend((*key as i32).to_le_bytes());
        }

        arrays.add(mz_array);
        arrays.add(intensity_array);
        arrays.add(ion_mobility_array);
        arrays.add(charge_array);
        arrays.add(marker_array);
        arrays
    }
}

/// [`ChargedFeature`] with an [`IonMobility`] dimension implement
/// [`BuildFromArrayMap`].
///
/// It requires a non-standard "feature identifier array" which assigns
/// a unique integer identifier to all points belonging to the same
/// [`Feature`].
///
/// Like [`DeconvolutedPeak`], [`ChargedFeature`] will compute its neutral mass from
/// m/z and assumes that the charge carrier is a proton.
///
/// While it explicitly requires [`ArrayType::RawIonMobilityArray`], it will accept
/// any raw ion mobility array type.
impl BuildFromArrayMap for ChargedFeature<Mass, IonMobility> {
    /// [`ChargedFeature`] requires [`ArrayType::MZArray`], [`ArrayType::ChargeArray`] and [`ArrayType::IntensityArray`],
    /// [`ArrayType::RawIonMobilityArray`], and a non-standard "feature identifier array"
    fn arrays_required() -> Option<Vec<ArrayType>> {
        Self::default().arrays_included()
    }

    fn try_from_arrays(arrays: &BinaryArrayMap) -> Result<Vec<Self>, ArrayRetrievalError> {
        let mz_array = arrays.mzs()?;
        let intensity_array = arrays.intensities()?;
        let im_array = arrays
            .get(&ArrayType::RawIonMobilityArray)
            .or_else(|| arrays.get(&ArrayType::RawDriftTimeArray))
            .or_else(|| arrays.get(&ArrayType::RawInverseReducedIonMobilityArray))
            .ok_or(ArrayRetrievalError::NotFound(
                ArrayType::RawIonMobilityArray,
            ))?
            .to_f64()?;

        let charge_array = arrays
            .get(&ArrayType::ChargeArray)
            .ok_or(ArrayRetrievalError::NotFound(ArrayType::ChargeArray))?
            .to_i32()?;

        let array_key = ArrayType::nonstandard("feature identifier array");
        let marker_array = arrays
            .get(&array_key)
            .ok_or(ArrayRetrievalError::NotFound(array_key))?
            .to_i32()?;

        let n = marker_array.iter().map(|i| *i as usize).max();

        let mut features = if let Some(n) = n {
            let mut features = Vec::with_capacity(n);
            features.resize(n, ChargedFeature::default());
            features
        } else {
            return Ok(Vec::new());
        };

        mz_array
            .iter()
            .zip(intensity_array.iter())
            .zip(im_array.iter().zip(charge_array.iter()))
            .zip(marker_array.iter())
            .for_each(|(((mz, inten), (im, charge)), key)| {
                let f = &mut features[(*key) as usize];
                if f.is_empty() {
                    f.charge = *charge;
                }
                f.push_raw(neutral_mass(*mz, *charge), *im, *inten);
            });

        Ok(features)
    }
}

/// [`Feature`] uses the default [`BuildArrayMap3DFrom`] implementation.
impl BuildArrayMap3DFrom for Feature<MZ, IonMobility> {}

impl BuildFromArrayMap3D for Feature<MZ, IonMobility> {
    fn try_from_arrays_3d(arrays: &BinaryArrayMap3D) -> Result<Vec<Self>, ArrayRetrievalError> {
        let key = ArrayType::nonstandard("feature identifier array");
        let mut n: usize = 0;
        for (_, arr) in arrays.iter() {
            if arr.is_empty() {
                continue;
            }
            if let Some(arr) = arr.get(&key) {
                if let Some(i) = arr.iter_i32()?.map(|i| i as usize).max() {
                    n = n.max(i);
                }
            }
        }

        if n == 0 {
            return Ok(Vec::new());
        }
        n += 1;
        let mut index = Vec::with_capacity(n);
        index.resize(n, Feature::default());

        for (im, arr) in arrays.iter() {
            if arr.is_empty() {
                continue;
            }

            let mz_array = arr.mzs()?;
            let intensity_array = arr.intensities()?;
            let marker_array = arr
                .get(&key)
                .ok_or_else(|| ArrayRetrievalError::NotFound(key.clone()))?
                .to_i32()?;

            for ((mz, inten), key_i) in mz_array
                .iter()
                .zip(intensity_array.iter())
                .zip(marker_array.iter())
            {
                index[(*key_i) as usize].push_raw(*mz, im, *inten);
            }
        }

        Ok(index)
    }
}

/// [`ChargedFeature`] uses the default [`BuildArrayMap3DFrom`] implementation.
impl BuildArrayMap3DFrom for ChargedFeature<Mass, IonMobility> {}

impl BuildFromArrayMap3D for ChargedFeature<Mass, IonMobility> {
    fn try_from_arrays_3d(arrays: &BinaryArrayMap3D) -> Result<Vec<Self>, ArrayRetrievalError> {
        let key = ArrayType::nonstandard("feature identifier array");
        let mut n: usize = 0;
        for (_, arr) in arrays.iter() {
            if arr.is_empty() {
                continue;
            }
            if let Some(arr) = arr.get(&key) {
                if let Some(i) = arr.iter_i32()?.map(|i| i as usize).max() {
                    n = n.max(i);
                }
            }
        }

        if n == 0 {
            return Ok(Vec::new());
        }
        n += 1;
        let mut index = Vec::with_capacity(n);
        index.resize(n, ChargedFeature::default());

        for (im, arr) in arrays.iter() {
            if arr.is_empty() {
                continue;
            }

            let mz_array = arr.mzs()?;
            let intensity_array = arr.intensities()?;
            let marker_array = arr
                .get(&key)
                .ok_or_else(|| ArrayRetrievalError::NotFound(key.clone()))?
                .to_i32()?;
            let charge_array = arr.charges()?;

            for ((mz, inten), (charge, key_i)) in mz_array
                .iter()
                .zip(intensity_array.iter())
                .zip(charge_array.iter().zip(marker_array.iter()))
            {
                let f = &mut index[(*key_i) as usize];
                f.push_raw(*mz, im, *inten);
                f.charge = *charge;
            }
        }

        Ok(index)
    }
}

// Ion mobility-aware peaks

/// TODO
impl BuildArrayMapFrom for IonMobilityAwareCentroidPeak {
    fn arrays_included(&self) -> Option<Vec<ArrayType>> {
        Some(vec![
            ArrayType::MZArray,
            ArrayType::IntensityArray,
            ArrayType::MeanIonMobilityArray,
        ])
    }

    fn as_arrays(source: &[Self]) -> BinaryArrayMap {
        let mut arrays = BinaryArrayMap::new();

        let mut mz_array = DataArray::from_name_type_size(
            &ArrayType::MZArray,
            BinaryDataArrayType::Float64,
            source.len() * BinaryDataArrayType::Float64.size_of(),
        );
        mz_array.unit = Unit::MZ;

        let mut intensity_array = DataArray::from_name_type_size(
            &ArrayType::IntensityArray,
            BinaryDataArrayType::Float32,
            source.len() * BinaryDataArrayType::Float32.size_of(),
        );
        intensity_array.unit = Unit::DetectorCounts;

        let mut im_array = DataArray::from_name_type_size(
            &ArrayType::MeanIonMobilityArray,
            BinaryDataArrayType::Float64,
            source.len() * BinaryDataArrayType::Float64.size_of(),
        );

        for p in source.iter() {
            let mz: f64 = p.mz();
            let inten: f32 = p.intensity();
            let im = p.ion_mobility();

            let raw_bytes: [u8; mem::size_of::<f64>()] = mz.to_le_bytes();
            mz_array.data.extend_from_slice(&raw_bytes);

            let raw_bytes: [u8; mem::size_of::<f32>()] = inten.to_le_bytes();
            intensity_array.data.extend_from_slice(&raw_bytes);

            let raw_bytes: [u8; mem::size_of::<f64>()] = im.to_le_bytes();
            im_array.data.extend_from_slice(&raw_bytes);
        }

        arrays.add(mz_array);
        arrays.add(intensity_array);
        arrays.add(im_array);
        arrays
    }
}

/// TODO
impl BuildFromArrayMap for IonMobilityAwareCentroidPeak {
    fn try_from_arrays(arrays: &BinaryArrayMap) -> Result<Vec<Self>, ArrayRetrievalError> {
        let mz_array = arrays.mzs()?;
        let intensity_array = arrays.intensities()?;
        let (im_array, _) = arrays.ion_mobility()?;
        let mut peaks = Vec::with_capacity(mz_array.len());

        for (i, (mz, (intensity, ion_mobility))) in mz_array
            .iter()
            .zip(intensity_array.iter().zip(im_array.iter()))
            .enumerate()
        {
            peaks.push(IonMobilityAwareCentroidPeak {
                mz: *mz,
                intensity: *intensity,
                index: i as u32,
                ion_mobility: *ion_mobility,
            })
        }
        Ok(peaks)
    }

    fn arrays_required() -> Option<Vec<ArrayType>> {
        Some(vec![
            ArrayType::MZArray,
            ArrayType::IntensityArray,
            ArrayType::IonMobilityArray,
        ])
    }
}

/// TODO
impl BuildArrayMapFrom for IonMobilityAwareDeconvolutedPeak {
    fn arrays_included(&self) -> Option<Vec<ArrayType>> {
        Some(vec![
            ArrayType::MZArray,
            ArrayType::IntensityArray,
            ArrayType::ChargeArray,
            ArrayType::MeanIonMobilityArray,
        ])
    }

    fn as_arrays(source: &[Self]) -> BinaryArrayMap {
        let mut arrays = BinaryArrayMap::new();

        let mut mz_array = DataArray::from_name_type_size(
            &ArrayType::MZArray,
            BinaryDataArrayType::Float64,
            source.len() * BinaryDataArrayType::Float64.size_of(),
        );
        mz_array.unit = Unit::MZ;

        let mut intensity_array = DataArray::from_name_type_size(
            &ArrayType::IntensityArray,
            BinaryDataArrayType::Float32,
            source.len() * BinaryDataArrayType::Float32.size_of(),
        );
        intensity_array.unit = Unit::DetectorCounts;

        let mut charge_array = DataArray::from_name_type_size(
            &ArrayType::ChargeArray,
            BinaryDataArrayType::Int32,
            source.len() * BinaryDataArrayType::Int32.size_of(),
        );

        let mut im_array = DataArray::from_name_type_size(
            &ArrayType::MeanIonMobilityArray,
            BinaryDataArrayType::Float64,
            source.len() * BinaryDataArrayType::Float64.size_of(),
        );

        for p in source.iter() {
            let mz: f64 = p.mz();
            let inten: f32 = p.intensity();
            let im = p.ion_mobility();

            mz_array.data.extend_from_slice(&mz.to_le_bytes());
            intensity_array.data.extend_from_slice(&inten.to_le_bytes());
            im_array.data.extend_from_slice(&im.to_le_bytes());
            charge_array
                .data
                .extend_from_slice(&p.charge().to_le_bytes());
        }

        arrays.add(mz_array);
        arrays.add(intensity_array);
        arrays.add(im_array);
        arrays.add(charge_array);
        arrays
    }
}

/// TODO
impl BuildFromArrayMap for IonMobilityAwareDeconvolutedPeak {
    fn try_from_arrays(arrays: &BinaryArrayMap) -> Result<Vec<Self>, ArrayRetrievalError> {
        let mz_array = arrays.mzs()?;
        let intensity_array = arrays.intensities()?;
        let charge_array = arrays.charges()?;
        let (im_array, _) = arrays.ion_mobility()?;
        let mut peaks = Vec::with_capacity(mz_array.len());

        for (i, (mz, (intensity, (ion_mobility, charge)))) in mz_array
            .iter()
            .zip(
                intensity_array
                    .iter()
                    .zip(im_array.iter().zip(charge_array.iter())),
            )
            .enumerate()
        {
            let mass = neutral_mass(*mz, *charge);
            peaks.push(IonMobilityAwareDeconvolutedPeak {
                neutral_mass: mass,
                intensity: *intensity,
                index: i as u32,
                ion_mobility: *ion_mobility,
                charge: *charge,
            })
        }
        Ok(peaks)
    }

    fn arrays_required() -> Option<Vec<ArrayType>> {
        Some(vec![
            ArrayType::MZArray,
            ArrayType::IntensityArray,
            ArrayType::ChargeArray,
            ArrayType::IonMobilityArray,
        ])
    }
}

#[cfg(feature = "mzsignal")]
/// [`mzsignal::FittedPeak`] implements [`BuildArrayMapFrom`], which will
/// produce [`ArrayType::MZArray`], [`ArrayType::IntensityArray`], and
/// [`ArrayType::SignalToNoiseArray`]. All other properties will be lost
impl BuildArrayMapFrom for mzsignal::FittedPeak {
    fn arrays_included(&self) -> Option<Vec<ArrayType>> {
        Some(vec![
            ArrayType::MZArray,
            ArrayType::IntensityArray,
            ArrayType::SignalToNoiseArray,
        ])
    }

    fn as_arrays(source: &[Self]) -> BinaryArrayMap {
        let mut arrays = BinaryArrayMap::new();

        let mut mz_array = DataArray::from_name_type_size(
            &ArrayType::MZArray,
            BinaryDataArrayType::Float64,
            source.len() * BinaryDataArrayType::Float64.size_of(),
        );
        mz_array.unit = Unit::MZ;

        let mut intensity_array = DataArray::from_name_type_size(
            &ArrayType::IntensityArray,
            BinaryDataArrayType::Float32,
            source.len() * BinaryDataArrayType::Float32.size_of(),
        );
        intensity_array.unit = Unit::DetectorCounts;

        let mut snr_array = DataArray::from_name_type_size(
            &ArrayType::SignalToNoiseArray,
            BinaryDataArrayType::Float32,
            source.len() * BinaryDataArrayType::Float32.size_of(),
        );

        for p in source.iter() {
            let mz: f64 = p.coordinate();
            let inten: f32 = p.intensity();
            mz_array.push(mz).unwrap();
            intensity_array.push(inten).unwrap();
            snr_array.push(p.signal_to_noise).unwrap();
        }

        arrays.add(mz_array);
        arrays.add(intensity_array);
        arrays.add(snr_array);
        arrays
    }
}

#[cfg(feature = "mzsignal")]
/// [`mzsignal::FittedPeak`] implements [`BuildFromArrayMap`], which will
/// use [`ArrayType::MZArray`], [`ArrayType::IntensityArray`], and
/// look for [`ArrayType::SignalToNoiseArray`] but will default to 0.0 if it
/// is absent. All other properties are forced to 0.0.
impl BuildFromArrayMap for mzsignal::FittedPeak {
    fn arrays_required() -> Option<Vec<ArrayType>> {
        CentroidPeak::arrays_required()
    }

    fn try_from_arrays(arrays: &BinaryArrayMap) -> Result<Vec<Self>, ArrayRetrievalError> {
        let mzs = arrays.mzs()?;
        let intens = arrays.intensities()?;
        let snrs = arrays
            .get(&ArrayType::SignalToNoiseArray)
            .and_then(|a| a.to_f32().ok());

        let n = mzs.len();

        let mut out = Vec::with_capacity(n);
        for (i, (mz, inten)) in mzs.iter().copied().zip(intens.iter().copied()).enumerate() {
            out.push(mzsignal::FittedPeak::new(
                mz,
                inten,
                0,
                snrs.as_ref()
                    .and_then(|v| v.get(i).copied())
                    .unwrap_or_default(),
                0.0,
            ))
        }
        Ok(out)
    }
}
