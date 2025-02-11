use std::{collections::HashSet, convert::TryInto, mem};

use mzpeaks::{
    feature::{ChargedFeature, Feature},
    CentroidLike, CentroidPeak, CoordinateLike, DeconvolutedPeak, DeconvolutedPeakSet,
    IntensityMeasurement, IonMobility, KnownCharge, MZPeakSetType, Mass, PeakCollection, PeakSet,
    MZ,
};

use crate::utils::{mass_charge_ratio, neutral_mass};

use super::encodings::{
    ArrayRetrievalError, ArrayType, BinaryCompressionType, BinaryDataArrayType,
};
use super::map::BinaryArrayMap;
use super::ByteArrayView;
use super::{array::DataArray, BinaryArrayMap3D};

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

        mz_array.compression = BinaryCompressionType::Decoded;
        intensity_array.compression = BinaryCompressionType::Decoded;

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

        mz_array.compression = BinaryCompressionType::Decoded;
        intensity_array.compression = BinaryCompressionType::Decoded;
        charge_array.compression = BinaryCompressionType::Decoded;

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

#[derive(Debug, Clone)]
pub enum ArraysAvailable {
    Unknown,
    Ok,
    MissingArrays(Vec<ArrayType>),
}

pub trait BuildFromArrayMap: Sized {
    fn arrays_required() -> Option<Vec<ArrayType>> {
        None
    }

    fn try_from_arrays(arrays: &BinaryArrayMap) -> Result<Vec<Self>, ArrayRetrievalError>;

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

pub trait BuildArrayMap3DFrom: BuildArrayMapFrom {
    fn as_arrays_3d(source: &[Self]) -> BinaryArrayMap3D {
        BuildArrayMapFrom::as_arrays(source).try_into().unwrap()
    }
}

pub trait BuildFromArrayMap3D: BuildFromArrayMap {
    fn try_from_arrays_3d(arrays: &BinaryArrayMap3D) -> Result<Vec<Self>, ArrayRetrievalError> {
        BuildFromArrayMap::try_from_arrays(&arrays.unstack()?)
    }

    fn from_arrays_3d(arrays: &BinaryArrayMap3D) -> Vec<Self> {
        Self::try_from_arrays_3d(arrays).unwrap()
    }

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

pub trait BuildArrayMapFrom: Sized {
    fn arrays_included(&self) -> Option<Vec<ArrayType>> {
        None
    }

    fn as_arrays(source: &[Self]) -> BinaryArrayMap;
}

impl BuildArrayMapFrom for CentroidPeak {
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

        let mut intensity_array = DataArray::from_name_type_size(
            &ArrayType::IntensityArray,
            BinaryDataArrayType::Float32,
            source.len() * BinaryDataArrayType::Float32.size_of(),
        );

        mz_array.compression = BinaryCompressionType::Decoded;
        intensity_array.compression = BinaryCompressionType::Decoded;

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

    fn arrays_required() -> Option<Vec<ArrayType>> {
        Some(vec![ArrayType::MZArray, ArrayType::IntensityArray])
    }
}

impl BuildArrayMapFrom for DeconvolutedPeak {
    fn as_arrays(source: &[Self]) -> BinaryArrayMap {
        let mut arrays = BinaryArrayMap::new();

        let mut mz_array = DataArray::from_name_type_size(
            &ArrayType::MZArray,
            BinaryDataArrayType::Float64,
            source.len() * BinaryDataArrayType::Float64.size_of(),
        );

        let mut intensity_array = DataArray::from_name_type_size(
            &ArrayType::IntensityArray,
            BinaryDataArrayType::Float32,
            source.len() * BinaryDataArrayType::Float32.size_of(),
        );

        let mut charge_array = DataArray::from_name_type_size(
            &ArrayType::ChargeArray,
            BinaryDataArrayType::Int32,
            source.len() * BinaryDataArrayType::Int32.size_of(),
        );

        mz_array.compression = BinaryCompressionType::Decoded;
        intensity_array.compression = BinaryCompressionType::Decoded;
        charge_array.compression = BinaryCompressionType::Decoded;

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

    fn arrays_included(&self) -> Option<Vec<ArrayType>> {
        Some(vec![
            ArrayType::MZArray,
            ArrayType::IntensityArray,
            ArrayType::ChargeArray,
        ])
    }
}

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

impl BuildArrayMapFrom for Feature<MZ, IonMobility> {
    fn arrays_included(&self) -> Option<Vec<ArrayType>> {
        Some(vec![
            ArrayType::MZArray,
            ArrayType::IntensityArray,
            ArrayType::RawIonMobilityArray,
            ArrayType::nonstandard("feature identifier array"),
        ])
    }

    fn as_arrays(source: &[Self]) -> BinaryArrayMap {
        let mut arrays = BinaryArrayMap::new();
        let n: usize = source.iter().map(|f| f.len()).sum();

        let mut mz_array = DataArray::from_name_type_size(
            &ArrayType::MZArray,
            BinaryDataArrayType::Float64,
            n * BinaryDataArrayType::Float64.size_of(),
        );

        let mut intensity_array = DataArray::from_name_type_size(
            &ArrayType::IntensityArray,
            BinaryDataArrayType::Float32,
            n * BinaryDataArrayType::Float32.size_of(),
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

        mz_array.compression = BinaryCompressionType::Decoded;
        intensity_array.compression = BinaryCompressionType::Decoded;
        ion_mobility_array.compression = BinaryCompressionType::Decoded;
        marker_array.compression = BinaryCompressionType::Decoded;

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
            marker_array.data.extend_from_slice(&(*key as i32).to_le_bytes());
        }

        arrays.add(mz_array);
        arrays.add(intensity_array);
        arrays.add(ion_mobility_array);
        arrays.add(marker_array);
        arrays
    }
}

impl BuildFromArrayMap for Feature<MZ, IonMobility> {
    fn try_from_arrays(arrays: &BinaryArrayMap) -> Result<Vec<Self>, ArrayRetrievalError> {
        let mz_array = arrays.mzs()?;
        let intensity_array = arrays.intensities()?;
        let im_array = arrays
            .get(&ArrayType::RawIonMobilityArray)
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

impl BuildArrayMapFrom for ChargedFeature<Mass, IonMobility> {
    fn arrays_included(&self) -> Option<Vec<ArrayType>> {
        Some(vec![
            ArrayType::MZArray,
            ArrayType::IntensityArray,
            ArrayType::RawIonMobilityArray,
            ArrayType::ChargeArray,
            ArrayType::nonstandard("feature identifier array"),
        ])
    }

    fn as_arrays(source: &[Self]) -> BinaryArrayMap {
        let mut arrays = BinaryArrayMap::new();
        let n: usize = source.iter().map(|f| f.len()).sum();

        let mut mz_array = DataArray::from_name_type_size(
            &ArrayType::MZArray,
            BinaryDataArrayType::Float64,
            n * BinaryDataArrayType::Float64.size_of(),
        );

        let mut intensity_array = DataArray::from_name_type_size(
            &ArrayType::IntensityArray,
            BinaryDataArrayType::Float32,
            n * BinaryDataArrayType::Float32.size_of(),
        );

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

        mz_array.compression = BinaryCompressionType::Decoded;
        intensity_array.compression = BinaryCompressionType::Decoded;
        ion_mobility_array.compression = BinaryCompressionType::Decoded;
        marker_array.compression = BinaryCompressionType::Decoded;

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

impl BuildFromArrayMap for ChargedFeature<Mass, IonMobility> {
    fn try_from_arrays(arrays: &BinaryArrayMap) -> Result<Vec<Self>, ArrayRetrievalError> {
        let mz_array = arrays.mzs()?;
        let intensity_array = arrays.intensities()?;
        let im_array = arrays
            .get(&ArrayType::RawIonMobilityArray)
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
