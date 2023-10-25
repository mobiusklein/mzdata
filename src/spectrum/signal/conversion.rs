use std::mem;

use mzpeaks::{
    CoordinateLike, DeconvolutedPeak, DeconvolutedPeakSet, IntensityMeasurement, KnownCharge,
    PeakCollection, PeakSet, CentroidLike, CentroidPeak, MZPeakSetType,
};

use crate::utils::neutral_mass;

use super::bindata::DataArray;
use super::encodings::{ArrayType, BinaryCompressionType, BinaryDataArrayType};
use super::map::BinaryArrayMap;

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

            let raw_bytes: [u8; mem::size_of::<f64>()] = unsafe { mem::transmute(mz) };
            mz_array.data.extend(raw_bytes);

            let raw_bytes: [u8; mem::size_of::<f32>()] = unsafe { mem::transmute(inten) };
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
        let mz_array = arrays.mzs();
        let intensity_array = arrays.intensities();
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
        let mz_array = arrays.mzs();
        let intensity_array = arrays.intensities();
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

            let raw_bytes: [u8; mem::size_of::<f64>()] = unsafe { mem::transmute(mz) };
            mz_array.data.extend(raw_bytes);

            let raw_bytes: [u8; mem::size_of::<f32>()] = unsafe { mem::transmute(inten) };
            intensity_array.data.extend(raw_bytes);

            let raw_bytes: [u8; mem::size_of::<i32>()] = unsafe { mem::transmute(charge) };
            charge_array.data.extend(raw_bytes);
        }

        arrays.add(mz_array);
        arrays.add(intensity_array);
        arrays.add(charge_array);
        arrays
    }
}
