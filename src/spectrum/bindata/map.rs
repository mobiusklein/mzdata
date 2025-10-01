use std::borrow::Cow;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::collections::hash_map::{Iter, IterMut};
use std::convert::TryFrom;

#[cfg(feature = "parallelism")]
use rayon::prelude::*;

use mzpeaks::Tolerance;

use crate::params::Unit;

use super::BinaryDataArrayType;
use super::array::DataArray;
use super::encodings::{ArrayRetrievalError, ArrayType, BinaryCompressionType};
use super::traits::{ByteArrayView, ByteArrayViewMut};

/// A collection of [`DataArray`]s that are identified by name.
#[derive(Debug, Default, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BinaryArrayMap {
    pub byte_buffer_map: HashMap<ArrayType, DataArray>,
}

impl BinaryArrayMap {
    pub fn new() -> BinaryArrayMap {
        BinaryArrayMap {
            ..Default::default()
        }
    }

    /// Get the number of arrays in the map
    pub fn len(&self) -> usize {
        self.byte_buffer_map.len()
    }

    pub fn is_empty(&self) -> bool {
        self.byte_buffer_map.is_empty()
    }

    /// Check if there is an ion mobility array present
    pub fn has_ion_mobility(&self) -> bool {
        self.byte_buffer_map.keys().any(|a| a.is_ion_mobility())
    }

    /// Iterate over references to the key-value pairs of this map
    pub fn iter(&self) -> Iter<'_, ArrayType, DataArray> {
        self.byte_buffer_map.iter()
    }

    #[cfg(feature = "parallelism")]
    pub fn par_iter(&self) -> rayon::collections::hash_map::Iter<'_, ArrayType, DataArray> {
        self.byte_buffer_map.par_iter()
    }

    /// Iterate over mutable references to the key-value pairs of this map
    pub fn iter_mut(&mut self) -> IterMut<'_, ArrayType, DataArray> {
        self.byte_buffer_map.iter_mut()
    }

    #[cfg(feature = "parallelism")]
    pub fn par_iter_mut(
        &mut self,
    ) -> rayon::collections::hash_map::IterMut<'_, ArrayType, DataArray> {
        self.byte_buffer_map.par_iter_mut()
    }

    /// Compress a specific [`DataArray`] with the `compression` scheme provided, if it is present.
    ///
    /// This method may fail if the compression fails or if the array type is missing.
    ///
    /// ## See Also
    /// [`DataArray::store_compressed`]
    pub fn encode_array(
        &mut self,
        array_type: &ArrayType,
        compression: BinaryCompressionType,
    ) -> Result<(), ArrayRetrievalError> {
        if let Some(arr) = self.get_mut(array_type) {
            arr.store_compressed(compression)?;
            Ok(())
        } else {
            Err(ArrayRetrievalError::NotFound(array_type.clone()))
        }
    }

    /// Decode all [`DataArray`] in this map if they have not been decoded already so
    /// that they are ready for use. If there are many arrays and the `parallelism` feature
    /// is enabled, arrays may be decoded in parallel.
    pub fn decode_all_arrays(&mut self) -> Result<(), ArrayRetrievalError> {
        #[cfg(not(feature = "parallelism"))]
        {
            self._decode_all_arrays()
        }
        #[cfg(feature = "parallelism")]
        {
            if self.len() > 2 {
                self._decode_all_arrays_parallel()
            } else {
                self._decode_all_arrays()
            }
        }
    }

    fn _decode_all_arrays(&mut self) -> Result<(), ArrayRetrievalError> {
        for (_key, value) in self.iter_mut() {
            match value.compression {
                BinaryCompressionType::Decoded => {}
                _ => {
                    value.decode_and_store()?;
                }
            }
        }
        Ok(())
    }

    #[cfg(feature = "parallelism")]
    fn _decode_all_arrays_parallel(&mut self) -> Result<(), ArrayRetrievalError> {
        let res: Result<(), ArrayRetrievalError> = self
            .iter_mut()
            .par_bridge()
            .map(|(_key, value)| {
                match value.compression {
                    BinaryCompressionType::Decoded => {}
                    _ => {
                        value.decode_and_store()?;
                    }
                }
                Ok(())
            })
            .collect::<Result<(), ArrayRetrievalError>>();
        res
    }

    /// Decode a specific [`DataArray`] if it is present.
    ///
    /// This method may fail if decoding fails or if the array type is missing.
    pub fn decode_array(&mut self, array_type: &ArrayType) -> Result<(), ArrayRetrievalError> {
        if let Some(array) = self.get_mut(array_type) {
            array.decode_and_store()?;
            Ok(())
        } else {
            Err(ArrayRetrievalError::NotFound(array_type.clone()))
        }
    }

    /// Add a [`DataArray`] to the map by its [`ArrayType`] name
    pub fn add(&mut self, array: DataArray) {
        self.byte_buffer_map.insert(array.name.clone(), array);
    }

    /// Get a reference to a specific [`DataArray`] if present
    pub fn get(&self, array_type: &ArrayType) -> Option<&DataArray> {
        self.byte_buffer_map.get(array_type)
    }

    /// Get a mutable reference to a specific [`DataArray`] if present
    pub fn get_mut(&mut self, array_type: &ArrayType) -> Option<&mut DataArray> {
        self.byte_buffer_map.get_mut(array_type)
    }

    /// Check whether a specific [`ArrayType`] is present
    pub fn has_array(&self, array_type: &ArrayType) -> bool {
        self.byte_buffer_map.contains_key(array_type)
    }

    /// Clear the map, discarding any array data
    pub fn clear(&mut self) {
        self.byte_buffer_map.clear();
    }

    /// Search for a specific m/z
    pub fn search(&self, query: f64, error_tolerance: Tolerance) -> Option<usize> {
        if let Ok(mzs) = self.mzs() {
            let (lower, _upper) = error_tolerance.bounds(query);
            match mzs[..].binary_search_by(|m| m.partial_cmp(&lower).unwrap()) {
                Ok(i) => {
                    let mut best_error = error_tolerance.call(query, mzs[i]).abs();
                    let mut best_index = i;
                    let mut index = i + 1;
                    while index < mzs.len() {
                        let error = error_tolerance.call(query, mzs[index]).abs();
                        if error < best_error {
                            best_index = index;
                            best_error = error;
                        }
                        index += 1;
                    }
                    if best_error < error_tolerance.tol() {
                        return Some(best_index);
                    }
                    None
                }
                Err(_err) => None,
            }
        } else {
            None
        }
    }

    /// Get a reference to the m/z array if it is present
    pub fn mzs(&'_ self) -> Result<Cow<'_, [f64]>, ArrayRetrievalError> {
        let mz_array = self
            .get(&ArrayType::MZArray)
            .ok_or(ArrayRetrievalError::NotFound(ArrayType::MZArray))?
            .to_f64()
            .inspect_err(|e| log::error!("Failed to decode m/z array: {e}"))?;
        Ok(mz_array)
    }

    /// Sort all the arrays in the map by `name` if they are the same length
    pub fn sort_by_array(&mut self, name: &ArrayType) -> Result<(), ArrayRetrievalError> {
        let query_axis = self
            .get(name)
            .ok_or_else(|| ArrayRetrievalError::NotFound(name.clone()))?;
        macro_rules! sort_mask {
            ($conv:expr, $cmp:expr) => {{
                let vals = $conv?;
                if vals.is_sorted() {
                    return Ok(());
                }
                let n = vals.len();
                let mut mask: Vec<usize> = (0..n).into_iter().collect();
                mask.sort_by(|i, j| {
                    let a = vals[*i];
                    let b = vals[*j];
                    $cmp(&a, &b)
                });
                (mask, n)
            }};
        }
        let (mut mask, n) = match query_axis.dtype() {
            BinaryDataArrayType::Float64 => {
                sort_mask!(query_axis.to_f64(), f64::total_cmp)
            }
            BinaryDataArrayType::Float32 => {
                sort_mask!(query_axis.to_f32(), f32::total_cmp)
            }
            BinaryDataArrayType::Int64 => {
                sort_mask!(query_axis.to_i64(), i64::cmp)
            }
            BinaryDataArrayType::Int32 => {
                sort_mask!(query_axis.to_i32(), i32::cmp)
            }
            BinaryDataArrayType::ASCII => todo!(),
            BinaryDataArrayType::Unknown => todo!(),
        };

        const TOMBSTONE: usize = usize::MAX;
        for idx in 0..n {
            if mask[idx] != TOMBSTONE {
                let mut current_idx = idx;
                loop {
                    let next_idx = mask[current_idx];
                    mask[current_idx] = TOMBSTONE;
                    if mask[next_idx] == TOMBSTONE {
                        break;
                    }
                    for (_, v) in self.iter_mut() {
                        if v.data_len()? != n {
                            continue;
                        }
                        match v.dtype {
                            BinaryDataArrayType::Float64 => {
                                let view = v.coerce_mut::<f64>()?;
                                view.swap(current_idx, next_idx);
                            }
                            BinaryDataArrayType::Float32 => {
                                let view = v.coerce_mut::<f32>()?;
                                view.swap(current_idx, next_idx);
                            }
                            BinaryDataArrayType::Int64 => {
                                let view = v.coerce_mut::<i64>()?;
                                view.swap(current_idx, next_idx);
                            }
                            BinaryDataArrayType::Int32 => {
                                let view = v.coerce_mut::<i32>()?;
                                view.swap(current_idx, next_idx);
                            }
                            BinaryDataArrayType::ASCII => todo!(),
                            BinaryDataArrayType::Unknown => todo!(),
                        }
                    }
                    current_idx = next_idx;
                }
            }
        }
        Ok(())
    }

    /// Get a mutable reference to the m/z array if it is present
    pub fn mzs_mut(&mut self) -> Result<&mut [f64], ArrayRetrievalError> {
        if let Some(mz_array) = self.get_mut(&ArrayType::MZArray) {
            mz_array
                .decode_and_store()
                .inspect_err(|e| log::error!("Failed to decode m/z array: {e}"))?;
            mz_array.store_as(BinaryDataArrayType::Float64)?;
            mz_array.coerce_mut()
        } else {
            Err(ArrayRetrievalError::NotFound(ArrayType::MZArray))
        }
    }

    /// Get a reference to the intensity array if it is present
    pub fn intensities(&'_ self) -> Result<Cow<'_, [f32]>, ArrayRetrievalError> {
        let intensities = self
            .get(&ArrayType::IntensityArray)
            .ok_or(ArrayRetrievalError::NotFound(ArrayType::IntensityArray))?
            .to_f32()
            .inspect_err(|e| log::error!("Failed to decode intensity array: {e:?}"))?;
        Ok(intensities)
    }

    /// Get a mutable reference to the intensity array if it is present
    pub fn intensities_mut(&mut self) -> Result<&mut [f32], ArrayRetrievalError> {
        if let Some(mz_array) = self.get_mut(&ArrayType::IntensityArray) {
            mz_array
                .decode_and_store()
                .inspect_err(|e| log::error!("Failed to decode intensity array: {e}"))?;
            mz_array
                .store_as(BinaryDataArrayType::Float32)
                .inspect_err(|e| log::error!("Failed to decode intensity array: {e}"))?;
            mz_array.coerce_mut()
        } else {
            Err(ArrayRetrievalError::NotFound(ArrayType::IntensityArray))
        }
    }

    /// Get a reference to the charge array if it is present
    pub fn charges(&'_ self) -> Result<Cow<'_, [i32]>, ArrayRetrievalError> {
        match self.get(&ArrayType::ChargeArray) {
            Some(data_array) => data_array.to_i32(),
            None => Err(ArrayRetrievalError::NotFound(ArrayType::ChargeArray)),
        }
    }

    /// Get a mutable reference to the charge array if it is present
    pub fn charge_mut(&mut self) -> Result<&mut [i32], ArrayRetrievalError> {
        if let Some(mz_array) = self.get_mut(&ArrayType::ChargeArray) {
            mz_array.decode_and_store()?;
            mz_array
                .store_as(BinaryDataArrayType::Int32)
                .inspect_err(|e| log::error!("Failed to decode charge array: {e}"))?;
            mz_array.coerce_mut()
        } else {
            Err(ArrayRetrievalError::NotFound(ArrayType::ChargeArray))
        }
    }

    /// Get a reference to the ion mobility array if it is present
    pub fn ion_mobility(&self) -> Result<(Cow<'_, [f64]>, ArrayType), ArrayRetrievalError> {
        if let Some((array_type, data_array)) = self
            .byte_buffer_map
            .iter()
            .find(|(a, _)| a.is_ion_mobility())
        {
            Ok((
                data_array
                    .to_f64()
                    .inspect_err(|e| log::error!("Failed to decode ion mobility array: {e}"))?,
                array_type.clone(),
            ))
        } else {
            Err(ArrayRetrievalError::NotFound(ArrayType::IonMobilityArray))
        }
    }

    /// Get a mutable reference to the ion mobility array if it is present
    pub fn ion_mobility_mut(&mut self) -> Result<(&mut [f64], ArrayType), ArrayRetrievalError> {
        if let Some((array_type, data_array)) = self
            .byte_buffer_map
            .iter_mut()
            .find(|(a, _)| a.is_ion_mobility())
        {
            data_array
                .decode_and_store()
                .inspect_err(|e| log::error!("Failed to decode ion mobility array: {e}"))?;
            data_array.store_as(BinaryDataArrayType::Float32)?;
            Ok((data_array.coerce_mut()?, array_type.clone()))
        } else {
            Err(ArrayRetrievalError::NotFound(ArrayType::IonMobilityArray))
        }
    }

    /// Split an array map along the ion mobility dimension if it is present.
    ///
    /// This function will fail if there is no ion mobility dimension.
    pub fn stack_ion_mobility(self) -> Result<BinaryArrayMap3D, ArrayRetrievalError> {
        BinaryArrayMap3D::try_from(self)
    }
}

impl IntoIterator for BinaryArrayMap {
    type Item = (ArrayType, DataArray);

    type IntoIter = <HashMap<ArrayType, DataArray> as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        self.byte_buffer_map.into_iter()
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
struct NonNaNF64(f64);

impl std::hash::Hash for NonNaNF64 {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        ((self.0 * 10000.0) as i64).hash(state);
    }
}

impl From<f64> for NonNaNF64 {
    fn from(value: f64) -> Self {
        Self::wrap(value)
            .unwrap_or_else(|| panic!("Expected an order-able f64 value, but found {}", value))
    }
}

impl NonNaNF64 {
    fn wrap(value: f64) -> Option<Self> {
        if value.is_nan() {
            None
        } else {
            Some(Self(value))
        }
    }
}

impl PartialOrd for NonNaNF64 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for NonNaNF64 {}

impl Ord for NonNaNF64 {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

macro_rules! _populate_stacked_array_from {
    ($im_dim:ident, $view:ident, $index_map:ident, $array_bins:ident, $array_type:ident, $array:ident) => {
        for (i_im, im) in $im_dim.iter() {
            let v = $view[*i_im];
            let i_axis = $index_map[&NonNaNF64(*im)];
            let bin = &mut $array_bins[i_axis];
            if let Some(bin_array) = bin.get_mut($array_type) {
                bin_array.push(v)?;
            } else {
                let mut bin_array = DataArray::from_name_and_type($array_type, $array.dtype());
                bin_array.push(v)?;
                bin.add(bin_array);
            }
        }
    };
}

/// Represent a set of [`BinaryArrayMap`] that has been split across the
/// ion mobility dimension.
#[derive(Debug, Default, Clone)]
#[cfg_attr(feature = "serde", serde_with::serde_as)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BinaryArrayMap3D {
    pub ion_mobility_dimension: Vec<f64>,
    pub ion_mobility_type: ArrayType,
    pub ion_mobility_unit: Unit,
    pub arrays: Vec<BinaryArrayMap>,
    pub additional_arrays: BinaryArrayMap,
    #[cfg_attr(feature = "serde", serde_as(as = "Vec<(_, _)>"))]
    ion_mobility_index: HashMap<NonNaNF64, usize>,
}

impl BinaryArrayMap3D {
    pub fn from_ion_mobility_dimension(
        ion_mobility_dimension: Vec<f64>,
        ion_mobility_type: ArrayType,
        ion_mobility_unit: Unit,
    ) -> BinaryArrayMap3D {
        let ion_mobility_index: HashMap<_, _> = ion_mobility_dimension
            .iter()
            .copied()
            .enumerate()
            .map(|(i, f)| {
                (
                    NonNaNF64::wrap(f)
                        .unwrap_or_else(|| panic!("Expected non-NaN value for ion mobility")),
                    i,
                )
            })
            .collect();
        let mut arrays = Vec::new();
        arrays.resize_with(ion_mobility_dimension.len(), BinaryArrayMap::default);
        Self {
            ion_mobility_dimension,
            ion_mobility_type,
            ion_mobility_unit,
            arrays,
            ion_mobility_index,
            additional_arrays: Default::default(),
        }
    }

    pub fn from_ion_mobility_dimension_and_arrays(
        ion_mobility_dimension: Vec<f64>,
        ion_mobility_type: ArrayType,
        ion_mobility_unit: Unit,
        arrays: Vec<BinaryArrayMap>,
    ) -> BinaryArrayMap3D {
        let ion_mobility_index: HashMap<_, _> = ion_mobility_dimension
            .iter()
            .copied()
            .enumerate()
            .map(|(i, f)| {
                (
                    NonNaNF64::wrap(f)
                        .unwrap_or_else(|| panic!("Expected non-NaN value for ion mobility")),
                    i,
                )
            })
            .collect();

        Self {
            ion_mobility_dimension,
            ion_mobility_type,
            ion_mobility_unit,
            arrays,
            ion_mobility_index,
            additional_arrays: Default::default(),
        }
    }

    /// Get the associated arrays at the requested ion mobility, if they exist
    pub fn get_ion_mobility(&self, ion_mobility: f64) -> Option<&BinaryArrayMap> {
        if let Some(i) = NonNaNF64::wrap(ion_mobility) {
            if let Some(i) = self.ion_mobility_index.get(&i) {
                self.arrays.get(*i)
            } else {
                None
            }
        } else {
            None
        }
    }

    pub fn search_ion_mobility(
        &self,
        ion_mobility: f64,
        error_tolerance: f64,
    ) -> Option<(&BinaryArrayMap, f64)> {
        match self
            .ion_mobility_dimension
            .binary_search_by(|x: &f64| x.total_cmp(&ion_mobility))
        {
            Ok(i) => {
                let delta = ion_mobility - self.ion_mobility_dimension[i];
                if delta.abs() <= error_tolerance {
                    self.arrays.get(i).map(|a| (a, delta))
                } else {
                    None
                }
            }
            Err(i) => {
                if self.arrays.is_empty() {
                    return None;
                }
                let delta = ion_mobility - self.ion_mobility_dimension[i];
                if delta.abs() <= error_tolerance {
                    self.arrays.get(i).map(|a| (a, delta))
                } else {
                    None
                }
            }
        }
    }

    /// Get the a mutable reference to the associated arrays at the requested ion mobility, if they exist
    pub fn get_ion_mobility_mut(&mut self, ion_mobility: f64) -> Option<&mut BinaryArrayMap> {
        if let Some(i) = NonNaNF64::wrap(ion_mobility) {
            if let Some(i) = self.ion_mobility_index.get(&i) {
                self.arrays.get_mut(*i)
            } else {
                None
            }
        } else {
            None
        }
    }

    /// Iterate over the ion mobility dimension and associated arrays
    /// at each point.
    pub fn iter(&self) -> impl Iterator<Item = (f64, &BinaryArrayMap)> {
        self.ion_mobility_dimension
            .iter()
            .copied()
            .zip(self.arrays.iter())
    }

    /// Iterate over the ion mobility dimension and a mutable reference to the associated arrays
    /// at each point.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = (f64, &mut BinaryArrayMap)> {
        self.ion_mobility_dimension
            .iter()
            .copied()
            .zip(self.arrays.iter_mut())
    }

    /// Flatten this array into a single [`BinaryArrayMap`].
    ///
    /// # Errors
    /// [`ArrayRetrievalError`] errors related to array decoding occur if
    /// any [`DataArray`] cannot be decoded, or if an expected array is absent.
    pub fn unstack(&self) -> Result<BinaryArrayMap, ArrayRetrievalError> {
        let mut destination = self.additional_arrays.clone();

        let mut im_dim =
            DataArray::from_name_and_type(&self.ion_mobility_type, BinaryDataArrayType::Float64);
        im_dim.unit = self.ion_mobility_unit;

        destination.add(im_dim);

        let mut current_mz = f64::INFINITY;
        let mut max_mz = f64::NEG_INFINITY;
        let mut mz_axes = Vec::new();
        let mut indices = Vec::new();

        // Prepare the destination to have a waiting array of all types in the layers of this
        // array stack.
        let mut total_points = 0;
        for layer in self.arrays.iter() {
            for (key, array) in layer.iter() {
                if !destination.has_array(key) {
                    let mut new_array = DataArray::from_name_and_type(key, array.dtype);
                    new_array.unit = array.unit;
                    new_array.params = array.params.clone();
                    new_array.compression = BinaryCompressionType::Decoded;
                    destination.add(new_array);
                }
            }
            let mzs = layer.mzs()?;
            total_points += mzs.len();
            if let Some(mz) = mzs.first() {
                current_mz = mz.min(current_mz);
            }
            if let Some(mz) = mzs.last() {
                max_mz = mz.max(max_mz);
            }
            mz_axes.push(mzs);
            indices.push(0usize);
        }

        let mut n_points_added = 0;
        while current_mz <= max_mz {
            let mut next_mz = f64::INFINITY;
            for (bin_i, (im, layer)) in self.iter().enumerate() {
                if let Some(i) = indices.get(bin_i).copied()
                    && let Some(mz) = mz_axes[bin_i].get(i).copied()
                {
                    if (mz - current_mz).abs() < 1e-12 {
                        n_points_added += 1;
                        destination
                            .get_mut(&ArrayType::MZArray)
                            .as_mut()
                            .unwrap()
                            .push(mz)?;
                        destination
                            .get_mut(&self.ion_mobility_type)
                            .as_mut()
                            .unwrap()
                            .push(im)?;
                        for (key, array) in layer.iter() {
                            if *key == ArrayType::MZArray {
                                continue;
                            }
                            match array.dtype() {
                                BinaryDataArrayType::Unknown => {
                                    panic!("Cannot re-sort opaque or unknown dimension data types")
                                }
                                BinaryDataArrayType::ASCII => {
                                    let val = array.decode()?[i];
                                    destination.get_mut(key).as_mut().unwrap().push(val)?;
                                }
                                BinaryDataArrayType::Float64 => {
                                    let val = array.to_f64()?[i];
                                    destination.get_mut(key).as_mut().unwrap().push(val)?;
                                }
                                BinaryDataArrayType::Float32 => {
                                    let val = array.to_f32()?[i];
                                    destination.get_mut(key).as_mut().unwrap().push(val)?;
                                }
                                BinaryDataArrayType::Int64 => {
                                    let val = array.to_i64()?[i];
                                    destination.get_mut(key).as_mut().unwrap().push(val)?;
                                }
                                BinaryDataArrayType::Int32 => {
                                    let val = array.to_i32()?[i];
                                    destination.get_mut(key).as_mut().unwrap().push(val)?;
                                }
                            }
                        }
                        indices[bin_i] += 1;
                        if let Some(mz) = mz_axes[bin_i].get(i + 1).copied() {
                            next_mz = mz.min(next_mz);
                        }
                    } else if mz > current_mz {
                        next_mz = mz.min(next_mz);
                    }
                }
            }
            current_mz = next_mz;
        }

        debug_assert_eq!(
            n_points_added,
            total_points,
            "Expected to have unstacked {total_points} from {} arrays, got {n_points_added} instead",
            mz_axes.len()
        );

        Ok(destination)
    }

    /// Convert a [`BinaryArrayMap`] into a [`BinaryArrayMap3D`] if it has an ion mobility dimension.
    ///
    /// Any arrays that aren't the same length as the ion mobility dimension will be in
    /// [`BinaryArrayMap3D::additional_arrays`].
    ///
    /// # Errors
    /// [`ArrayRetrievalError`] errors related to array decoding occur if any [`DataArray`]
    /// cannot be decoded, or if an expected array is absent.
    pub fn stack(source: &BinaryArrayMap) -> Result<Self, ArrayRetrievalError> {
        let mut this = Self::default();
        if !source.has_ion_mobility() {
            return Err(ArrayRetrievalError::NotFound(ArrayType::IonMobilityArray));
        }
        let (im_dim, im_type) = source.ion_mobility()?;
        this.ion_mobility_type = im_type;
        if im_dim.is_empty() {
            return Ok(this);
        }
        let mut im_dim: Vec<(usize, f64)> = im_dim.iter().copied().enumerate().collect();
        im_dim.sort_by(|(_, va), (_, vb)| va.total_cmp(vb));

        let mut im_axis = Vec::with_capacity(200);
        let mut last_v = im_dim.first().unwrap().1 - 1.0;
        let mut index_map = HashMap::new();
        for (_, v) in im_dim.iter() {
            if v.total_cmp(&last_v).is_gt() {
                last_v = *v;
                index_map.insert(NonNaNF64::from(*v), im_axis.len());
                im_axis.push(*v);
            }
        }

        let mut array_bins: Vec<BinaryArrayMap> = Vec::with_capacity(im_axis.len());
        array_bins.resize(im_axis.len(), BinaryArrayMap::default());
        for (array_type, array) in source.iter() {
            if array_type.is_ion_mobility() {
                continue;
            }
            if array.data_len()? != im_dim.len() {
                this.additional_arrays.add(array.clone());
                continue;
            }
            match array.dtype() {
                BinaryDataArrayType::Unknown => {
                    panic!("Cannot re-sort opaque or unknown dimension data types")
                }
                BinaryDataArrayType::Float64 => {
                    let view = array.to_f64()?;
                    _populate_stacked_array_from!(
                        im_dim, view, index_map, array_bins, array_type, array
                    );
                }
                BinaryDataArrayType::Float32 => {
                    let view = array.to_f32()?;
                    _populate_stacked_array_from!(
                        im_dim, view, index_map, array_bins, array_type, array
                    );
                }
                BinaryDataArrayType::Int64 => {
                    let view = array.to_i64()?;
                    _populate_stacked_array_from!(
                        im_dim, view, index_map, array_bins, array_type, array
                    );
                }
                BinaryDataArrayType::Int32 => {
                    let view = array.to_i32()?;
                    _populate_stacked_array_from!(
                        im_dim, view, index_map, array_bins, array_type, array
                    );
                }
                BinaryDataArrayType::ASCII => {
                    let view = array.decode()?;
                    _populate_stacked_array_from!(
                        im_dim, view, index_map, array_bins, array_type, array
                    );
                }
            }
        }
        this.ion_mobility_dimension = im_axis;
        this.ion_mobility_index = index_map;
        this.arrays = array_bins;
        Ok(this)
    }
}

impl TryFrom<BinaryArrayMap> for BinaryArrayMap3D {
    type Error = ArrayRetrievalError;

    fn try_from(value: BinaryArrayMap) -> Result<Self, Self::Error> {
        Self::stack(&value)
    }
}

impl TryFrom<&BinaryArrayMap> for BinaryArrayMap3D {
    type Error = ArrayRetrievalError;

    fn try_from(value: &BinaryArrayMap) -> Result<Self, Self::Error> {
        Self::stack(value)
    }
}

#[cfg(test)]
mod test {
    use crate::prelude::*;
    use crate::spectrum::BinaryDataArrayType;

    use super::*;
    use std::fs;
    use std::io;

    fn make_array_from_file() -> io::Result<DataArray> {
        let mut fh = fs::File::open("./test/data/mz_f64_zlib_bas64.txt")?;
        let mut buf = String::new();
        fh.read_to_string(&mut buf)?;
        let bytes: Vec<u8> = buf.into();
        let mut da = DataArray::wrap(&ArrayType::MZArray, BinaryDataArrayType::Float64, bytes);
        da.compression = BinaryCompressionType::Zlib;
        Ok(da)
    }

    #[test]
    fn test_construction() -> io::Result<()> {
        let da = make_array_from_file()?;
        let mut map = BinaryArrayMap::new();
        assert!(!map.has_array(&ArrayType::MZArray));
        map.add(da);
        assert!(map.has_array(&ArrayType::MZArray));
        Ok(())
    }

    #[test]
    fn test_decode() -> io::Result<()> {
        let da = make_array_from_file()?;
        let mut map = BinaryArrayMap::new();
        map.add(da);
        assert_eq!(
            map.get(&ArrayType::MZArray).unwrap().compression,
            BinaryCompressionType::Zlib
        );
        map.decode_all_arrays()?;
        assert_eq!(
            map.get(&ArrayType::MZArray).unwrap().compression,
            BinaryCompressionType::Decoded
        );
        Ok(())
    }

    #[test]
    fn test_3d_stack_unstack() -> io::Result<()> {
        let mut reader = crate::MZReader::open_gzipped_read(io::BufReader::new(fs::File::open(
            "test/data/20200204_BU_8B8egg_1ug_uL_7charges_60_min_Slot2-11_1_244.mzML.gz",
        )?))?;

        let spec = reader
            .get_spectrum_by_id("merged=42869 frame=9717 scanStart=1 scanEnd=705")
            .unwrap();
        let mut arrays = spec.arrays.unwrap();
        let mzs = arrays.mzs()?;
        assert!(!mzs.is_sorted());
        drop(mzs);
        arrays.sort_by_array(&ArrayType::MZArray)?;
        let mzs = arrays.mzs()?;
        assert!(mzs.is_sorted());
        let n = mzs.len();

        let arrays_3d = BinaryArrayMap3D::stack(&arrays)?;
        let stacked_n: usize = arrays_3d
            .iter()
            .map(|(_, va)| va.mzs().unwrap().len())
            .sum();

        assert_eq!(n, stacked_n);

        let unstacked = arrays_3d.unstack()?;
        let unstacked_n = unstacked.mzs()?.len();

        assert_eq!(unstacked_n, n);

        Ok(())
    }
}
