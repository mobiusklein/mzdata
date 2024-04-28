use std::borrow::Cow;
use std::collections::hash_map::{Iter, IterMut};
use std::collections::HashMap;
use std::convert::TryFrom;

#[cfg(feature = "parallelism")]
use rayon::prelude::*;

use mzpeaks::Tolerance;

use super::array::DataArray;
use super::encodings::{ArrayRetrievalError, ArrayType, BinaryCompressionType};
use super::traits::{ByteArrayView, ByteArrayViewMut};

#[derive(Debug, Default, Clone)]
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
    pub fn iter(&self) -> Iter<ArrayType, DataArray> {
        self.byte_buffer_map.iter()
    }

    #[cfg(feature = "parallelism")]
    pub fn par_iter(&self) -> rayon::collections::hash_map::Iter<'_, ArrayType, DataArray> {
        self.byte_buffer_map.par_iter()
    }

    /// Iterate over mutable references to the key-value pairs of this map
    pub fn iter_mut(&mut self) -> IterMut<ArrayType, DataArray> {
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
    pub fn encode_array(&mut self, array_type: &ArrayType, compression: BinaryCompressionType) -> Result<(), ArrayRetrievalError> {
        if let Some(arr) = self.get_mut(array_type) {
            arr.store_compressed(compression)?;
            Ok(())
        } else {
            Err(ArrayRetrievalError::NotFound(array_type.clone()))
        }
    }

    /// Decode all [`DataArray`] in this map. If there are many arrays and the
    /// `parallelism` feature is enabled, each array will be decoded on a separate
    /// thread.
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
            .to_f64()?;
        Ok(mz_array)
    }

    /// Get a mutable reference to the m/z array if it is present
    pub fn mzs_mut(&mut self) -> Result<&mut [f64], ArrayRetrievalError> {
        if let Some(mz_array) = self.get_mut(&ArrayType::MZArray) {
            mz_array.decode_and_store()?;
            mz_array.store_as(super::BinaryDataArrayType::Float64)?;
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
            .to_f32()?;
        Ok(intensities)
    }

    /// Get a mutable reference to the intensity array if it is present
    pub fn intensities_mut(&mut self) -> Result<&mut [f32], ArrayRetrievalError> {
        if let Some(mz_array) = self.get_mut(&ArrayType::IntensityArray) {
            mz_array.decode_and_store()?;
            mz_array.store_as(super::BinaryDataArrayType::Float32)?;
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
            mz_array.store_as(super::BinaryDataArrayType::Int32)?;
            mz_array.coerce_mut()
        } else {
            Err(ArrayRetrievalError::NotFound(ArrayType::ChargeArray))
        }
    }

    /// Get a reference to the ion mobility array if it is present
    pub fn ion_mobility(&self) -> Result<(Cow<'_, [f32]>, ArrayType), ArrayRetrievalError> {
        if let Some((array_type, data_array)) = self
            .byte_buffer_map
            .iter()
            .find(|(a, _)| a.is_ion_mobility())
        {
            Ok((data_array.to_f32()?, array_type.clone()))
        } else {
            Err(ArrayRetrievalError::NotFound(ArrayType::IonMobilityArray))
        }
    }

    /// Get a mutable reference to the ion mobility array if it is present
    pub fn ion_mobility_mut(&mut self) -> Result<(&mut [f32], ArrayType), ArrayRetrievalError> {
        if let Some((array_type, data_array)) = self
            .byte_buffer_map
            .iter_mut()
            .find(|(a, _)| a.is_ion_mobility())
        {
            data_array.decode_and_store()?;
            data_array.store_as(super::BinaryDataArrayType::Float32)?;
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

#[derive(Debug, Default, Clone)]
pub struct BinaryArrayMap3D {
    pub ion_mobility_dimension: Vec<f32>,
    pub ion_mobility_type: ArrayType,
    pub arrays: Vec<BinaryArrayMap>,
    ion_mobility_index: HashMap<OrderedF32, usize>,
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
struct OrderedF32(f32);

impl std::hash::Hash for OrderedF32 {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        ((self.0 * 10000.0) as i32).hash(state);
    }
}

impl From<f32> for OrderedF32 {
    fn from(value: f32) -> Self {
        Self::wrap(value)
            .unwrap_or_else(|| panic!("Expected an order-able f32 value, but found {}", value))
    }
}

impl OrderedF32 {
    fn wrap(value: f32) -> Option<Self> {
        if value.is_nan() {
            None
        } else {
            Some(Self(value))
        }
    }
}

impl PartialOrd for OrderedF32 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for OrderedF32 {}

impl Ord for OrderedF32 {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0)
    }
}

macro_rules! _populate_stacked_array_from {
    ($im_dim:ident, $view:ident, $index_map:ident, $array_bins:ident, $array_type:ident, $array:ident) => {
        for (i_im, im) in $im_dim.iter() {
            let v = $view[*i_im];
            let i_axis = $index_map[&OrderedF32(*im)];
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

#[allow(unused)]
impl BinaryArrayMap3D {
    pub fn get_ion_mobility(&self, ion_mobility: f32) -> Option<&BinaryArrayMap> {
        if let Some(i) = OrderedF32::wrap(ion_mobility) {
            if let Some(i) = self.ion_mobility_index.get(&i) {
                self.arrays.get(*i)
            } else {
                None
            }
        } else {
            None
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = (f32, &BinaryArrayMap)> {
        self.ion_mobility_dimension
            .iter()
            .copied()
            .zip(self.arrays.iter())
    }

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
        let mut im_dim: Vec<(usize, f32)> = im_dim.iter().copied().enumerate().collect();
        im_dim.sort_by(|(_, va), (_, vb)| va.total_cmp(vb));

        let mut im_axis = Vec::with_capacity(200);
        let mut last_v = im_dim.first().unwrap().1 - 1.0;
        let mut index_map = HashMap::new();
        for (_, v) in im_dim.iter() {
            if v.total_cmp(&last_v).is_gt() {
                last_v = *v;
                index_map.insert(OrderedF32::from(*v), im_axis.len());
                im_axis.push(*v);
            }
        }

        let mut array_bins: Vec<BinaryArrayMap> = Vec::with_capacity(im_axis.len());
        array_bins.resize(im_axis.len(), BinaryArrayMap::default());
        for (array_type, array) in source.iter() {
            if array_type.is_ion_mobility() {
                continue;
            }
            match array.dtype() {
                super::BinaryDataArrayType::Unknown => {
                    panic!("Cannot re-sort opaque or unknown dimension data types")
                }
                super::BinaryDataArrayType::Float64 => {
                    let view = array.to_f64()?;
                    _populate_stacked_array_from!(
                        im_dim, view, index_map, array_bins, array_type, array
                    );
                }
                super::BinaryDataArrayType::Float32 => {
                    let view = array.to_f32()?;
                    _populate_stacked_array_from!(
                        im_dim, view, index_map, array_bins, array_type, array
                    );
                }
                super::BinaryDataArrayType::Int64 => {
                    let view = array.to_i64()?;
                    _populate_stacked_array_from!(
                        im_dim, view, index_map, array_bins, array_type, array
                    );
                }
                super::BinaryDataArrayType::Int32 => {
                    let view = array.to_i32()?;
                    _populate_stacked_array_from!(
                        im_dim, view, index_map, array_bins, array_type, array
                    );
                }
                super::BinaryDataArrayType::ASCII => {
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

#[cfg(test)]
mod test {
    use crate::spectrum::BinaryDataArrayType;

    use super::*;
    use std::fs;
    use std::io::{self, prelude::*};

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
}
