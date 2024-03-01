use std::borrow::Cow;
use std::collections::hash_map::{Iter, IterMut};
use std::collections::HashMap;

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

    /// Iterate over mutable references to the key-value pairs of this map
    pub fn iter_mut(&mut self) -> IterMut<ArrayType, DataArray> {
        self.byte_buffer_map.iter_mut()
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

    pub fn mzs(&'_ self) -> Result<Cow<'_, [f64]>, ArrayRetrievalError> {
        let mz_array = self
            .get(&ArrayType::MZArray)
            .ok_or(ArrayRetrievalError::NotFound(ArrayType::MZArray))?
            .to_f64()?;
        Ok(mz_array)
    }

    pub fn mzs_mut(&mut self) -> Result<&mut [f64], ArrayRetrievalError> {
        if let Some(mz_array) = self.get_mut(&ArrayType::MZArray) {
            mz_array.decode_and_store()?;
            mz_array.store_as(super::BinaryDataArrayType::Float64)?;
            mz_array.coerce_mut()
        } else {
            Err(ArrayRetrievalError::NotFound(ArrayType::MZArray))
        }
    }

    pub fn intensities(&'_ self) -> Result<Cow<'_, [f32]>, ArrayRetrievalError> {
        let intensities = self
            .get(&ArrayType::IntensityArray)
            .ok_or(ArrayRetrievalError::NotFound(ArrayType::IntensityArray))?
            .to_f32()?;
        Ok(intensities)
    }

    pub fn intensities_mut(&mut self) -> Result<&mut [f32], ArrayRetrievalError> {
        if let Some(mz_array) = self.get_mut(&ArrayType::IntensityArray) {
            mz_array.decode_and_store()?;
            mz_array.store_as(super::BinaryDataArrayType::Float32)?;
            mz_array.coerce_mut()
        } else {
            Err(ArrayRetrievalError::NotFound(ArrayType::IntensityArray))
        }
    }

    pub fn charges(&'_ self) -> Result<Cow<'_, [i32]>, ArrayRetrievalError> {
        match self.get(&ArrayType::ChargeArray) {
            Some(data_array) => data_array.to_i32(),
            None => Err(ArrayRetrievalError::NotFound(ArrayType::ChargeArray)),
        }
    }

    pub fn charge_mut(&mut self) -> Result<&mut [i32], ArrayRetrievalError> {
        if let Some(mz_array) = self.get_mut(&ArrayType::ChargeArray) {
            mz_array.decode_and_store()?;
            mz_array.store_as(super::BinaryDataArrayType::Int32)?;
            mz_array.coerce_mut()
        } else {
            Err(ArrayRetrievalError::NotFound(ArrayType::ChargeArray))
        }
    }

    pub fn ion_mobility(&self) -> Result<(Cow<'_, [f32]>, ArrayType), ArrayRetrievalError> {
        if let Some((array_type, data_array)) = self
            .byte_buffer_map
            .iter()
            .filter(|(a, _)| a.is_ion_mobility())
            .next()
        {
            Ok((data_array.to_f32()?, array_type.clone()))
        } else {
            Err(ArrayRetrievalError::NotFound(ArrayType::IonMobilityArray))
        }
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
