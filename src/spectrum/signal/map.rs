
use std::borrow::Cow;
use std::collections::HashMap;
use std::collections::hash_map::{Iter, IterMut};

#[cfg(feature = "parallelism")]
use rayon::prelude::*;

use log::warn;
use mzpeaks::Tolerance;

use super::bindata::DataArray;
use super::encodings::{ArrayType, ArrayRetrievalError, BinaryCompressionType};
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

    pub fn len(&self) -> usize {
        self.byte_buffer_map.len()
    }

    pub fn is_empty(&self) -> bool {
        self.byte_buffer_map.is_empty()
    }

    pub fn iter(&self) -> Iter<ArrayType, DataArray> {
        self.byte_buffer_map.iter()
    }

    pub fn iter_mut(&mut self) -> IterMut<ArrayType, DataArray> {
        self.byte_buffer_map.iter_mut()
    }

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
    fn _decode_all_arrays_parallel(&mut self)  -> Result<(), ArrayRetrievalError> {
        let res: Result<(), ArrayRetrievalError> = self.iter_mut().par_bridge().map(|(_key, value)| {
            match value.compression {
                BinaryCompressionType::Decoded => {}
                _ => {
                    value.decode_and_store()?;
                }
            }
            Ok(())
        }).collect::<Result<(), ArrayRetrievalError>>();
        res
    }

    pub fn decode_array(&mut self, array_type: &ArrayType) -> Result<(), ArrayRetrievalError> {
        if let Some(array) = self.get_mut(array_type) {
            array.decode_and_store()?;
            Ok(())
        } else {
            Err(ArrayRetrievalError::NotFound)
        }
    }

    pub fn add(&mut self, array: DataArray) {
        self.byte_buffer_map.insert(array.name.clone(), array);
    }

    pub fn get(&self, array_type: &ArrayType) -> Option<&DataArray> {
        self.byte_buffer_map.get(array_type)
    }

    pub fn get_mut(&mut self, array_type: &ArrayType) -> Option<&mut DataArray> {
        self.byte_buffer_map.get_mut(array_type)
    }

    pub fn has_array(&self, array_type: &ArrayType) -> bool {
        self.byte_buffer_map.contains_key(array_type)
    }

    pub fn clear(&mut self) {
        self.byte_buffer_map.clear();
    }

    pub fn search(&self, query: f64, error_tolerance: Tolerance) -> Option<usize> {
        let mzs = self.mzs();
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
    }

    pub fn mzs(&'_ self) -> Cow<'_, [f64]> {
        let mz_array = self
            .get(&ArrayType::MZArray)
            .expect("Did not find m/z array")
            .to_f64()
            .expect("Failed to decode m/z array");
        mz_array
    }

    pub fn mzs_mut(&mut self) -> Result<&mut [f64], ArrayRetrievalError> {
        if let Some(mz_array) = self.get_mut(&ArrayType::MZArray) {
            mz_array.decode_and_store()?;
            mz_array.store_as(super::BinaryDataArrayType::Float64)?;
            mz_array.coerce_mut()
        } else {
            Err(ArrayRetrievalError::NotFound)
        }
    }

    pub fn intensities(&'_ self) -> Cow<'_, [f32]> {
        let intensities = self
            .get(&ArrayType::IntensityArray)
            .expect("Did not find intensity array")
            .to_f32()
            .expect("Failed to decode intensity array");
        intensities
    }

    pub fn intensities_mut(&mut self) -> Result<&mut [f32], ArrayRetrievalError> {
        if let Some(mz_array) = self.get_mut(&ArrayType::IntensityArray) {
            mz_array.decode_and_store()?;
            mz_array.store_as(super::BinaryDataArrayType::Float32)?;
            mz_array.coerce_mut()
        } else {
            Err(ArrayRetrievalError::NotFound)
        }
    }

    pub fn charges(&'_ self) -> Option<Cow<'_, [i32]>> {
        match self.get(&ArrayType::ChargeArray) {
            Some(data_array) => match data_array.to_i32() {
                Ok(array) => Some(array),
                Err(err) => {
                    warn!("Failed to decode charge state array: {:?}", err);
                    None
                }
            },
            None => None,
        }
    }

    pub fn charge_mut(&mut self) -> Result<&mut [i32], ArrayRetrievalError> {
        if let Some(mz_array) = self.get_mut(&ArrayType::ChargeArray) {
            mz_array.decode_and_store()?;
            mz_array.store_as(super::BinaryDataArrayType::Int32)?;
            mz_array.coerce_mut()
        } else {
            Err(ArrayRetrievalError::NotFound)
        }
    }

}


#[cfg(test)]
mod test {
    use crate::spectrum::BinaryDataArrayType;

    use super::*;
    use std::io::{self, prelude::*};
    use std::fs;

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
        assert_eq!(map.get(&ArrayType::MZArray).unwrap().compression, BinaryCompressionType::Zlib);
        map.decode_all_arrays()?;
        assert_eq!(map.get(&ArrayType::MZArray).unwrap().compression, BinaryCompressionType::Decoded);
        Ok(())
    }
}