
use std::borrow::Cow;
use std::collections::HashMap;

use log::warn;
use mzpeaks::Tolerance;

use super::bindata::DataArray;
use super::encodings::{ArrayType, ArrayRetrievalError, BinaryCompressionType};
use super::traits::ByteArrayView;


#[derive(Debug, Default, Clone)]
pub struct BinaryArrayMap {
    pub byte_buffer_map: HashMap<ArrayType, DataArray>,
}

impl<'transient, 'lifespan: 'transient> BinaryArrayMap {
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

    pub fn iter(&self) -> std::collections::hash_map::Iter<ArrayType, DataArray> {
        self.byte_buffer_map.iter()
    }

    pub fn iter_mut(&mut self) -> std::collections::hash_map::IterMut<ArrayType, DataArray> {
        self.byte_buffer_map.iter_mut()
    }

    pub fn decode(&mut self) -> Result<(), ArrayRetrievalError> {
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

    pub fn add(&mut self, array: DataArray) {
        self.byte_buffer_map.insert(array.name.clone(), array);
    }

    pub fn get(&'transient self, array_type: &ArrayType) -> Option<&'transient DataArray> {
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

    pub fn mzs(&'transient self) -> Cow<'transient, [f64]> {
        let mz_array = self
            .get(&ArrayType::MZArray)
            .expect("Did not find m/z array")
            .to_f64()
            .expect("Failed to decode m/z array");
        mz_array
    }

    pub fn intensities(&'transient self) -> Cow<'transient, [f32]> {
        let intensities = self
            .get(&ArrayType::IntensityArray)
            .expect("Did not find intensity array")
            .to_f32()
            .expect("Failed to decode intensity array");
        intensities
    }

    pub fn charges(&'transient self) -> Option<Cow<'transient, [i32]>> {
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
}
