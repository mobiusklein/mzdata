#[allow(unused)]
use std::io::prelude::*;

use indexmap::IndexMap;
use indexmap::map::{Iter, Keys};

/**
An ordered mapping from entity ID to byte offset into the source
file it resides in.

A wrapper around [`indexmap::IndexMap`].
*/
#[derive(Default, Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct OffsetIndex {
    /// The name of the index. There may potentially be more than one
    /// index per file
    pub name: String,

    /// The mapping from ID to byte offset, ordered by occurrence
    // If using serde_json to save this, use
    #[cfg_attr(feature = "serde", serde(with = "indexmap::map::serde_seq"))]
    pub offsets: IndexMap<Box<str>, u64>,

    /// Whether the index has been initalized explicitly or not, as
    /// it may be initially empty or read as empty.
    pub init: bool,
}

impl OffsetIndex {
    pub fn new(name: String) -> OffsetIndex {
        OffsetIndex {
            name,
            ..Default::default()
        }
    }

    /// Get the offset of the specified key
    #[inline]
    pub fn get(&self, key: &str) -> Option<u64> {
        self.offsets.get(key).copied()
    }

    /// Get the associated key and offset for the specified index position
    #[inline]
    pub fn get_index(&self, index: usize) -> Option<(&str, u64)> {
        if let Some((key, offset)) = self.offsets.get_index(index) {
            Some((key, *offset))
        } else {
            None
        }
    }

    /// Get the position in the index for a specific key
    #[inline]
    pub fn index_of(&self, key: &str) -> Option<usize> {
        self.offsets.get_index_of(key)
    }

    /// Insert `key` into the index with an offset value
    #[inline]
    pub fn insert<T: Into<Box<str>>>(&mut self, key: T, offset: u64) -> Option<u64> {
        self.offsets.insert(key.into(), offset)
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.offsets.len()
    }

    pub fn is_empty(&self) -> bool {
        self.offsets.is_empty()
    }

    pub fn keys(&self) -> Keys<'_, Box<str>, u64> {
        self.offsets.keys()
    }

    pub fn clear(&mut self) {
        self.offsets.clear();
    }

    /// Iterate over the keys and indices
    pub fn iter(&self) -> Iter<'_, Box<str>, u64> {
        self.offsets.iter()
    }

    /// Check if the key is in the index
    #[inline]
    pub fn contains_key(&self, key: &str) -> bool {
        self.offsets.contains_key(key)
    }

    #[cfg(feature = "serde")]
    /// Write the index out in JSON format to `writer`
    pub fn to_writer<W: Write>(&self, writer: W) -> serde_json::Result<()> {
        serde_json::to_writer(writer, self)
    }

    #[cfg(feature = "serde")]
    /// Read an index in JSON format from `reader`
    pub fn from_reader<R: Read>(reader: R) -> serde_json::Result<Self> {
        serde_json::from_reader(reader)
    }
}
