use std::io::prelude::*;

use serde::{Deserialize, Serialize};
use serde_json;

use indexmap::map::{Iter, Keys};
use indexmap::IndexMap;


/**
An ordered mapping from entity ID to byte offset into the source
file it resides in.

A wrapper around [`indexmap::IndexMap`].
*/
#[derive(Default, Debug, Serialize, Deserialize, Clone)]
pub struct OffsetIndex {
    /// The name of the index. There may potentially be more than one
    /// index per file
    pub name: String,

    /// The mapping from ID to byte offset, ordered by occurrence
    // If using serde_json to save this, use
    #[serde(with = "indexmap::map::serde_seq")]
    pub offsets: IndexMap<String, u64>,

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
    pub fn get_index(&self, index: usize) -> Option<(&String, u64)> {
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
    pub fn insert(&mut self, key: String, offset: u64) -> Option<u64> {
        self.offsets.insert(key, offset)
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.offsets.len()
    }

    pub fn is_empty(&self) -> bool {
        self.offsets.is_empty()
    }

    pub fn keys(&self) -> Keys<String, u64> {
        self.offsets.keys()
    }

    /// Iterate over the keys and indices
    pub fn iter(&self) -> Iter<String, u64> {
        self.offsets.iter()
    }

    /// Check if the key is in the index
    #[inline]
    pub fn contains_key(&self, key: &str) -> bool {
        self.offsets.contains_key(key)
    }

    /// Write the index out in JSON format to `writer`
    pub fn to_writer<W: Write>(&self, writer: W) -> serde_json::Result<()> {
        serde_json::to_writer(writer, self)
    }

    /// Read an index in JSON format from `reader`
    pub fn from_reader<R: Read>(reader: R) -> serde_json::Result<Self> {
        serde_json::from_reader(reader)
    }
}
