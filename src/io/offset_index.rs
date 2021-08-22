use std::io::prelude::*;

use serde::{Deserialize, Serialize};
use serde_json;

use indexmap::map::{Iter, Keys};
use indexmap::IndexMap;

#[derive(Default, Debug, Serialize, Deserialize, Clone)]
pub struct OffsetIndex {
    pub name: String,
    // If using serde_json to save this, use
    // https://docs.rs/indexmap/1.7.0/indexmap/serde_seq/index.html
    #[serde(with = "indexmap::serde_seq")]
    pub offsets: IndexMap<String, u64>,
    pub init: bool,
}

impl OffsetIndex {
    pub fn new(name: String) -> OffsetIndex {
        OffsetIndex {
            name,
            ..Default::default()
        }
    }

    #[inline]
    pub fn get(&self, key: &str) -> Option<u64> {
        self.offsets.get(key).map(|offset| *offset)
    }

    #[inline]
    pub fn get_index(&self, index: usize) -> Option<(&String, u64)> {
        if let Some((key, offset)) = self.offsets.get_index(index) {
            Some((key, *offset))
        } else {
            None
        }
    }

    #[inline]
    pub fn index_of(&self, key: &str) -> Option<usize> {
        self.offsets.get_index_of(key)
    }

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

    pub fn iter(&self) -> Iter<String, u64> {
        self.offsets.iter()
    }

    #[inline]
    pub fn contains_key(&self, key: &str) -> bool {
        self.offsets.contains_key(key)
    }

    pub fn to_writer<W: Write>(&self, writer: W) -> serde_json::Result<()> {
        serde_json::to_writer(writer, self)
    }

    pub fn from_reader<R: Read>(reader: R) -> serde_json::Result<Self> {
        serde_json::from_reader(reader)
    }
}
