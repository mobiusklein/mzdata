use indexmap::map::Keys;
use indexmap::IndexMap;

#[derive(Default, Debug)]
pub struct OffsetIndex {
    pub name: String,
    // If using serde_json to save this, use
    // https://docs.rs/indexmap/1.7.0/indexmap/serde_seq/index.html
    // #[serde(with="indexmap::serde_seq")]
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

    pub fn get(&self, key: &str) -> Option<u64> {
        if let Some(offset) = self.offsets.get(key) {
            Some(*offset)
        } else {
            None
        }
    }

    pub fn get_index(&self, index: usize) -> Option<(&String, u64)> {
        if let Some((key, offset)) = self.offsets.get_index(index) {
            Some((key, *offset))
        } else {
            None
        }
    }

    pub fn insert(&mut self, key: String, offset: u64) -> Option<u64> {
        self.offsets.insert(key, offset)
    }

    pub fn len(&self) -> usize {
        self.offsets.len()
    }

    pub fn keys(&self) -> Keys<String, u64> {
        self.offsets.keys()
    }

    pub fn contains_key(&self, key: &str) -> bool {
        self.offsets.contains_key(key)
    }
}
