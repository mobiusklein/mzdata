use std::collections::{hash_map::Entry, HashMap, HashSet, VecDeque};

#[derive(Debug, Default)]
pub(crate) enum GroupIterState {
    #[default]
    Precursor,
    Product(usize),
    Done,
}

#[derive(Default, Debug)]
pub(crate) struct GenerationTracker {
    generation_to_id: HashMap<usize, HashSet<String>>,
    id_to_generation: HashMap<String, usize>,
    generations: VecDeque<usize>,
}

impl GenerationTracker {
    fn add_generation(&mut self, generation: usize) {
        match self.generations.binary_search(&generation) {
            Ok(i) => {
                self.generations.insert(i, generation);
            }
            Err(i) => {
                self.generations.insert(i, generation);
            }
        }
    }

    pub fn clear(&mut self) {
        self.generation_to_id.clear();
        self.id_to_generation.clear();
        self.generations.clear();
    }

    pub fn add(&mut self, identifier: String, generation: usize) {
        if !self.generation_to_id.contains_key(&generation) {
            self.add_generation(generation);
        }
        self.generation_to_id
            .entry(generation)
            .or_default()
            .insert(identifier.clone());
        match self.id_to_generation.entry(identifier) {
            Entry::Occupied(mut e) => {
                e.insert(generation);
            }
            Entry::Vacant(e) => {
                e.insert(generation);
            }
        }
    }

    #[allow(unused)]
    pub fn remove(&mut self, identifier: String) -> bool {
        match self.id_to_generation.entry(identifier.clone()) {
            Entry::Occupied(ent) => {
                let generation = ent.get();
                let id_set = self
                    .generation_to_id
                    .get_mut(generation)
                    .unwrap_or_else(|| {
                        panic!("Generation {} did not contain {}", generation, identifier)
                    });
                id_set.remove(&identifier);
                if id_set.is_empty() {
                    self.generations.remove(*generation);
                }
                true
            }
            Entry::Vacant(_ent) => false,
        }
    }

    pub fn older_than(&mut self, generation: usize) -> Vec<String> {
        let mut result = Vec::new();
        for gen in self.generations.iter() {
            if *gen < generation {
                if let Some(members) = self.generation_to_id.remove(gen) {
                    result.extend(members);
                }
            } else {
                break;
            }
        }
        for r in result.iter() {
            self.id_to_generation.remove(r);
        }
        result
    }
}
