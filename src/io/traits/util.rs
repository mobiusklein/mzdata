use std::io;

use crate::io::DetailLevel;

pub trait SeekRead: io::Read + io::Seek {}
impl<T: io::Read + io::Seek> SeekRead for T {}

/// Restore a temporary detail setting on early returns and dropped futures.
pub(crate) struct DetailLevelGuard<'a, T: ?Sized> {
    pub source: &'a mut T,
    saved: DetailLevel,
    restore: fn(&mut T, DetailLevel),
}

impl<'a, T: ?Sized> DetailLevelGuard<'a, T> {
    pub fn new(source: &'a mut T, saved: DetailLevel, restore: fn(&mut T, DetailLevel)) -> Self {
        Self {
            source,
            saved,
            restore,
        }
    }
}

impl<T: ?Sized> Drop for DetailLevelGuard<'_, T> {
    fn drop(&mut self) {
        (self.restore)(self.source, self.saved);
    }
}
