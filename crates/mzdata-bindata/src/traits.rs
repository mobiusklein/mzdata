use std::marker::PhantomData;
use std::mem;
use std::borrow::Cow;

use bytemuck::Pod;
use num_traits::{AsPrimitive, Num};
use mzdata_param::Unit;

use super::encodings::{ArrayRetrievalError, BinaryDataArrayType, Bytes};
use super::ArrayType;

/// A byte array that may be encoded and compressed as [`BinaryCompressionType`], or unpacked in
/// memory and re-interpreted as native types appropriate to [`BinaryDataArrayType`].
///
/// The raw bytes must be made available through the [`ByteArrayView::view`] method,
/// but all other behavior is added on top of that
pub trait ByteArrayView<'transient, 'lifespan: 'transient> {
    /// The core method that exposes the decoded byte array that the type *must* provide.
    fn view(&'lifespan self) -> Result<Cow<'lifespan, [u8]>, ArrayRetrievalError>;

    /// This is a helper method for the various `to_<type>` methods like [`Self::to_f64`]. It should not
    /// need to be called directly. This method relies on [`bytemuck::try_cast_slice`] to do the heavy lifting.
    fn coerce_from<T: Pod>(
        buffer: Cow<'transient, [u8]>,
    ) -> Result<Cow<'transient, [T]>, ArrayRetrievalError> {
        let n = buffer.len();
        if n == 0 {
            return Ok(Cow::Owned(Vec::new()))
        }
        let z = mem::size_of::<T>();
        if n % z != 0 {
            return Err(ArrayRetrievalError::DataTypeSizeMismatch);
        }
        match buffer {
            Cow::Borrowed(c) => {
                Ok(Cow::Borrowed(bytemuck::try_cast_slice(c)?))
            },
            Cow::Owned(v) => {
                let size_type = n / z;
                let mut buf = Vec::with_capacity(size_type);
                v.chunks_exact(z).try_for_each(|c| {
                    buf.extend(bytemuck::try_cast_slice(c)?);
                    Ok::<(), bytemuck::PodCastError>(())
                })?;
                Ok(Cow::Owned(buf))
            },
        }
    }

    /// This is a helper method for the various `to_<type>` methods like [`Self::to_f64`]. It should not
    /// need to be called directly. This calls into [`Self::coerce_from`] internally after decoding.
    fn coerce<T: Pod>(
        &'lifespan self,
    ) -> Result<Cow<'transient, [T]>, ArrayRetrievalError> {
        match self.view() {
            Ok(data) => {
                #[cfg(target_endian = "big")]
                {
                    let mut data: Cow<'_, [u8]> = Cow::Owned(data.to_vec());
                    self.dtype().swap_bytes(data.to_mut())?;
                    Self::coerce_from(data)
                }
                #[cfg(not(target_endian = "big"))]
                {
                    Self::coerce_from(data)
                }
            },
            Err(err) => Err(err),
        }
    }

    /// Decode the array, then copy it to a new array, converting each element from type `D` to to type `S`
    ///
    /// This makes a copy of the data, forcibly coercing the values using [`AsPrimitive`]
    fn convert<S: Num + Clone + AsPrimitive<D> + Pod, D: Num + Clone + Copy + 'static>(
        &'lifespan self,
    ) -> Result<Cow<'transient, [D]>, ArrayRetrievalError> {
        match self.coerce::<S>() {
            Ok(view) => {
                match view {
                    Cow::Borrowed(view) => {
                        Ok(Cow::Owned(view.iter().map(|a| a.as_()).collect()))
                    }
                    Cow::Owned(owned) => {
                        let res = owned.iter().map(|a| a.as_()).collect();
                        Ok(Cow::Owned(res))
                    }
                }
            }
            Err(err) => Err(err),
        }
    }

    /// The kind of array this is
    fn name(&self) -> &ArrayType;

    /// The real data type encoded in bytes
    fn dtype(&self) -> BinaryDataArrayType;

    /// The unit of measurement each data point is in
    fn unit(&self) -> Unit;

    /// Get the identifier referencing a [`DataProcessing`](crate::meta::DataProcessing)
    fn data_processing_reference(&self) -> Option<&str> {
        None
    }

    /// Get a view of the data as [`f32`]. If the data are stored that way already, no
    /// copy is required.
    fn to_f32(&'lifespan self) -> Result<Cow<'transient, [f32]>, ArrayRetrievalError> {
        type D = f32;
        match self.dtype() {
            BinaryDataArrayType::Float32 | BinaryDataArrayType::ASCII => self.coerce::<D>(),
            BinaryDataArrayType::Float64 => {
                type S = f64;
                self.convert::<S, D>()
            }
            BinaryDataArrayType::Int32 => {
                type S = i32;
                self.convert::<S, D>()
            }
            BinaryDataArrayType::Int64 => {
                type S = i64;
                self.convert::<S, D>()
            }
            _ => Err(ArrayRetrievalError::DataTypeSizeMismatch),
        }
    }

    /// Get a view of the data as [`f64`]. If the data are stored that way already, no
    /// copy is required.
    fn to_f64(&'lifespan self) -> Result<Cow<'transient, [f64]>, ArrayRetrievalError> {
        type D = f64;
        match self.dtype() {
            BinaryDataArrayType::Float32 => {
                type S = f32;
                self.convert::<S, D>()
            }
            BinaryDataArrayType::Float64 | BinaryDataArrayType::ASCII => self.coerce(),
            BinaryDataArrayType::Int32 => {
                type S = i32;
                self.convert::<S, D>()
            }
            BinaryDataArrayType::Int64 => {
                type S = i64;
                self.convert::<S, D>()
            }
            _ => Err(ArrayRetrievalError::DataTypeSizeMismatch),
        }
    }

    /// Get a view of the data as [`i32`]. If the data are stored that way already, no
    /// copy is required.
    fn to_i32(&'lifespan self) -> Result<Cow<'transient, [i32]>, ArrayRetrievalError> {
        type D = i32;
        match self.dtype() {
            BinaryDataArrayType::Float32 => {
                type S = f32;
                self.convert::<S, D>()
            }
            BinaryDataArrayType::Float64 => {
                type S = f64;
                self.convert::<S, D>()
            }
            BinaryDataArrayType::Int32 | BinaryDataArrayType::ASCII => self.coerce::<D>(),
            BinaryDataArrayType::Int64 => {
                type S = i64;
                self.convert::<S, D>()
            }
            _ => Err(ArrayRetrievalError::DataTypeSizeMismatch),
        }
    }

    /// Get a view of the data as [`i64`]. If the data are stored that way already, no
    /// copy is required.
    fn to_i64(&'lifespan self) -> Result<Cow<'transient, [i64]>, ArrayRetrievalError> {
        type D = i64;
        match self.dtype() {
            BinaryDataArrayType::Float32 => {
                type S = f32;
                self.convert::<S, D>()
            }
            BinaryDataArrayType::Float64 => {
                type S = f64;
                self.convert::<S, D>()
            }
            BinaryDataArrayType::Int64 | BinaryDataArrayType::ASCII => self.coerce::<D>(),
            BinaryDataArrayType::Int32 => {
                type S = i32;
                self.convert::<S, D>()
            }
            _ => Err(ArrayRetrievalError::DataTypeSizeMismatch),
        }
    }

    /// The size of encoded array in terms of # of elements of the [`BinaryDataArrayType`] given by [`ByteArrayView::dtype`]
    fn data_len(&'lifespan self) -> Result<usize, ArrayRetrievalError> {
        let view = self.view()?;
        let n = view.len();
        Ok(n / self.dtype().size_of())
    }

    fn iter_type<T: Pod>(&'lifespan self) -> Result<DataSliceIter<'lifespan, T>, ArrayRetrievalError> {
        Ok(DataSliceIter::new(self.view()?))
    }

    fn iter_u8(&'lifespan self) -> Result<DataSliceIter<'lifespan, u8>, ArrayRetrievalError> {
        Ok(DataSliceIter::new(self.view()?))
    }

    fn iter_f32(&'lifespan self) -> Result<DataSliceIter<'lifespan, f32>, ArrayRetrievalError> {
        Ok(DataSliceIter::new(self.view()?))
    }

    fn iter_f64(&'lifespan self) -> Result<DataSliceIter<'lifespan, f64>, ArrayRetrievalError> {
        Ok(DataSliceIter::new(self.view()?))
    }

    fn iter_i32(&'lifespan self) -> Result<DataSliceIter<'lifespan, i32>, ArrayRetrievalError> {
        Ok(DataSliceIter::new(self.view()?))
    }

    fn iter_i64(&'lifespan self) -> Result<DataSliceIter<'lifespan, i64>, ArrayRetrievalError> {
        Ok(DataSliceIter::new(self.view()?))
    }
}

/// A mutable byte array
pub trait ByteArrayViewMut<'transient, 'lifespan: 'transient>:
    ByteArrayView<'transient, 'lifespan>
{

    /// Specify the unit of the data array
    fn unit_mut(&mut self) -> &mut Unit;

    /// Get a mutable view of the bytes backing this data array.
    ///
    /// This is in turn used by [`ByteArrayViewMut::coerce_mut`] to produce a typed array
    fn view_mut(&'transient mut self) -> Result<&'transient mut Bytes, ArrayRetrievalError>;

    /// Reinterpret decoded bytes as a mutable slice in native byte order.
    ///
    /// The returned slice borrows the input buffer. `T` must be [`Pod`] so that
    /// all bit patterns are valid and writing values cannot introduce padding.
    /// Empty buffers return an empty slice.
    ///
    /// # Errors
    /// Returns [`ArrayRetrievalError::DataTypeSizeMismatch`] if the buffer is
    /// not aligned for `T` or its length is not a multiple of the element size.
    /// A nonempty buffer cannot be cast to a zero-sized type.
    ///
    /// The returned reference cannot outlive its buffer:
    ///
    /// ```compile_fail
    /// use mzdata_bindata::{ByteArrayViewMut, DataArray};
    ///
    /// fn escape_buffer() -> &'static mut [u8] {
    ///     let mut bytes = vec![1u8, 2, 3];
    ///     <DataArray as ByteArrayViewMut<'static, 'static>>::coerce_from_mut::<u8>(
    ///         &mut bytes,
    ///     ).unwrap()
    /// }
    /// ```
    ///
    /// Types such as `bool` cannot represent arbitrary bytes:
    ///
    /// ```compile_fail
    /// use mzdata_bindata::{ByteArrayViewMut, DataArray};
    ///
    /// let mut bytes = [0u8];
    /// let _ = <DataArray as ByteArrayViewMut>::coerce_from_mut::<bool>(&mut bytes);
    /// ```
    fn coerce_from_mut<T: Pod>(
        buffer: &'transient mut [u8],
    ) -> Result<&'transient mut [T], ArrayRetrievalError> {
        if buffer.is_empty() {
            return Ok(&mut []);
        }
        Ok(bytemuck::try_cast_slice_mut(buffer)?)
    }

    fn coerce_mut<T: Pod>(
        &'lifespan mut self,
    ) -> Result<&'transient mut [T], ArrayRetrievalError> {
        let view = self.view_mut()?;
        #[cfg(target_endian = "big")]
        {
            log::error!("Mutable view of raw bytes on big endian system is only partially supported.")
        }
        Self::coerce_from_mut(view)
    }

    #[allow(unused)]
    /// Set the identifier referencing a [`DataProcessing`](crate::meta::DataProcessing)
    fn set_data_processing_reference(&mut self, data_processing_reference: Option<Box<str>>) {}
}

#[derive(Debug)]
pub struct DataSliceIter<'a, T: Pod> {
    buffer: Cow<'a, [u8]>,
    i: usize,
    _t: PhantomData<T>
}

impl<T: Pod> ExactSizeIterator for DataSliceIter<'_, T> {
    fn len(&self) -> usize {
        let z = mem::size_of::<T>();
        self.buffer.len() / z - self.i
    }
}

impl<'a, T: Pod> DataSliceIter<'a, T> {
    pub fn new(buffer: Cow<'a, [u8]>) -> Self {
        Self { buffer, i: 0, _t: PhantomData }
    }

    pub fn next_value(&mut self) -> Option<T> {
        let z = mem::size_of::<T>();
        let offset = z * self.i;
        if (offset + z) > self.buffer.len() {
            None
        } else {
            let data = &self.buffer[offset..offset + z];
            #[cfg(target_endian = "big")]
            {
                let mut data = data.to_vec();
                data.reverse();
                let val = bytemuck::pod_read_unaligned(&data);
                self.i += 1;
                Some(val)
            }
            #[cfg(not(target_endian = "big"))]
            {
                let val = bytemuck::pod_read_unaligned(data);
                self.i += 1;
                Some(val)
            }
        }
    }
}

impl<T: Pod> Iterator for DataSliceIter<'_, T> {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_value()
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        if mem::size_of::<T>() == 0 {
            return (0, None);
        }
        let n = self.len();
        (n, Some(n))
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn check_remaining<T: Pod + PartialEq + std::fmt::Debug>(values: &[T]) {
        let mut it = DataSliceIter::<T>::new(Cow::Borrowed(bytemuck::cast_slice(values)));
        for (i, value) in values.iter().enumerate() {
            let n = values.len() - i;
            assert_eq!(it.len(), n);
            assert_eq!(it.size_hint(), (n, Some(n)));
            assert_eq!(it.next(), Some(*value));
        }
        assert_eq!(it.len(), 0);
        assert_eq!(it.size_hint(), (0, Some(0)));
        assert_eq!(it.next(), None);
        assert_eq!(it.next(), None);
    }

    #[test]
    fn test_iterator_remaining() {
        check_remaining::<u8>(&[]);
        check_remaining(&[1u8]);
        check_remaining(&[1u8, 2, 3]);
        check_remaining(&[1i32, -2, 3]);
        check_remaining(&[1i64, -2, 3]);
        check_remaining(&[1.25f32, -2.5, 3.0]);
        check_remaining(&[1.25f64, -2.5, 3.0]);

        let units = DataSliceIter::<()>::new(Cow::Borrowed(&[]));
        assert_eq!(units.take(3).collect::<Vec<_>>(), vec![(); 3]);
    }
}
