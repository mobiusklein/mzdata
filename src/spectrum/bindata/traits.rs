use std::marker::PhantomData;
use std::slice;
use std::mem;
use std::borrow::Cow;

use bytemuck::Pod;
use num_traits::{AsPrimitive, Num};
use crate::params::Unit;

use super::encodings::{ArrayRetrievalError, BinaryDataArrayType, Bytes};
use super::ArrayType;


pub trait ByteArrayView<'transient, 'lifespan: 'transient> {
    fn view(&'lifespan self) -> Result<Cow<'lifespan, [u8]>, ArrayRetrievalError>;

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
        return match buffer {
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
        };
    }

    fn coerce<T: Pod>(
        &'lifespan self,
    ) -> Result<Cow<'transient, [T]>, ArrayRetrievalError> {
        match self.view() {
            Ok(data) => Self::coerce_from(data),
            Err(err) => Err(err),
        }
    }

    /// Decode the array, then copy it to a new array, converting each element from type `D` to to type `S`
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

pub trait ByteArrayViewMut<'transient, 'lifespan: 'transient>:
    ByteArrayView<'transient, 'lifespan>
{
    fn view_mut(&'transient mut self) -> Result<&'transient mut Bytes, ArrayRetrievalError>;

    fn coerce_from_mut<T: Clone + Sized>(
        buffer: &mut [u8],
    ) -> Result<&'transient mut [T], ArrayRetrievalError> {
        let n = buffer.len();
        if n == 0 {
            return Ok(&mut [])
        }
        let z = mem::size_of::<T>();
        if n % z != 0 {
            return Err(ArrayRetrievalError::DataTypeSizeMismatch);
        }
        let m = n / z;
        unsafe { Ok(slice::from_raw_parts_mut(buffer.as_mut_ptr() as *mut T, m)) }
    }

    fn coerce_mut<T: Clone + Sized>(
        &'lifespan mut self,
    ) -> Result<&'transient mut [T], ArrayRetrievalError> {
        let view = match self.view_mut() {
            Ok(data) => data,
            Err(err) => return Err(err),
        };
        Self::coerce_from_mut(view)
    }
}

#[derive(Debug)]
pub struct DataSliceIter<'a, T: Pod> {
    buffer: Cow<'a, [u8]>,
    i: usize,
    _t: PhantomData<T>
}

impl<'a, T: Pod> ExactSizeIterator for DataSliceIter<'a, T> {
    fn len(&self) -> usize {
        let z = mem::size_of::<T>();
        (self.buffer.len() / z) as usize
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
            let val = bytemuck::from_bytes(&self.buffer[offset..offset + z]);
            self.i += 1;
            Some(*val)
        }
    }
}

impl<'a, T: Pod> Iterator for DataSliceIter<'a, T> {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_value()
    }
}