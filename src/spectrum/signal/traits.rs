use std::slice;
use std::mem;
use std::borrow::Cow;

use bytemuck::Pod;
use num_traits::{AsPrimitive, Num};
use super::encodings::{ArrayRetrievalError, BinaryDataArrayType, Bytes};


pub trait ByteArrayView<'transient, 'lifespan: 'transient> {
    fn view(&'lifespan self) -> Result<Cow<'lifespan, [u8]>, ArrayRetrievalError>;

    fn coerce_from<T: Pod>(
        buffer: Cow<'transient, [u8]>,
    ) -> Result<Cow<'transient, [T]>, ArrayRetrievalError> {
        let n = buffer.len();
        let z = mem::size_of::<T>();
        if n % z != 0 {
            return Err(ArrayRetrievalError::DataTypeSizeMismatch);
        }
        let m = n / z;
        unsafe {
            Ok(Cow::Borrowed(slice::from_raw_parts(
                buffer.as_ptr() as *const T,
                m,
            )))
        }
    }

    fn coerce<T: Pod>(
        &'lifespan self,
    ) -> Result<Cow<'transient, [T]>, ArrayRetrievalError> {
        match self.view() {
            Ok(data) => Self::coerce_from(data),
            Err(err) => Err(err),
        }
    }

    fn convert<S: Num + Clone + AsPrimitive<D> + Pod, D: Num + Clone + Copy + 'static>(
        &'lifespan self,
    ) -> Result<Cow<'transient, [D]>, ArrayRetrievalError> {
        match self.coerce::<S>() {
            Ok(view) => {
                match view {
                    Cow::Borrowed(view) => {
                        Ok(Cow::Owned(view.into_iter().map(|a| a.as_()).collect()))
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

    fn dtype(&self) -> BinaryDataArrayType;

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
}

pub trait ByteArrayViewMut<'transient, 'lifespan: 'transient>:
    ByteArrayView<'transient, 'lifespan>
{
    fn view_mut(&'transient mut self) -> Result<&'transient mut Bytes, ArrayRetrievalError>;

    fn coerce_from_mut<T: Clone + Sized>(
        buffer: &mut [u8],
    ) -> Result<&'transient mut [T], ArrayRetrievalError> {
        let n = buffer.len();
        let z = mem::size_of::<T>();
        if n % z != 0 {
            return Err(ArrayRetrievalError::DataTypeSizeMismatch);
        }
        let m = n / z;
        unsafe { Ok(slice::from_raw_parts_mut(buffer.as_ptr() as *mut T, m)) }
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
