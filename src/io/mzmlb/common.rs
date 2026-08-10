use hdf5;
use hdf5::types::{TypeDescriptor, IntSize, FloatSize};

use crate::spectrum::BinaryDataArrayType;

pub(crate) fn type_descr_to_binary_data_array_type(value: &TypeDescriptor) -> BinaryDataArrayType {
    match value {
        TypeDescriptor::Integer(z) => {
            match z {
                IntSize::U1 => BinaryDataArrayType::Unknown,
                IntSize::U2 => BinaryDataArrayType::Unknown,
                IntSize::U4 => BinaryDataArrayType::Int32,
                IntSize::U8 => BinaryDataArrayType::Int64,
            }
        },
        TypeDescriptor::Unsigned(z) => {
            match z {
                IntSize::U1 => BinaryDataArrayType::Unknown,
                IntSize::U2 => BinaryDataArrayType::Unknown,
                IntSize::U4 => BinaryDataArrayType::Int32,
                IntSize::U8 => BinaryDataArrayType::Int64,
            }
        },
        TypeDescriptor::Float(z) => {
            match z {
                FloatSize::U4 => BinaryDataArrayType::Float32,
                FloatSize::U8 => BinaryDataArrayType::Float64,
            }
        },
        TypeDescriptor::Boolean => BinaryDataArrayType::Unknown,
        TypeDescriptor::Enum(_) => BinaryDataArrayType::Unknown,
        TypeDescriptor::Compound(_) => BinaryDataArrayType::Unknown,
        TypeDescriptor::FixedArray(_, _) => BinaryDataArrayType::Unknown,
        TypeDescriptor::FixedAscii(_) => {
            BinaryDataArrayType::ASCII
        },
        TypeDescriptor::FixedUnicode(_) => {
            BinaryDataArrayType::ASCII
        },
        TypeDescriptor::VarLenArray(_) => todo!(),
        TypeDescriptor::VarLenAscii => BinaryDataArrayType::Unknown,
        TypeDescriptor::VarLenUnicode => BinaryDataArrayType::Unknown,
    }
}

pub(crate) fn binary_array_data_type_to_descriptor(value: &BinaryDataArrayType) -> TypeDescriptor {
    match value {
        BinaryDataArrayType::Unknown => TypeDescriptor::Unsigned(IntSize::U1),
        BinaryDataArrayType::Float64 => TypeDescriptor::Float(FloatSize::U8),
        BinaryDataArrayType::Float32 => TypeDescriptor::Float(FloatSize::U4),
        BinaryDataArrayType::Int64 => TypeDescriptor::Integer(IntSize::U8),
        BinaryDataArrayType::Int32 => TypeDescriptor::Integer(IntSize::U4),
        BinaryDataArrayType::ASCII => TypeDescriptor::Unsigned(IntSize::U1),
    }
}


pub(crate) fn h5_type_to_binary_array_data_type(value: hdf5::Datatype) -> BinaryDataArrayType {
    match value.size() {
        1 => Self::ASCII,
        4 => {
            if value.is::<i32>() {
                Self::Int32
            } else if value.is::<f32>() {
                Self::Float32
            } else {
                Self::Unknown
            }
        }
        8 => {
            if value.is::<i64>() {
                Self::Int64
            } else if value.is::<f64>() {
                Self::Float64
            } else {
                Self::Unknown
            }
        }
        _ => Self::Unknown,
    }
}
