use bytemuck::{self, Pod};

use std::{
    fmt::Display,
    io,
    ops::{AddAssign, Mul},
};

use thiserror::{self, Error};
use num_traits::{Num, ToBytes};

#[cfg(feature = "numpress")]
use numpress;

use mzdata_param::{curie, Param, ControlledVocabulary, ParamCow, Unit, CURIE, ValueRef};

pub type Bytes = Vec<u8>;

/// Convert the values in `data` into an owned buffer of little-endian bytes
///
/// # Examples
/// ```
/// use mzdata_bindata::to_bytes;
///
/// let values: [f64; 3] = [1.0, 2.5, -3.25];
/// let bytes = to_bytes(&values);
/// assert_eq!(bytes.len(), 3 * 8);
/// assert_eq!(&bytes[0..8], &1.0f64.to_le_bytes());
/// ```
pub fn to_bytes<T: Pod + ToBytes>(data: &[T]) -> Bytes {
    let n = data.len();
    let mut buf = Vec::with_capacity(n * size_of::<T>());
    for v in data {
        buf.extend_from_slice(v.to_le_bytes().as_ref());
    }
    buf
}

/// Reinterpret a slice of a [`bytemuck::Pod`] type as raw bytes without copying or
/// converting endianness.
///
/// Unlike [`to_bytes`], this does not convert byte order: the returned bytes are in the
/// host's native byte order, which only matches [`to_bytes`]'s little-endian output on a
/// little-endian host.
///
/// # Examples
/// ```
/// use mzdata_bindata::as_bytes;
///
/// let values: [i32; 2] = [7, -1];
/// let bytes = as_bytes(&values);
/// assert_eq!(bytes.len(), 8);
/// ```
pub fn as_bytes<T: Pod>(data: &[T]) -> &[u8] {
    bytemuck::cast_slice(data)
}

/// Convert an owned [`Vec`] of a [`bytemuck::Pod`] type into an owned buffer of bytes.
///
/// Unlike [`to_bytes`], this does not convert byte order: the returned bytes are in the
/// host's native byte order, which only matches [`to_bytes`]'s little-endian output on a
/// little-endian host.
///
/// # Examples
/// ```
/// use mzdata_bindata::vec_as_bytes;
///
/// let values: Vec<f32> = vec![1.0, 2.0, 3.0];
/// let bytes = vec_as_bytes(values);
/// assert_eq!(bytes.len(), 3 * 4);
/// ```
pub fn vec_as_bytes<T: Pod>(data: Vec<T>) -> Bytes {
    let mut buf = Bytes::with_capacity(data.len() * std::mem::size_of::<T>());
    for val in data {
        buf.extend_from_slice(bytemuck::bytes_of(&val));
    }
    buf
}

mod byte_rotation {
    use super::*;

    pub fn transpose_bytes_into<T: Pod, const N: usize>(data: &[T], buffer: &mut Vec<u8>) {
        assert_eq!(core::mem::size_of::<T>(), N);
        let bytes = bytemuck::cast_slice::<T, [u8; N]>(data);
        buffer.clear();
        let delta = (data.len() * N).saturating_sub(buffer.capacity());
        if delta > 0 {
            buffer.reserve(delta);
        }

        #[cfg(target_endian = "little")]
        {
            for i in 0..N {
                buffer.extend(bytes.iter().map(|b| b[i]))
            }
        }
        #[cfg(target_endian = "big")]
        {
            for i in (0..N).rev() {
                buffer.extend(bytes.iter().map(|b| b[i]))
            }
        }
    }

    pub fn transpose_bytes<T: Pod, const N: usize>(data: &[T]) -> Bytes {
        let mut result = Bytes::with_capacity(data.len() * N);
        transpose_bytes_into::<T, N>(data, &mut result);
        result
    }

    pub fn transpose_4bytes<T: Pod>(data: &[T]) -> Bytes {
        assert_eq!(std::mem::size_of::<T>(), 4);
        transpose_bytes::<_, 4>(data)
    }

    pub fn transpose_8bytes<T: Pod>(data: &[T]) -> Bytes {
        assert_eq!(std::mem::size_of::<T>(), 8);
        transpose_bytes::<_, 8>(data)
    }

    pub fn transpose_i32(data: &[i32]) -> Bytes {
        transpose_4bytes(data)
    }

    pub fn transpose_f32(data: &[f32]) -> Bytes {
        transpose_4bytes(data)
    }

    pub fn transpose_i64(data: &[i64]) -> Bytes {
        transpose_8bytes(data)
    }

    pub fn transpose_f64(data: &[f64]) -> Bytes {
        transpose_8bytes(data)
    }

    pub fn reverse_transpose_bytes_into<const N: usize>(data: &[u8], buffer: &mut Vec<u8>) {
        let rem = data.len() % N;
        assert_eq!(rem, 0);
        let n_entries = data.len() / N;
        buffer.clear();
        buffer.resize(data.len(), 0);

        #[cfg(target_endian = "little")]
        {
            for (i, band) in data.chunks_exact(n_entries).enumerate() {
                for (j, byte) in band.iter().copied().enumerate() {
                    bytemuck::cast_slice_mut::<_, [u8; N]>(buffer)[j][i] = byte;
                }
            }
        }
        #[cfg(target_endian = "big")]
        {
            for (i, band) in data.chunks_exact(n_entries).enumerate() {
                for (j, byte) in band.iter().copied().enumerate() {
                    bytemuck::cast_slice_mut::<_, [u8; N]>(buffer)[j][(N - 1) - i] = byte;
                }
            }
        }
    }

    pub fn reverse_transpose_bytes<const N: usize>(data: &[u8]) -> Bytes {
        let mut result: Bytes = vec![0; data.len()];
        reverse_transpose_bytes_into::<N>(data, &mut result);
        result
    }

    pub fn reverse_transpose_4bytes<T: Pod>(data: &[u8]) -> Bytes {
        assert_eq!(std::mem::size_of::<T>(), 4);
        reverse_transpose_bytes::<4>(data)
    }

    pub fn reverse_transpose_8bytes<T: Pod>(data: &[u8]) -> Bytes {
        assert_eq!(std::mem::size_of::<T>(), 8);
        reverse_transpose_bytes::<8>(data)
    }

    pub fn reverse_transpose_i32(data: &[u8]) -> Vec<u8> {
        reverse_transpose_4bytes::<i32>(data)
    }

    pub fn reverse_transpose_f32(data: &[u8]) -> Vec<u8> {
        reverse_transpose_4bytes::<f32>(data)
    }

    pub fn reverse_transpose_i64(data: &[u8]) -> Vec<u8> {
        reverse_transpose_8bytes::<i64>(data)
    }

    pub fn reverse_transpose_f64(data: &[u8]) -> Vec<u8> {
        reverse_transpose_8bytes::<f64>(data)
    }
}

mod dictionary_encoding {
    use super::*;
    use io::prelude::*;
    use num_traits::ops::bytes::{FromBytes, ToBytes};
    use std::{
        borrow::Cow,
        collections::{HashMap, HashSet},
        hash::Hash,
        io::BufWriter,
    };

    trait DictValue<const W: usize>:
        Pod + ToBytes<Bytes = [u8; W]> + Hash + Eq + Ord + FromBytes<Bytes = [u8; W]>
    {
    }

    macro_rules! impl_dict_value {
        ($val:ty, $size:literal) => {
            impl DictValue<$size> for $val {}
        };
    }

    impl_dict_value!(u8, 1);
    impl_dict_value!(u16, 2);
    impl_dict_value!(u32, 4);
    impl_dict_value!(u64, 8);

    trait DictIndex<const W: usize>:
        Pod + ToBytes<Bytes = [u8; W]> + FromBytes<Bytes = [u8; W]>
    {
        fn from_usize(index: usize) -> Self;
        fn to_usize(&self) -> usize;
    }

    macro_rules! impl_dict_index {
        ($idx:ty, $size:literal) => {
            impl DictIndex<$size> for $idx {
                fn from_usize(index: usize) -> Self {
                    index as Self
                }

                fn to_usize(&self) -> usize {
                    *self as usize
                }
            }
        };
    }

    impl_dict_index!(u8, 1);
    impl_dict_index!(u16, 2);
    impl_dict_index!(u32, 4);
    impl_dict_index!(u64, 8);

    #[derive(Default, Debug)]
    pub struct DictionaryEncoder {
        shuffle: bool,
        buffer: Vec<u8>,
    }

    impl DictionaryEncoder {
        pub fn new(shuffle: bool) -> Self {
            Self {
                shuffle,
                buffer: Vec::new(),
            }
        }

        fn build_value_map<T: Pod, const W1: usize, V: DictValue<W1>>(
            &self,
            data: &[T],
        ) -> (Vec<V>, HashMap<V, usize>) {
            debug_assert_eq!(core::mem::size_of::<T>(), core::mem::size_of::<V>());
            debug_assert_eq!(core::mem::size_of::<T>(), W1);
            let mut value_codes = HashSet::new();
            for v in data {
                let k: V = *bytemuck::from_bytes(bytemuck::bytes_of(v));
                value_codes.insert(k);
            }

            let mut value_codes: Vec<_> = value_codes.into_iter().collect();
            value_codes.sort();

            let byte_map: HashMap<V, usize> = value_codes
                .iter()
                .enumerate()
                .map(|(i, k)| (*k, i))
                .collect();
            (value_codes, byte_map)
        }

        fn create_writer<
            T: Pod,
            const W1: usize,
            V: Pod + Ord + Hash + ToBytes<Bytes = [u8; W1]> + Eq,
            const W2: usize,
            K: DictIndex<W2>,
        >(
            &self,
            data: &[T],
            value_codes: &[V],
        ) -> BufWriter<Vec<u8>> {
            let data_offset = 16 + std::mem::size_of_val(value_codes);
            let dict_buffer: Vec<u8> =
                Vec::with_capacity(data_offset + data.len() * core::mem::size_of::<K>());
            BufWriter::new(dict_buffer)
        }

        fn encode_dict_indices<
            T: Pod,
            const W1: usize,
            V: Pod + Ord + Hash + ToBytes<Bytes = [u8; W1]> + Eq,
            const W2: usize,
            K: DictIndex<W2>,
        >(
            &mut self,
            data: &[T],
            value_codes: &[V],
            byte_map: HashMap<V, usize>,
        ) -> io::Result<Vec<u8>> {
            let data_offset = 16 + std::mem::size_of_val(value_codes);
            let mut writer = self.create_writer::<T, W1, V, W2, K>(data, value_codes);
            writer.write_all(&(data_offset as u64).to_le_bytes())?;
            writer.write_all(&(value_codes.len() as u64).to_le_bytes())?;

            if self.shuffle {
                // This isn't endian-correct yet, see note in `transpose_bytes_into`
                byte_rotation::transpose_bytes_into::<V, W1>(value_codes, &mut self.buffer);
                writer.write_all(&self.buffer)?;
                // for v in self.buffer.chunks(W1) {
                //     writer.write_all(v)?;
                // }
            } else {
                for v in value_codes.iter() {
                    let bts = v.to_le_bytes();
                    writer.write_all(&bts)?;
                }
            }

            if self.shuffle {
                let mut buf = Vec::with_capacity(data.len());
                for v in data {
                    let i = *byte_map
                        .get(bytemuck::from_bytes(bytemuck::bytes_of(v)))
                        .unwrap();
                    let ik: K = K::from_usize(i);
                    buf.push(ik);
                }
                // This isn't endian-correct
                byte_rotation::transpose_bytes_into::<K, W2>(&buf, &mut self.buffer);
                writer.write_all(&self.buffer)?;
            } else {
                for v in data {
                    let i = *byte_map
                        .get(bytemuck::from_bytes(bytemuck::bytes_of(v)))
                        .unwrap();
                    let ik: K = K::from_usize(i);
                    writer.write_all(&ik.to_le_bytes())?;
                }
            }

            writer.flush()?;
            let val = writer.into_inner().unwrap();
            Ok(val)
        }

        fn encode_values<T: Pod, const W1: usize, V: DictValue<W1>>(
            &mut self,
            data: &[T],
        ) -> Result<Vec<u8>, io::Error> {
            let (value_codes, byte_map) = self.build_value_map(data);
            let n_value_codes = value_codes.len();

            if n_value_codes <= 2usize.pow(8) {
                self.encode_dict_indices::<T, W1, V, 1, u8>(data, &value_codes, byte_map)
            } else if n_value_codes <= 2usize.pow(16) {
                self.encode_dict_indices::<T, W1, V, 2, u16>(data, &value_codes, byte_map)
            } else if n_value_codes <= 2usize.pow(32) {
                self.encode_dict_indices::<T, W1, V, 4, u32>(data, &value_codes, byte_map)
            } else if n_value_codes <= 2usize.pow(64) {
                self.encode_dict_indices::<T, W1, V, 8, u64>(data, &value_codes, byte_map)
            } else {
                Err(io::Error::new(
                    io::ErrorKind::Unsupported,
                    "Cannot encode a dictionary with more than 2 ** 64 values",
                ))
            }
        }

        pub fn encode<T: Pod>(&mut self, data: &[T]) -> io::Result<Bytes> {
            if data.is_empty() {
                return Ok(Vec::new());
            }
            let z_val = core::mem::size_of::<T>();
            if z_val <= 1 {
                self.encode_values::<T, 1, u8>(data)
            } else if z_val <= 2 {
                self.encode_values::<T, 2, u16>(data)
            } else if z_val <= 4 {
                self.encode_values::<T, 4, u32>(data)
            } else if z_val <= 8 {
                self.encode_values::<T, 8, u64>(data)
            } else {
                Err(io::Error::new(
                    io::ErrorKind::Unsupported,
                    "Cannot encode a dictionary with more than 2 ** 64 keys",
                ))
            }
        }
    }

    #[derive(Default, Debug)]
    pub struct DictionaryDecoder {
        shuffle: bool,
        buffer: Vec<u8>,
    }

    impl DictionaryDecoder {
        pub fn new(shuffle: bool) -> Self {
            Self {
                shuffle,
                buffer: Default::default(),
            }
        }

        fn make_reader<'a>(&self, buffer: &'a [u8]) -> io::BufReader<&'a [u8]> {
            io::BufReader::new(buffer)
        }

        fn decode_value_buffer<T: Pod, const W1: usize, V: DictValue<W1>>(
            &mut self,
            buffer: &[u8],
            n_values: usize,
        ) -> Vec<T> {
            macro_rules! decode_chunk {
                ($chunk:ident) => {{
                    let chunk_a: [u8; W1] = $chunk.try_into().unwrap();
                    let val = V::from_le_bytes(&chunk_a);
                    let val: T = *bytemuck::from_bytes(bytemuck::bytes_of(&val));
                    val
                }};
            }
            let mut value_buffer = Vec::with_capacity(n_values);
            if self.shuffle {
                let blocks = match core::mem::size_of::<T>() {
                    1 => Cow::Borrowed(buffer),
                    2 => {
                        byte_rotation::reverse_transpose_bytes_into::<2>(buffer, &mut self.buffer);
                        Cow::Borrowed(self.buffer.as_slice())
                    }
                    4 => {
                        byte_rotation::reverse_transpose_bytes_into::<4>(buffer, &mut self.buffer);
                        Cow::Borrowed(self.buffer.as_slice())
                    }
                    8 => {
                        byte_rotation::reverse_transpose_bytes_into::<8>(buffer, &mut self.buffer);
                        Cow::Borrowed(self.buffer.as_slice())
                    }
                    x => {
                        panic!("Unsupported size {x}");
                    }
                };
                for chunk in blocks.chunks_exact(W1) {
                    let val = decode_chunk!(chunk);
                    value_buffer.push(val);
                }
            } else {
                for chunk in buffer.chunks_exact(W1) {
                    let val = decode_chunk!(chunk);
                    value_buffer.push(val);
                }
            };
            value_buffer
        }

        fn decode_index_buffer<T: Pod, const W2: usize, K: DictIndex<W2>>(
            &mut self,
            value_codes: &[T],
            index_buffer: &[u8],
        ) -> Vec<T> {
            let mut result = Vec::with_capacity(index_buffer.len() / W2);
            if self.shuffle {
                byte_rotation::reverse_transpose_bytes_into::<W2>(index_buffer, &mut self.buffer);
                for chunk in self.buffer.chunks_exact(W2) {
                    let b: [u8; W2] = chunk.try_into().unwrap();
                    let k: usize = K::from_le_bytes(&b).to_usize();
                    result.push(value_codes[k])
                }
            } else {
                for chunk in index_buffer.chunks_exact(W2) {
                    let b: [u8; W2] = chunk.try_into().unwrap();
                    let k: usize = K::from_le_bytes(&b).to_usize();
                    result.push(value_codes[k])
                }
            }
            result
        }

        pub fn decode<T: Pod>(&mut self, buffer: &[u8]) -> io::Result<Vec<T>> {
            if buffer.is_empty() {
                return Ok(Vec::new());
            }
            let mut reader = self.make_reader(buffer);
            let mut z_buf = [0u8; 8];
            reader.read_exact(&mut z_buf)?;
            let data_offset = u64::from_le_bytes(z_buf);
            if data_offset == 0 {
                return Ok(Vec::new());
            }
            let mut z_buf = [0u8; 8];
            reader.read_exact(&mut z_buf)?;
            let n_value_codes = u64::from_le_bytes(z_buf);
            if n_value_codes == 0 {
                return Ok(Vec::new());
            }
            let value_buffer = &buffer[16..(data_offset as usize)];
            let value_width = (data_offset - 16) / n_value_codes;
            let index_buffer = &buffer[data_offset as usize..];

            let n_value_codes = n_value_codes as usize;

            macro_rules! decode_indices {
                ($values:ident) => {
                    if n_value_codes <= 2usize.pow(8) {
                        self.decode_index_buffer::<T, 1, u8>(&$values, index_buffer)
                    } else if n_value_codes <= 2usize.pow(16) {
                        self.decode_index_buffer::<T, 2, u16>(&$values, index_buffer)
                    } else if n_value_codes <= 2usize.pow(32) {
                        self.decode_index_buffer::<T, 4, u32>(&$values, index_buffer)
                    } else if n_value_codes <= 2usize.pow(64) {
                        self.decode_index_buffer::<T, 8, u64>(&$values, index_buffer)
                    } else {
                        return Err(io::Error::new(
                            io::ErrorKind::Unsupported,
                            "Cannot decode a dictionary with more than 2 ** 64 indices",
                        ));
                    }
                };
            }

            let values = if value_width <= 1 {
                let values = self.decode_value_buffer::<T, 1, u8>(value_buffer, n_value_codes);
                decode_indices!(values)
            } else if value_width <= 2 {
                let values = self.decode_value_buffer::<T, 2, u16>(value_buffer, n_value_codes);
                decode_indices!(values)
            } else if value_width <= 4 {
                let values = self.decode_value_buffer::<T, 4, u32>(value_buffer, n_value_codes);
                decode_indices!(values)
            } else if value_width <= 8 {
                let values = self.decode_value_buffer::<T, 8, u64>(value_buffer, n_value_codes);
                decode_indices!(values)
            } else {
                return Err(io::Error::new(
                    io::ErrorKind::Unsupported,
                    "Cannot decode dictionary with value byte width greater than 8",
                ));
            };

            Ok(values)
        }
    }

    pub fn dictionary_encoding<T: Pod>(data: &[T]) -> Result<Vec<u8>, io::Error> {
        let mut encoder = DictionaryEncoder::new(true);
        encoder.encode(data)
    }

    pub fn dictionary_decoding<T: Pod>(buffer: &[u8]) -> io::Result<Vec<T>> {
        let mut decoder = DictionaryDecoder::new(true);
        decoder.decode(buffer)
    }
}

#[allow(unused)]
pub use byte_rotation::*;

#[allow(unused)]
pub use dictionary_encoding::{dictionary_decoding, dictionary_encoding};

/// The kinds of data arrays found in mass spectrometry data files governed
/// by the PSI-MS controlled vocabulary.
#[derive(Debug, Clone, PartialEq, Hash, Eq, PartialOrd, Ord, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ArrayType {
    /// No information is given about this array, not even a name
    #[default]
    Unknown,
    /// The m/z array, the coordinate values of a mass spectrum. The default unit is [`Unit::MZ`].
    MZArray,
    /// The intensity or abundance of a signal at a given point.
    /// The default unit is [`Unit::DetectorCounts`].
    IntensityArray,
    /// The deconvolved charge state of an ion
    ChargeArray,
    /// The estimated signal to noise ratio of a point observation
    SignalToNoiseArray,
    /// The time coordinate an observation occured such as for a feature, trace, or chromatgam.
    /// The default unit is [`Unit::Minute`]
    TimeArray,
    /// The wavelength coordinate of an observation in a photospectrometry measurement.
    /// The default unit is [`Unit::Nanometer`].
    WavelengthArray,

    /// A basic ion mobility array that contains no additional information about the unit
    /// or level of processing or aggregation. **Note: It has no implicit unit information!**
    IonMobilityArray,

    /// The mean ion mobility measurement of an observation. This is the centroid of an ion mobility
    /// measurement, distinct from a profile array view. **Note: It has no implicit unit information!**
    MeanIonMobilityArray,
    /// The mean drift time of an observation. This is the centroid of an ion mobility measurement
    /// expressed in terms of drift time. The default unit is [`Unit::Millisecond`].
    MeanDriftTimeArray,
    /// The mean drift time of an observation. This is the centroid of an ion mobility measurement
    /// expressed in terms of $V*s/cm^2$. The default unit is [`Unit::VoltSecondPerSquareCentimeter`].
    MeanInverseReducedIonMobilityArray,

    /// The raw ion mobility measurement of an observation. This is unprocessed and may be profile
    /// or centroid ion mobility measurements.
    RawIonMobilityArray,
    /// The raw drift time of an observation. This is may be a profile or centroid representation
    /// of an ion mobility measurement expressed in terms of drift time. The default unit is [`Unit::Millisecond`].
    RawDriftTimeArray,
    /// The raw drift time of an observation. This is may be a profile or centroid representation
    /// of an ion mobility measurement expressed in terms of $V*s/cm^2$. The default unit is [`Unit::VoltSecondPerSquareCentimeter`].
    RawInverseReducedIonMobilityArray,

    /// The mean ion mobility measurement of an observation after deconvolution of cross-cutting properties like
    /// charge state, adduct, or some other relevant property. This is the centroid of an ion mobility
    /// measurement, distinct from a profile array view. **Note: It has no implicit unit information!**
    DeconvolutedIonMobilityArray,
    /// The mean ion mobility drift time measurement of an observation expressed in drift time after deconvolution
    /// of cross-cutting properties like charge state, adduct, or some other relevant property. This is the centroid of
    /// an ion mobility measurement. The default unit is [`Unit::Millisecond`].
    DeconvolutedDriftTimeArray,
    /// The mean ion mobility drift time measurement of an observation expressed in $V*s/cm^2$ after deconvolution
    /// of cross-cutting properties like charge state, adduct, or some other relevant property. This is the centroid of
    /// an ion mobility measurement. The default unit is [`Unit::Millisecond`].
    DeconvolutedInverseReducedIonMobilityArray,

    /// Array of m/z values representing the lower bound m/z of the quadrupole position at each point in the spectrum.
    /// The default unit is [`Unit::MZ`]
    ScanningQuadrupolePositionLowerBoundMZ,
    /// Array of m/z values representing the upper bound m/z of the quadrupole position at each point in the spectrum.
    /// The default unit is [`Unit::MZ`]
    ScanningQuadrupolePositionUpperBoundMZ,

    /// A data array whose values are index values for another collection or array.
    IndexArray,

    /// A data array of signal baseline values (the signal in the absence of analytes).
    BaselineArray,
    /// A data array of resolution values, which influence the accuracy of measurements.
    ResolutionArray,
    /// A data array of pressure measurements.
    PressureArray,
    /// A data array of temperature measurements.
    TemperatureArray,
    /// A data array of flow rate measurements.
    FlowRateArray,

    /// A data array that is not covered by any of the existing types, but given a distinct name.
    /// This is used for user- or application-defined arrays.
    NonStandardDataArray {
        name: Box<String>,
    },
}

impl Display for ArrayType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl ArrayType {
    /// Get the data type that the array is compatible with in the
    /// `mzdata` type expectations.
    ///
    /// By default, the m/z array is encoded using `Float64`,
    /// the charge state array is encoded using `Int32`, and
    /// all other arrays are encoded using `Float32`.
    ///
    /// # Examples
    /// ```
    /// use mzdata_bindata::{ArrayType, BinaryDataArrayType};
    ///
    /// assert_eq!(ArrayType::MZArray.preferred_dtype(), BinaryDataArrayType::Float64);
    /// assert_eq!(ArrayType::IntensityArray.preferred_dtype(), BinaryDataArrayType::Float32);
    /// assert_eq!(ArrayType::ChargeArray.preferred_dtype(), BinaryDataArrayType::Int32);
    /// ```
    pub const fn preferred_dtype(&self) -> BinaryDataArrayType {
        match self {
            ArrayType::MZArray => BinaryDataArrayType::Float64,
            ArrayType::IntensityArray => BinaryDataArrayType::Float32,
            ArrayType::ChargeArray => BinaryDataArrayType::Int32,
            ArrayType::IndexArray => BinaryDataArrayType::Int64,
            _ => BinaryDataArrayType::Float32,
        }
    }

    /// Convert an ion mobility array to its mean variant
    ///
    /// # Examples
    /// ```
    /// use mzdata_bindata::ArrayType;
    ///
    /// assert_eq!(
    ///     ArrayType::RawDriftTimeArray.as_mean_ion_mobility(),
    ///     Some(ArrayType::MeanDriftTimeArray)
    /// );
    /// assert_eq!(ArrayType::MZArray.as_mean_ion_mobility(), None);
    /// ```
    pub const fn as_mean_ion_mobility(&self) -> Option<ArrayType> {
        Some(match self {
            Self::RawDriftTimeArray
            | Self::DeconvolutedDriftTimeArray
            | Self::MeanDriftTimeArray => Self::MeanDriftTimeArray,
            Self::RawInverseReducedIonMobilityArray
            | Self::DeconvolutedInverseReducedIonMobilityArray
            | Self::MeanInverseReducedIonMobilityArray => Self::MeanInverseReducedIonMobilityArray,
            Self::RawIonMobilityArray
            | Self::DeconvolutedIonMobilityArray
            | Self::MeanIonMobilityArray => Self::MeanIonMobilityArray,
            _ => return None,
        })
    }

    /// Convert an ion mobility array to its raw variant
    pub const fn as_raw_ion_mobility(&self) -> Option<ArrayType> {
        Some(match self {
            Self::RawDriftTimeArray
            | Self::DeconvolutedDriftTimeArray
            | Self::MeanDriftTimeArray => Self::RawDriftTimeArray,
            Self::RawInverseReducedIonMobilityArray
            | Self::DeconvolutedInverseReducedIonMobilityArray
            | Self::MeanInverseReducedIonMobilityArray => Self::RawInverseReducedIonMobilityArray,
            Self::RawIonMobilityArray
            | Self::DeconvolutedIonMobilityArray
            | Self::MeanIonMobilityArray => Self::RawIonMobilityArray,
            _ => return None,
        })
    }

    /// Convert an ion mobility array to its deconvoluted variant
    pub const fn as_deconvoluted_ion_mobility(&self) -> Option<ArrayType> {
        Some(match self {
            Self::RawDriftTimeArray
            | Self::DeconvolutedDriftTimeArray
            | Self::MeanDriftTimeArray => Self::DeconvolutedDriftTimeArray,
            Self::RawInverseReducedIonMobilityArray
            | Self::DeconvolutedInverseReducedIonMobilityArray
            | Self::MeanInverseReducedIonMobilityArray => {
                Self::DeconvolutedInverseReducedIonMobilityArray
            }
            Self::RawIonMobilityArray
            | Self::DeconvolutedIonMobilityArray
            | Self::MeanIonMobilityArray => Self::DeconvolutedIonMobilityArray,
            _ => return None,
        })
    }

    /// Create a [`ArrayType::NonStandardDataArray`] with the provided name.
    pub fn nonstandard<S: ToString>(name: S) -> ArrayType {
        ArrayType::NonStandardDataArray {
            name: name.to_string().into(),
        }
    }

    /// The preferred unit for this array type.
    pub const fn preferred_unit(&self) -> Option<Unit> {
        let unit = match self {
            Self::MZArray | Self::ScanningQuadrupolePositionLowerBoundMZ | Self::ScanningQuadrupolePositionUpperBoundMZ => Unit::MZ,
            Self::IntensityArray => Unit::DetectorCounts,
            Self::WavelengthArray => Unit::Nanometer,
            Self::TimeArray | Self::FlowRateArray => Unit::Minute,
            Self::PressureArray => Unit::Pascal,
            Self::RawDriftTimeArray | Self::MeanDriftTimeArray | Self::DeconvolutedDriftTimeArray => Unit::Millisecond,
            Self::RawInverseReducedIonMobilityArray | Self::MeanInverseReducedIonMobilityArray | Self::DeconvolutedInverseReducedIonMobilityArray => Unit::VoltSecondPerSquareCentimeter,
            _ => return None
        };
        Some(unit)
    }

    /// Test if the the array describes an ion mobility quantity.
    ///
    /// # Examples
    /// ```
    /// use mzdata_bindata::ArrayType;
    ///
    /// assert!(!ArrayType::MZArray.is_ion_mobility());
    /// assert!(ArrayType::MeanInverseReducedIonMobilityArray.is_ion_mobility());
    /// ```
    pub const fn is_ion_mobility(&self) -> bool {
        matches!(
            self,
            Self::IonMobilityArray
                | Self::MeanIonMobilityArray
                | Self::MeanDriftTimeArray
                | Self::MeanInverseReducedIonMobilityArray
                | Self::DeconvolutedIonMobilityArray
                | Self::DeconvolutedDriftTimeArray
                | Self::DeconvolutedInverseReducedIonMobilityArray
                | Self::RawIonMobilityArray
                | Self::RawDriftTimeArray
                | Self::RawInverseReducedIonMobilityArray,
        )
    }

    /// Create a [`Param`] for this array type.
    ///
    /// If a unit is provided, that unit will be specified, otherwise a default unit may
    /// be used instead.
    ///
    /// Prefer this to [`Self::as_param_const`] when you do not need to worry about const-ness
    /// and the
    pub fn as_param(&self, unit: Option<Unit>) -> Param {
        const CV: ControlledVocabulary = ControlledVocabulary::MS;
        let unit = unit.or_else(|| self.preferred_unit());

        if let Self::NonStandardDataArray { name } = self {
            let mut p = CV.param_val(
                1000786,
                "non-standard data array",
                name.to_string()
            );
            p.unit = unit.unwrap_or_default();
            p
        }
        else {
            match unit {
                Some(unit) => self.as_param_with_unit_const(unit).into(),
                None => self.as_param_const().into()
            }
        }

        // match self {
        //     ArrayType::MZArray => CV
        //         .const_param_ident_unit("m/z array", 1000514, unit.unwrap_or_default())
        //         .into(),
        //     ArrayType::IntensityArray => CV
        //         .const_param_ident_unit(
        //             "intensity array",
        //             1000515,
        //             unit.unwrap_or_default(),
        //         )
        //         .into(),
        //     ArrayType::ChargeArray => CV.const_param_ident("charge array", 1000516).into(),
        //     ArrayType::TimeArray => CV
        //         .const_param_ident_unit("time array", 1000595, unit.unwrap_or_default())
        //         .into(),
        //     ArrayType::WavelengthArray => CV
        //         .const_param_ident_unit("wavelength array", 1000617, unit.unwrap_or_default())
        //         .into(),
        //     ArrayType::SignalToNoiseArray => CV
        //         .const_param_ident("signal to noise array", 1000517)
        //         .into(),
        //     ArrayType::IonMobilityArray => CV
        //         .const_param_ident_unit(
        //             "ion mobility array", 1002893, unit.unwrap_or_default())
        //         .into(),

        //     ArrayType::RawDriftTimeArray => CV
        //         .const_param_ident_unit(
        //             "raw ion mobility drift time array",
        //             1003153,
        //             unit.unwrap_or_default(),
        //         )
        //         .into(),
        //     ArrayType::RawInverseReducedIonMobilityArray => CV
        //         .const_param_ident_unit(
        //             "raw inverse reduced ion mobility array",
        //             1003008,
        //             unit.unwrap_or_default(),
        //         )
        //         .into(),
        //     ArrayType::RawIonMobilityArray => CV
        //         .const_param_ident_unit("raw ion mobility array", 1003007, unit.unwrap_or_default())
        //         .into(),

        //     ArrayType::MeanIonMobilityArray => CV
        //         .const_param_ident_unit(
        //             "mean ion mobility array",
        //             1002816,
        //             unit.unwrap_or_default(),
        //         )
        //         .into(),
        //     ArrayType::MeanDriftTimeArray => CV
        //         .const_param_ident_unit(
        //             "mean ion mobility drift time array",
        //             1002477,
        //             unit.unwrap_or_default(),
        //         )
        //         .into(),
        //     ArrayType::MeanInverseReducedIonMobilityArray => CV
        //         .const_param_ident_unit(
        //             "mean inverse reduced ion mobility array",
        //             1003006,
        //             unit.unwrap_or_default(),
        //         )
        //         .into(),

        //     ArrayType::DeconvolutedIonMobilityArray => CV
        //         .const_param_ident_unit(
        //             "deconvoluted ion mobility array",
        //             1003154,
        //             unit.unwrap_or_default(),
        //         )
        //         .into(),
        //     ArrayType::DeconvolutedDriftTimeArray => CV
        //         .const_param_ident_unit(
        //             "deconvoluted ion mobility drift time array",
        //             1003156,
        //             unit.unwrap_or_default(),
        //         )
        //         .into(),
        //     ArrayType::DeconvolutedInverseReducedIonMobilityArray => CV
        //         .const_param_ident_unit(
        //             "deconvoluted inverse reduced ion mobility array",
        //             1003155,
        //             unit.unwrap_or_default(),
        //         )
        //         .into(),

        //     ArrayType::NonStandardDataArray { name } => {
        //         let mut p = CV.param_val(1000786, "non-standard data array", name.to_string());
        //         p.unit = unit.unwrap_or_default();
        //         p
        //     }
        //     ArrayType::BaselineArray => CV.const_param_ident("baseline array", 1002530).into(),
        //     ArrayType::ResolutionArray => CV.const_param_ident("resolution array", 1002529).into(),
        //     ArrayType::PressureArray => {
        //         let mut p = CV.const_param_ident("pressure array", 1000821);
        //         p.unit = unit.unwrap_or_default();
        //         p.into()
        //     }
        //     ArrayType::TemperatureArray => {
        //         let mut p = CV.const_param_ident("temperature array", 1000822);
        //         p.unit = unit.unwrap_or_default();
        //         p.into()
        //     }
        //     ArrayType::FlowRateArray => {
        //         let mut p = CV.const_param_ident("flow rate array", 1000820);
        //         p.unit = unit.unwrap_or_default();
        //         p.into()
        //     }
        //     ArrayType::ScanningQuadrupolePositionLowerBoundMZ => {
        //         let mut p = CV.const_param_ident("scanning quadrupole position lower bound m/z array", 1003157);
        //         p.unit = unit.unwrap_or_default();
        //         p.into()
        //     }
        //     ArrayType::ScanningQuadrupolePositionUpperBoundMZ => {
        //         let mut p = CV.const_param_ident("scanning quadrupole position upper bound m/z array", 1003158);
        //         p.unit = unit.unwrap_or_default();
        //         p.into()
        //     }
        //     ArrayType::IndexArray => {
        //         CV.const_param_ident("index array", 1003870).into()
        //     }
        //     _ => {
        //         panic!("Could not determine how to name for array {}", self);
        //     }
        // }
    }

    /// Create a [`ParamCow`] for this array type in a `const` context using the default unit.
    ///
    /// **NOTE**: If this is a [`Self::NonStandardDataArray`], this function *will* panic as
    /// this variant requires allocation.
    ///
    /// Prefer this method over [`Self::as_param_with_unit_const`] when the default unit is
    /// acceptable, and you are confident you are not dealing with a non-standard array.
    ///
    /// # Panics
    /// When `self` is either [`Self::NonStandardDataArray`] or an unimplemented variant.
    pub const fn as_param_const(&self) -> ParamCow<'static> {
        const CV: ControlledVocabulary = ControlledVocabulary::MS;

        match self {
            ArrayType::MZArray => CV.const_param_ident_unit("m/z array", 1000514, self.preferred_unit().unwrap()),
            ArrayType::IntensityArray => {
                CV.const_param_ident_unit("intensity array", 1000515, self.preferred_unit().unwrap())
            }
            ArrayType::ChargeArray => CV.const_param_ident("charge array", 1000516),
            ArrayType::TimeArray => CV.const_param_ident_unit("time array", 1000595, self.preferred_unit().unwrap()),
            ArrayType::WavelengthArray => {
                CV.const_param_ident_unit("wavelength array", 1000617, self.preferred_unit().unwrap())
            }
            ArrayType::SignalToNoiseArray => CV.const_param_ident("signal to noise array", 1000517),
            ArrayType::IonMobilityArray => CV.const_param_ident("ion mobility array", 1002893),
            ArrayType::RawIonMobilityArray => {
                CV.const_param_ident("raw ion mobility array", 1003007)
            }
            ArrayType::MeanIonMobilityArray => {
                CV.const_param_ident("mean ion mobility array", 1002816)
            }
            ArrayType::DeconvolutedIonMobilityArray => {
                CV.const_param_ident("deconvoluted ion mobility array", 1003154)
            }
            ArrayType::RawDriftTimeArray => CV.const_param_ident_unit(
                "raw ion mobility drift time array",
                1003153,
                self.preferred_unit().unwrap(),
            ),
            ArrayType::RawInverseReducedIonMobilityArray => CV.const_param_ident_unit(
                "raw inverse reduced ion mobility array",
                1003008,
                self.preferred_unit().unwrap(),
            ),

            ArrayType::MeanDriftTimeArray => CV.const_param_ident_unit(
                "mean ion mobility drift time array",
                1002477,
                self.preferred_unit().unwrap(),
            ),
            ArrayType::MeanInverseReducedIonMobilityArray => CV.const_param_ident_unit(
                "mean inverse reduced ion mobility array",
                1003006,
                self.preferred_unit().unwrap(),
            ),

            ArrayType::DeconvolutedDriftTimeArray => CV.const_param_ident_unit(
                "deconvoluted ion mobility drift time array",
                1003156,
                self.preferred_unit().unwrap(),
            ),
            ArrayType::DeconvolutedInverseReducedIonMobilityArray => CV.const_param_ident_unit(
                "deconvoluted inverse reduced ion mobility array",
                1003155,
                self.preferred_unit().unwrap(),
            ),

            ArrayType::NonStandardDataArray { name: _name } => {
                panic!(
                    "Cannot format NonStandardDataArray in a const context, please use `as_param`"
                );
            }
            ArrayType::BaselineArray => CV.const_param_ident("baseline array", 1002530),
            ArrayType::ResolutionArray => CV.const_param_ident("resolution array", 1002529),
            ArrayType::PressureArray => CV.const_param_ident_unit("pressure array", 1000821, self.preferred_unit().unwrap()),
            ArrayType::TemperatureArray => CV.const_param_ident("temperature array", 1000822),
            ArrayType::FlowRateArray => CV.const_param_ident_unit("flow rate array", 1000820, self.preferred_unit().unwrap()),
            ArrayType::ScanningQuadrupolePositionLowerBoundMZ => {
                let mut p = CV.const_param_ident("scanning quadrupole position lower bound m/z array", 1003157);
                p.unit = self.preferred_unit().unwrap();
                p
            }
            ArrayType::IndexArray => {
                CV.const_param_ident("index array", 1003870)
            }
            ArrayType::ScanningQuadrupolePositionUpperBoundMZ => {
                let mut p = CV.const_param_ident("scanning quadrupole position upper bound m/z array", 1003158);
                p.unit = self.preferred_unit().unwrap();
                p
            }
            _ => {
                panic!("Could not determine how to name for array");
            }
        }
    }

    /// Create a [`ParamCow`] for this array type in a `const` context using the provided unit.
    ///
    /// **NOTE**: If this is a [`Self::NonStandardDataArray`], this function *will* panic as
    /// this variant requires allocation.
    ///
    /// Prefer this method over [`Self::as_param_const`] when you want to control the unit
    /// directly rather than accepting the default unit and possibly updating it in the calling
    /// code.
    ///
    /// # Panics
    /// When `self` is either [`Self::NonStandardDataArray`] or an unimplemented variant.
    pub const fn as_param_with_unit_const(&self, unit: Unit) -> ParamCow<'static> {
        let mut inner = self.as_param_const();
        inner.unit = unit;
        inner
    }

    /// Translate from [`CURIE`] to an [`ArrayType`]
    pub fn from_accession(x: CURIE) -> Option<Self> {
        let tp = if x == Self::MZArray.as_param_const().curie().unwrap() {
            Self::MZArray
        } else if x == Self::IntensityArray.as_param_const().curie().unwrap() {
            Self::IntensityArray
        } else if x == Self::ChargeArray.as_param_const().curie().unwrap() {
            Self::ChargeArray
        } else if x == Self::SignalToNoiseArray.as_param_const().curie().unwrap() {
            Self::SignalToNoiseArray
        } else if x == Self::TimeArray.as_param_const().curie().unwrap() {
            Self::TimeArray
        } else if x == Self::WavelengthArray.as_param_const().curie().unwrap() {
            Self::WavelengthArray
        } else if x == Self::IonMobilityArray.as_param_const().curie().unwrap() {
            Self::IonMobilityArray
        } else if x == Self::MeanIonMobilityArray.as_param_const().curie().unwrap() {
            Self::MeanIonMobilityArray
        } else if x == Self::MeanDriftTimeArray.as_param_const().curie().unwrap() {
            Self::MeanDriftTimeArray
        } else if x
            == Self::MeanInverseReducedIonMobilityArray
                .as_param_const()
                .curie()
                .unwrap()
        {
            Self::MeanInverseReducedIonMobilityArray
        } else if x == Self::RawIonMobilityArray.as_param_const().curie().unwrap() {
            Self::RawIonMobilityArray
        } else if x == Self::RawDriftTimeArray.as_param_const().curie().unwrap() {
            Self::RawDriftTimeArray
        } else if x
            == Self::RawInverseReducedIonMobilityArray
                .as_param_const()
                .curie()
                .unwrap()
        {
            Self::RawInverseReducedIonMobilityArray
        } else if x
            == Self::DeconvolutedIonMobilityArray
                .as_param_const()
                .curie()
                .unwrap()
        {
            Self::DeconvolutedIonMobilityArray
        } else if x
            == Self::DeconvolutedDriftTimeArray
                .as_param_const()
                .curie()
                .unwrap()
        {
            Self::DeconvolutedDriftTimeArray
        } else if x
            == Self::DeconvolutedInverseReducedIonMobilityArray
                .as_param_const()
                .curie()
                .unwrap()
        {
            Self::DeconvolutedInverseReducedIonMobilityArray
        } else if x == Self::BaselineArray.as_param_const().curie().unwrap() {
            Self::BaselineArray
        } else if x == Self::ResolutionArray.as_param_const().curie().unwrap() {
            Self::ResolutionArray
        } else if x == Self::PressureArray.as_param_const().curie().unwrap() {
            Self::PressureArray
        } else if x == Self::TemperatureArray.as_param_const().curie().unwrap() {
            Self::TemperatureArray
        } else if x == Self::FlowRateArray.as_param_const().curie().unwrap() {
            Self::FlowRateArray
        } else if x == Self::ScanningQuadrupolePositionLowerBoundMZ.as_param_const().curie().unwrap() {
            Self::ScanningQuadrupolePositionLowerBoundMZ
        } else if x == Self::ScanningQuadrupolePositionUpperBoundMZ.as_param_const().curie().unwrap() {
            Self::ScanningQuadrupolePositionUpperBoundMZ
        } else if x == Self::IndexArray.as_param_const().curie().unwrap() {
            Self::IndexArray
        }
        else if x
            == (Self::NonStandardDataArray {
                name: "".to_string().into(),
            })
            .as_param(None)
            .curie()
            .unwrap()
        {
            Self::NonStandardDataArray {
                name: "".to_string().into(),
            }
        } else {
            return None;
        };
        Some(tp)
    }
}

/// The canonical primitive data types found in MS data file formats
/// supported by the PSI-MS controlled vocabulary
#[derive(Debug, Clone, Copy, PartialEq, Hash, Eq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum BinaryDataArrayType {
    /// Unspecified, opaque data
    #[default]
    Unknown,
    /// 64-bit floating point numbers
    Float64,
    /// 32-bit floating point numbers
    Float32,
    /// 64-bit signed integers
    Int64,
    /// 32-bit signed integers
    Int32,
    /// Null-terminated ASCII strings
    ASCII,
}

impl Display for BinaryDataArrayType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl BinaryDataArrayType {
    /// Get the size in bytes of a single value of this type
    ///
    /// # Examples
    /// ```rust
    /// use mzdata_bindata::BinaryDataArrayType;
    ///
    /// assert_eq!(BinaryDataArrayType::Float32.size_of(), 4);
    /// assert_eq!(BinaryDataArrayType::Float64.size_of(), 8);
    /// assert_eq!(BinaryDataArrayType::ASCII.size_of(), 1);
    /// ```
    pub const fn size_of(&self) -> usize {
        match self {
            BinaryDataArrayType::Unknown | BinaryDataArrayType::ASCII => 1,
            BinaryDataArrayType::Float32 | BinaryDataArrayType::Int32 => 4,
            BinaryDataArrayType::Float64 | BinaryDataArrayType::Int64 => 8,
        }
    }

    /// Convert the data type to a [`ParamCow`] if possible.
    pub const fn as_param_const(&self) -> Option<ParamCow<'static>> {
        let name = match self {
            BinaryDataArrayType::Unknown => return None,
            BinaryDataArrayType::Float64 => "64-bit float",
            BinaryDataArrayType::Float32 => "32-bit float",
            BinaryDataArrayType::Int64 => "64-bit integer",
            BinaryDataArrayType::Int32 => "32-bit integer",
            BinaryDataArrayType::ASCII => "null-terminated ASCII string",
        };
        if let Some(curie) = self.curie() {
            Some(ParamCow::const_new(
                name,
                ValueRef::Empty,
                Some(curie.accession),
                Some(curie.controlled_vocabulary),
                Unit::Unknown,
            ))
        } else {
            None
        }
    }

    /// Convert the data type to a [`CURIE`] if one exists for it
    pub const fn curie(&self) -> Option<CURIE> {
        match self {
            Self::Float32 => Some(curie!(MS:1000521)),
            Self::Float64 => Some(curie!(MS:1000523)),
            Self::Int32 => Some(curie!(MS:1000519)),
            Self::Int64 => Some(curie!(MS:1000522)),
            Self::ASCII => Some(curie!(MS:1001479)),
            _ => None,
        }
    }

    pub fn from_accession(accession: CURIE) -> Option<Self> {
        match accession {
            x if Some(x) == Self::Float32.curie() => Some(Self::Float32),
            x if Some(x) == Self::Float64.curie() => Some(Self::Float64),
            x if Some(x) == Self::Int32.curie() => Some(Self::Int32),
            x if Some(x) == Self::Int64.curie() => Some(Self::Int64),
            x if Some(x) == Self::ASCII.curie() => Some(Self::ASCII),
            _ => None,
        }
    }

    /// Byte-order swap the data stored in `data`, treating it as a sequence of
    /// fixed-width elements of `self`'s size.
    ///
    /// # Examples
    /// ```rust
    /// use mzdata_bindata::BinaryDataArrayType;
    ///
    /// let mut bytes = 42u32.to_le_bytes();
    /// BinaryDataArrayType::Int32.swap_bytes(&mut bytes).unwrap();
    /// assert_eq!(u32::from_be_bytes(bytes), 42);
    /// ```
    pub fn swap_bytes(&self, data: &mut [u8]) -> Result<(), ArrayRetrievalError> {
        let z = self.size_of();
        if !(data.len() % z == 0) {
            return Err(ArrayRetrievalError::DataTypeSizeMismatch)
        }
        match z {
            1 => {
                data.reverse();
            }
            4 => {
                data.as_chunks_mut::<4>().0.into_iter().for_each(|c| {
                    *c = u32::from_ne_bytes(*c).swap_bytes().to_ne_bytes();
                });
            }
            8 => {
                data.as_chunks_mut::<8>().0.into_iter().for_each(|c| {
                    *c = u64::from_ne_bytes(*c).swap_bytes().to_ne_bytes();
                });
                // data.chunks_exact_mut(4).for_each(|c| {
                //     c.swap(0, 7);
                //     c.swap(1, 6);
                //     c.swap(2, 5);
                //     c.swap(3, 4);
                // });
            }
            x => {
                data.chunks_exact_mut(x).for_each(|c| c.reverse());
            }
        }
        Ok(())
    }
}

/// The range of compression and encoding states that a raw byte buffer
/// might be in during different stages of decoding. Other than `Decoded`,
/// these states may or may not include intermediate base64 encoding.
#[derive(Debug, Clone, Copy, PartialEq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum BinaryCompressionType {
    /// The data are not compressed, but are base64 encoded and not ready for native use.
    #[default]
    NoCompression,

    /// A well known, widely available general-purpose compressor codec with many implementations. Always
    /// available, but several backends with differing speed and size tradeoffs are available under different
    /// features.
    Zlib,

    /// **Warning: Lossy** A combination simple linear model fit and byte packing codec that stores the model
    /// coefficients as float64 and the residuals as unsigned 16 or 32-bit integers. With an optimal fixed point, this will produce a median loss of accuracy below
    /// of about 2e-8 absolute units, but may vary widely from dataset to dataset. This codec **MUST** be used with sorted data only, and is only suitable for
    /// float64 data. It may be combined with `zlib` [`BinaryCompressionType::NumpressLinearZlib`] or `zstd` [`BinaryCompressionType::NumpressLinearZstd`] for even
    /// greater space savings. See details at <https://doi.org/10.1074/mcp.O114.037879>
    ///
    /// Requires feature `numpress`
    NumpressLinear,

    /// **Warning: Lossy** A lossy numerical approximation that log transforms and scales values to make them very
    /// close to 16-bit integers. This codec is lossier the larger the value is, and is only appropriate for count data like intensities. It can also combine with
    /// `zlib` [`BinaryCompressionType::NumpressSLOFZlib`] and `zstd` [`BinaryCompressionType::NumpressSLOFZstd`]. See details at <https://doi.org/10.1074/mcp.O114.037879>
    ///
    /// Requires feature `numpress`
    NumpressSLOF,

    /// **Warning: Lossy** This simply rounds a floating point value and then packs it
    /// into a truncated integer. May be suitable for intensity values, but not consistently better than [`BinaryCompressionType::NumpressSLOF`].
    /// See details at <https://doi.org/10.1074/mcp.O114.037879>
    ///
    /// Requires feature `numpress`
    NumpressPIC,
    /// **Warning: Lossy** Combines [`Self::Zlib`] with [`Self::NumpressLinear`]
    NumpressLinearZlib,
    /// **Warning: Lossy** Combines [`Self::Zlib`] with [`Self::NumpressSLOF`]
    NumpressSLOFZlib,
    /// **Warning: Lossy** Combines [`Self::Zlib`] with [`Self::NumpressPIC`]
    NumpressPICZlib,
    /// **Warning: Lossy** The same concept as [`BinaryCompressionType::NumpressLinear`] without the truncation
    /// All of the same concerns apply though it loses less precision. This does not perform any form of compression. It was uses as part of mzMLb which applies
    /// container-level compression. Only suitable for float64 data.
    LinearPrediction,
    /// **Warning: Lossy** This computes first order difference of a sorted array of values and stores the difference
    /// and the starting point. It is sometimes more accurate than [`BinaryCompressionType::LinearPrediction`]. Like [`BinaryCompressionType::LinearPrediction`],
    /// it does not do any actual compression. It was meant for use with mzMLb, which compresses the entire container. Only suitable for float64 data.
    DeltaPrediction,
    /// The data are completely decoded in native byte order
    Decoded,
    /// A newer general-purpose codec that is not widely used in mzML but very popular in other applications.
    /// Much faster and often better compression compared to `zlib`.
    ///
    /// Unlike [`Self::Zlib`], [`Self::Zstd`]'s default compression level is controlled by [`zstd::DEFAULT_COMPRESSION_LEVEL`],
    /// but the `MZDATA_ZSTD_LEVEL` environment variable can be used to raise or lower it as desired.
    ///
    /// Requires feature `zstd`
    Zstd,
    /// This first applies a byte shuffle transform to the data to compress
    /// the `n`th byte of each value contiguously, and then apply [`Self::Zstd`]. This makes sorted data compress more effectively.
    ///
    /// Requires feature `zstd`
    ShuffleZstd,
    /// **Warning: Lossy** This computes first order difference of a sorted array of values and stores the difference
    /// and the starting point. It then byte shuffles the differences and [`Self::Zstd`] compresses them.
    /// Only suitable for sorted float64 data.
    ///
    /// Requires feature `zstd`
    DeltaShuffleZstd,
    /// Another elaboration of `byte-shuffled zstd` where the values are first dictionary-encoded
    /// so that repeated values are denoted by indices rather than the full size value itself. The dictionary is sorted and then both the dictionary and the
    /// indices are byte shuffle encoded and then [`Self::Zstd`] compressed. This is even more effective than previous codecs for data with repeated, semi-sorted values
    /// like m/z and ion mobility.
    ///
    /// Requires feature `zstd`
    ZstdDict,
    /// **Warning: Lossy** Combines [`Self::Zstd`] with [`Self::NumpressLinear`]
    NumpressLinearZstd,
    /// **Warning: Lossy** Combines [`Self::Zstd`] with [`Self::NumpressSLOF`]
    NumpressSLOFZstd,
    /// **Warning: Lossy** Combines [`Self::Zstd`] with [`Self::NumpressPIC`]
    NumpressPICZstd,
}

impl BinaryCompressionType {
    /// All compression methods in an array
    pub const COMPRESSION_METHODS: &[Self] = &[
        Self::NoCompression,
        Self::Zlib,

        #[cfg(feature = "numpress")]
        Self::NumpressLinear,
        #[cfg(feature = "numpress")]
        Self::NumpressLinearZlib,
        #[cfg(feature = "numpress")]
        Self::NumpressSLOF,
        #[cfg(feature = "numpress")]
        Self::NumpressSLOFZlib,

        #[cfg(feature = "zstd")]
        Self::Zstd,
        #[cfg(feature = "zstd")]
        Self::ShuffleZstd,
        #[cfg(feature = "zstd")]
        Self::DeltaShuffleZstd,
        #[cfg(feature = "zstd")]
        Self::ZstdDict,

        #[cfg(all(feature = "zstd", feature = "numpress"))]
        Self::NumpressLinearZstd,
        #[cfg(all(feature = "zstd", feature = "numpress"))]
        Self::NumpressSLOFZstd,
    ];

    /// Generate a user-understandable message about why a compression conversion operation failed
    pub fn unsupported_msg(&self, context: Option<&str>) -> String {
        match context {
            Some(ctx) => format!("Cannot decode array compressed with {:?} ({})", self, ctx),
            None => format!("Cannot decode array compressed with {:?}", self),
        }
    }

    /// Check if the compressor handles endianess during decoding. Otherwise the data needs to be
    /// byteswapped to little endianess on big endian systems before use
    pub const fn is_endian_aware(&self) -> bool {
        match self {
            BinaryCompressionType::Zlib |
            BinaryCompressionType::NoCompression |
            BinaryCompressionType::DeltaPrediction |
            BinaryCompressionType::LinearPrediction |
            BinaryCompressionType::Zstd => false,
            _ => true,

        }
    }

    /// Get the accession number in the [`ControlledVocabulary::MS`] vocabulary.
    pub const fn accession(&self) -> Option<u32> {
        let acc = match self {
            BinaryCompressionType::NoCompression => 1000576,
            BinaryCompressionType::Zlib => 1000574,
            BinaryCompressionType::NumpressLinear => 1002312,
            BinaryCompressionType::NumpressSLOF => 1002314,
            BinaryCompressionType::NumpressPIC => 1002313,
            BinaryCompressionType::NumpressLinearZlib => 1002746,
            BinaryCompressionType::NumpressSLOFZlib => 1002748,
            BinaryCompressionType::NumpressPICZlib => 1002747,
            BinaryCompressionType::DeltaPrediction => 1003089,
            BinaryCompressionType::LinearPrediction => 1003090,
            BinaryCompressionType::NumpressSLOFZstd => 1003785,
            BinaryCompressionType::NumpressLinearZstd => 1003783,
            BinaryCompressionType::NumpressPICZstd => 1003784,
            BinaryCompressionType::ZstdDict => 1003782,
            BinaryCompressionType::Zstd => 1003780,
            BinaryCompressionType::ShuffleZstd => 1003781,
            BinaryCompressionType::DeltaShuffleZstd => 9999999,
            BinaryCompressionType::Decoded => return None,
        };
        Some(acc)
    }

    /// Convert a [`CURIE`] into a compression method, if one exists.
    pub fn from_accession(accession: CURIE) -> Option<Self> {
        match accession {
            CURIE {
                controlled_vocabulary: ControlledVocabulary::MS,
                accession: 1000576,
            } => Some(Self::NoCompression),
            CURIE {
                controlled_vocabulary: ControlledVocabulary::MS,
                accession: 1000574,
            } => Some(BinaryCompressionType::Zlib),
            CURIE {
                controlled_vocabulary: ControlledVocabulary::MS,
                accession: 1002312,
            } => Some(BinaryCompressionType::NumpressLinear),
            CURIE {
                controlled_vocabulary: ControlledVocabulary::MS,
                accession: 1002314,
            } => Some(BinaryCompressionType::NumpressSLOF),
            CURIE {
                controlled_vocabulary: ControlledVocabulary::MS,
                accession: 1002313,
            } => Some(BinaryCompressionType::NumpressPIC),
            CURIE {
                controlled_vocabulary: ControlledVocabulary::MS,
                accession: 1002746,
            } => Some(BinaryCompressionType::NumpressLinearZlib),
            CURIE {
                controlled_vocabulary: ControlledVocabulary::MS,
                accession: 1002748,
            } => Some(BinaryCompressionType::NumpressSLOFZlib),
            CURIE {
                controlled_vocabulary: ControlledVocabulary::MS,
                accession: 1002747,
            } => Some(BinaryCompressionType::NumpressPICZlib),
            CURIE {
                controlled_vocabulary: ControlledVocabulary::MS,
                accession: 1003089,
            } => Some(BinaryCompressionType::DeltaPrediction),
            CURIE {
                controlled_vocabulary: ControlledVocabulary::MS,
                accession: 1003090,
            } => Some(BinaryCompressionType::LinearPrediction),
            x if x
                == CURIE {
                    controlled_vocabulary: ControlledVocabulary::MS,
                    accession: BinaryCompressionType::NumpressSLOFZstd.accession().unwrap(),
                } =>
            {
                Some(BinaryCompressionType::NumpressSLOFZstd)
            }
            x if x
                == CURIE {
                    controlled_vocabulary: ControlledVocabulary::MS,
                    accession: BinaryCompressionType::NumpressPICZstd.accession().unwrap(),
                } =>
            {
                Some(BinaryCompressionType::NumpressPICZstd)
            }
            x if x
                == CURIE {
                    controlled_vocabulary: ControlledVocabulary::MS,
                    accession: BinaryCompressionType::NumpressLinearZstd
                        .accession()
                        .unwrap(),
                } =>
            {
                Some(BinaryCompressionType::NumpressLinearZstd)
            }
            x if x
                == CURIE {
                    controlled_vocabulary: ControlledVocabulary::MS,
                    accession: BinaryCompressionType::ZstdDict.accession().unwrap(),
                } =>
            {
                Some(BinaryCompressionType::ZstdDict)
            }
            x if x
                == CURIE {
                    controlled_vocabulary: ControlledVocabulary::MS,
                    accession: BinaryCompressionType::Zstd.accession().unwrap(),
                } =>
            {
                Some(BinaryCompressionType::Zstd)
            }
            x if x
                == CURIE {
                    controlled_vocabulary: ControlledVocabulary::MS,
                    accession: BinaryCompressionType::ShuffleZstd.accession().unwrap(),
                } =>
            {
                Some(BinaryCompressionType::ShuffleZstd)
            }
            x if x
                == CURIE {
                    controlled_vocabulary: ControlledVocabulary::MS,
                    accession: BinaryCompressionType::DeltaShuffleZstd.accession().unwrap(),
                } =>
            {
                Some(BinaryCompressionType::DeltaShuffleZstd)
            }
            _ => None,
        }
    }

    /// Convert the compression type to a [`ParamCow`].
    ///
    /// Most compression methods have a controlled vocabulary
    /// term.
    pub const fn as_param(&self) -> Option<ParamCow<'static>> {
        let (name, accession) = match self {
            BinaryCompressionType::Decoded => return None,
            BinaryCompressionType::NoCompression => ("no compression", self.accession()),
            BinaryCompressionType::Zlib => ("zlib compression", self.accession()),
            BinaryCompressionType::NumpressLinear => (
                "MS-Numpress linear prediction compression",
                self.accession(),
            ),
            BinaryCompressionType::NumpressSLOF => (
                "MS-Numpress short logged float compression",
                self.accession(),
            ),
            BinaryCompressionType::NumpressPIC => {
                ("MS-Numpress positive integer compression", self.accession())
            }
            BinaryCompressionType::NumpressLinearZlib => (
                "MS-Numpress linear prediction compression followed by zlib compression",
                self.accession(),
            ),
            BinaryCompressionType::NumpressSLOFZlib => (
                "MS-Numpress short logged float compression followed by zlib compression",
                self.accession(),
            ),
            BinaryCompressionType::NumpressPICZlib => (
                "MS-Numpress positive integer compression followed by zlib compression",
                self.accession(),
            ),
            BinaryCompressionType::DeltaPrediction => (
                "truncation, delta prediction and zlib compression",
                self.accession(),
            ),
            BinaryCompressionType::LinearPrediction => (
                "truncation, linear prediction and zlib compression",
                self.accession(),
            ),
            BinaryCompressionType::NumpressSLOFZstd => {
                return Some(ParamCow::const_new(
                    "MS-Numpress short logged float compression followed by zstd compression",
                    ValueRef::Empty,
                    self.accession(),
                    Some(ControlledVocabulary::MS),
                    Unit::Unknown,
                ))
            }
            BinaryCompressionType::NumpressLinearZstd => {
                return Some(ParamCow::const_new(
                    "MS-Numpress linear prediction compression followed by zstd compression",
                    ValueRef::Empty,
                    self.accession(),
                    Some(ControlledVocabulary::MS),
                    Unit::Unknown,
                ))
            }
            BinaryCompressionType::NumpressPICZstd => {
                return Some(ParamCow::const_new(
                    "MS-Numpress positive integer compression followed by zstd compression",
                    ValueRef::Empty,
                    self.accession(),
                    Some(ControlledVocabulary::MS),
                    Unit::Unknown,
                ))
            }
            BinaryCompressionType::ZstdDict => {
                return Some(ParamCow::const_new(
                    "dict-zstd compression",
                    ValueRef::Empty,
                    self.accession(),
                    Some(ControlledVocabulary::MS),
                    Unit::Unknown,
                ))
            }
            BinaryCompressionType::Zstd => {
                return Some(ParamCow::const_new(
                    "zstd compression",
                    ValueRef::Empty,
                    self.accession(),
                    Some(ControlledVocabulary::MS),
                    Unit::Unknown,
                ))
            }
            BinaryCompressionType::ShuffleZstd => {
                return Some(ParamCow::const_new(
                    "byte-shuffle-zstd compression",
                    ValueRef::Empty,
                    self.accession(),
                    Some(ControlledVocabulary::MS),
                    Unit::Unknown,
                ))
            }
            BinaryCompressionType::DeltaShuffleZstd => {
                return Some(ParamCow::const_new(
                    "delta-byte-shuffle-zstd compression",
                    ValueRef::Empty,
                    self.accession(),
                    Some(ControlledVocabulary::MS),
                    Unit::Unknown,
                ))
            }
        };
        Some(
            ControlledVocabulary::MS
                .const_param_ident(name, unsafe { accession.unwrap_unchecked() }),
        )
    }
}

impl Display for BinaryCompressionType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

/// A high level set of failure modes that an operation to retrieve a typed memory buffer
/// from a [`BinaryArrayMap`](crate::BinaryArrayMap) might encounter. May also be used to represented conversion
/// during reading or writing.
#[derive(Debug, Clone, Error, PartialEq)]
pub enum ArrayRetrievalError {
    #[error("Array type {0:?} not found")]
    NotFound(ArrayType),
    #[error("An error occurred while decompressing: {0}")]
    DecompressionError(String),
    #[error("The requested data type does not match the number of bytes available in the buffer, or the conversion is unsupported")]
    DataTypeSizeMismatch,
    #[error("The content was not properly Base64 encoded")]
    MalformedBase64Encoding,
}

impl From<bytemuck::PodCastError> for ArrayRetrievalError {
    fn from(value: bytemuck::PodCastError) -> Self {
        match value {
            bytemuck::PodCastError::TargetAlignmentGreaterAndInputNotAligned => {
                Self::DataTypeSizeMismatch
            }
            bytemuck::PodCastError::OutputSliceWouldHaveSlop => Self::DataTypeSizeMismatch,
            bytemuck::PodCastError::SizeMismatch => Self::DataTypeSizeMismatch,
            bytemuck::PodCastError::AlignmentMismatch => Self::DataTypeSizeMismatch,
        }
    }
}

impl From<ArrayRetrievalError> for io::Error {
    fn from(value: ArrayRetrievalError) -> Self {
        match value {
            ArrayRetrievalError::NotFound(_) => io::Error::new(io::ErrorKind::NotFound, value),
            ArrayRetrievalError::DecompressionError(e) => {
                io::Error::new(io::ErrorKind::InvalidData, e)
            }
            ArrayRetrievalError::DataTypeSizeMismatch => {
                io::Error::new(io::ErrorKind::InvalidData, value)
            },
            ArrayRetrievalError::MalformedBase64Encoding => {
                io::Error::new(io::ErrorKind::InvalidData, value)
            }
        }
    }
}

#[cfg(feature = "numpress")]
impl From<numpress::Error> for ArrayRetrievalError {
    fn from(value: numpress::Error) -> Self {
        ArrayRetrievalError::DecompressionError(value.to_string())
    }
}

/// The decoder for [`BinaryCompressionType::LinearPrediction`].
pub fn linear_prediction_decoding<F: Num + Copy + Mul + AddAssign>(values: &mut [F]) -> &mut [F] {
    if values.len() < 2 {
        return values;
    }

    let two = F::one() + F::one();

    let prev2 = values[1];
    let prev1 = values[2];
    let offset = values[1];

    values
        .iter_mut()
        .skip(2)
        .fold((prev1, prev2), |(prev1, prev2), current| {
            let tmp = *current + two * prev1 - prev2 - offset;
            let prev1 = *current;
            let prev2 = prev1;
            *current = tmp;
            (prev1, prev2)
        });

    for i in 0..values.len() {
        if i < 2 {
            continue;
        }
        let v = values[i] + two * values[i - 1] - values[i - 2] - values[1];
        values[i] = v;
    }
    values
}

/// The encoder for [`BinaryCompressionType::LinearPrediction`].
pub fn linear_prediction_encoding<F: Num + Copy + Mul<F> + AddAssign>(
    values: &mut [F],
) -> &mut [F] {
    let n = values.len();
    if n < 3 {
        return values;
    }
    let offset = values[1];
    let prev2 = values[0];
    let prev1 = values[1];
    let two = F::one() + F::one();

    values
        .iter_mut()
        .fold((prev1, prev2), |(prev1, prev2), val| {
            *val += offset - two * prev1 + prev2;
            let tmp = prev1;
            let prev1 = *val + two * prev1 - prev2 - offset;
            let prev2 = tmp;
            (prev1, prev2)
        });
    values
}

/// The decoder for [`BinaryCompressionType::DeltaPrediction`].
pub fn delta_decoding<F: Num + Copy + Mul + AddAssign>(values: &mut [F]) -> &mut [F] {
    if values.len() < 2 {
        return values;
    }

    let offset = values[0];
    let prev = values[1];

    values.iter_mut().skip(2).fold(prev, |prev, current| {
        *current += prev - offset;
        *current
    });
    values
}

/// The encoder for [`BinaryCompressionType::DeltaPrediction`].
pub fn delta_encoding<F: Num + Copy + Mul + AddAssign>(values: &mut [F]) -> &mut [F] {
    let n = values.len();
    if n < 2 {
        return values;
    }
    let prev = values[0];
    let offset = values[0];

    let it = values.iter_mut();
    it.skip(1).fold(prev, |prev, current| {
        let tmp = *current;
        *current += offset - prev;
        tmp
    });
    values
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_dtype_size() {
        assert_eq!(BinaryDataArrayType::ASCII.size_of(), 1);
        assert_eq!(BinaryDataArrayType::Float32.size_of(), 4);
        assert_eq!(BinaryDataArrayType::Int32.size_of(), 4);
        assert_eq!(BinaryDataArrayType::Float64.size_of(), 8);
        assert_eq!(BinaryDataArrayType::Int64.size_of(), 8);
    }

    #[test]
    fn test_array_type_param() {
        let array_types = [
            ArrayType::MZArray,
            ArrayType::IntensityArray,
            ArrayType::ChargeArray,
            ArrayType::SignalToNoiseArray,
            ArrayType::TimeArray,
            ArrayType::WavelengthArray,
            ArrayType::IonMobilityArray,
            ArrayType::MeanIonMobilityArray,
            ArrayType::RawIonMobilityArray,
            ArrayType::DeconvolutedIonMobilityArray,
        ];

        for at in array_types {
            assert_eq!(at.as_param_const().name, at.as_param(None).name)
        }
    }

    #[test]
    fn test_binary_encoding_conv() {
        let encodings = [
            BinaryCompressionType::Decoded,
            BinaryCompressionType::NoCompression,
            BinaryCompressionType::NumpressLinear,
            BinaryCompressionType::NumpressLinearZlib,
            BinaryCompressionType::NumpressPIC,
            BinaryCompressionType::NumpressPICZlib,
            BinaryCompressionType::NumpressSLOF,
            BinaryCompressionType::NumpressSLOFZlib,
            BinaryCompressionType::Zlib,
        ];

        for enc in encodings {
            let reps = match enc {
                BinaryCompressionType::NoCompression => ("no compression", 1000576),
                BinaryCompressionType::Zlib => ("zlib compression", 1000574),
                BinaryCompressionType::NumpressLinear => {
                    ("MS-Numpress linear prediction compression", 1002312)
                }
                BinaryCompressionType::NumpressSLOF => {
                    ("MS-Numpress short logged float compression", 1002314)
                }
                BinaryCompressionType::NumpressPIC => {
                    ("MS-Numpress positive integer compression", 1002313)
                }
                BinaryCompressionType::NumpressLinearZlib => (
                    "MS-Numpress linear prediction compression followed by zlib compression",
                    1002746,
                ),
                BinaryCompressionType::NumpressPICZlib => (
                    "MS-Numpress positive integer compression followed by zlib compression",
                    1002747,
                ),
                BinaryCompressionType::NumpressSLOFZlib => (
                    "MS-Numpress short logged float compression followed by zlib compression",
                    1002748,
                ),
                _ => ("", 0),
            };
            if let Some(p) = enc.as_param() {
                assert_eq!(p.name, reps.0);
                assert_eq!(p.accession.unwrap(), reps.1);
            }
        }
    }

    #[test]
    fn test_transpose() {
        let data: Vec<_> = (0..128i32).map(|i| i.pow(2u32) as f64).collect();
        let flip = transpose_f64(&data);
        let rev = reverse_transpose_f64(&flip);
        let rev_cast: &[f64] = bytemuck::cast_slice(&rev);
        assert_eq!(data, rev_cast);
    }

    #[test]
    fn test_dict() {
        let data: Vec<_> = (0..127i32).map(|i| i.pow(2u32) as f64).collect();
        let encoded = dictionary_encoding(&data).unwrap();
        let decoded: Vec<f64> = dictionary_decoding(&encoded).unwrap();

        assert_eq!(data, decoded);
    }

    #[test]
    fn test_byteswap() {
        let x = 42u32;
        let mut bytes_of = x.to_le_bytes();
        BinaryDataArrayType::Int32.swap_bytes(&mut bytes_of).unwrap();
        let y1 = u32::from_le_bytes(bytes_of);
        let y2 = x.swap_bytes();
        assert_eq!(y1, y2);

        bytes_of = x.to_be_bytes();
        BinaryDataArrayType::Int32.swap_bytes(&mut bytes_of).unwrap();
        let y1 = u32::from_le_bytes(bytes_of);
        assert_eq!(x, y1);
    }
}
