use std::io;
use std::path::PathBuf;

use futures::stream;
use mzpeaks::{CentroidPeak, DeconvolutedPeak};
use tokio;

use super::ThermoRawReaderType as SyncThermoRawReaderType;
use crate::{
    io::{
        DetailLevel,
        traits::{AsyncMZFileReader, AsyncRandomAccessSpectrumIterator, SpectrumStream},
    },
    prelude::*,
    spectrum::MultiLayerSpectrum,
};

pub struct ThermoRawReaderType<
    C: CentroidLike + From<CentroidPeak> + Send = CentroidPeak,
    D: DeconvolutedCentroidLike + Send = DeconvolutedPeak,
> {
    inner: Option<SyncThermoRawReaderType<C, D>>,
}

#[cfg(feature = "async")]
impl<
    C: CentroidLike + From<CentroidPeak> + Send + 'static,
    D: DeconvolutedCentroidLike + Send + 'static,
> AsyncMZFileReader<C, D, MultiLayerSpectrum<C, D>> for ThermoRawReaderType<C, D>
{
    async fn construct_index_from_stream(&mut self) -> u64 {
        self.len() as u64
    }

    /// The underlying Thermo library requires an explicit file system path to open
    /// the file, and as such this method always fails.
    ///
    /// [`open_path`](Self::open_path) works as normal.
    #[allow(unused)]
    async fn open_file(source: tokio::fs::File) -> io::Result<Self> {
        Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "Cannot read a Thermo RAW file from an open file handle, only directly from a path",
        ))
    }

    async fn open_path<P>(path: P) -> io::Result<Self>
    where
        P: Into<std::path::PathBuf>,
    {
        Self::new(path.into()).await
    }
}

impl<C: CentroidLike + From<CentroidPeak> + Send, D: DeconvolutedCentroidLike + Send>
    MSDataFileMetadata for ThermoRawReaderType<C, D>
{
    crate::delegate_impl_metadata_trait!(expr, this => { this.inner.as_ref().unwrap() }, &mut => { this.inner.as_mut().unwrap() });
}

impl<
    C: CentroidLike + From<CentroidPeak> + Send + 'static,
    D: DeconvolutedCentroidLike + Send + 'static,
> AsyncSpectrumSource<C, D, MultiLayerSpectrum<C, D>> for ThermoRawReaderType<C, D>
{
    fn reset(&mut self) -> impl std::future::Future<Output = ()> {
        self.inner.as_mut().unwrap().reset();
        futures::future::ready(())
    }

    fn detail_level(&self) -> &DetailLevel {
        &self.inner.as_ref().unwrap().detail_level
    }

    fn set_detail_level(&mut self, detail_level: DetailLevel) {
        self.inner.as_mut().unwrap().set_detail_level(detail_level)
    }

    fn get_spectrum_by_id(
        &mut self,
        id: &str,
    ) -> impl std::future::Future<Output = Option<MultiLayerSpectrum<C, D>>> {
        self.get_spectrum_by_id(id)
    }

    fn get_spectrum_by_index(
        &mut self,
        index: usize,
    ) -> impl std::future::Future<Output = Option<MultiLayerSpectrum<C, D>>> {
        self.get_spectrum_by_index(index)
    }

    fn get_index(&self) -> &crate::io::OffsetIndex {
        self.get_index()
    }

    fn set_index(&mut self, index: crate::io::OffsetIndex) {
        self.inner.as_mut().unwrap().set_index(index);
    }

    fn read_next(&mut self) -> impl std::future::Future<Output = Option<MultiLayerSpectrum<C, D>>> {
        self.read_next()
    }

    async fn get_spectrum_by_time(&mut self, time: f64) -> Option<MultiLayerSpectrum<C, D>> {
        self.get_spectrum_by_time(time).await
    }
}

impl<
    C: CentroidLike + From<CentroidPeak> + Send + 'static,
    D: DeconvolutedCentroidLike + Send + 'static,
> ThermoRawReaderType<C, D>
{
    /// Create a new [`ThermoRawReaderType`] from a path.
    /// This may trigger an expensive I/O operation to checksum the file
    pub async fn new<P: Into<PathBuf> + 'static + Send>(path: P) -> io::Result<Self> {
        Self::new_with_detail_level_and_centroiding(path, DetailLevel::Full, false).await
    }

    /// Create a new [`ThermoRawReaderType`] from a path.
    /// This may trigger an expensive I/O operation to checksum the file
    pub async fn new_with_detail_level_and_centroiding<P: Into<PathBuf> + Send + 'static>(
        path: P,
        detail_level: DetailLevel,
        centroiding: bool,
    ) -> io::Result<Self> {
        tokio::task::spawn_blocking(move || {
            let inner = SyncThermoRawReaderType::new_with_detail_level_and_centroiding(
                path,
                detail_level,
                centroiding,
            )?;
            let this = Self { inner: Some(inner) };
            Ok(this)
        })
        .await
        .unwrap()
    }

    pub fn len(&self) -> usize {
        self.inner.as_ref().unwrap().len()
    }

    pub fn is_empty(&self) -> bool {
        self.inner.as_ref().unwrap().is_empty()
    }

    pub fn get_centroiding(&self) -> bool {
        self.inner.as_ref().unwrap().get_centroiding()
    }

    pub fn set_centroiding(&mut self, value: bool) {
        self.inner.as_mut().unwrap().set_centroiding(value)
    }

    /// Get whether or not to load extended spectrum signal information for the spectrum.
    ///
    /// The loaded data isn't incorporated into a peak list, instead access them under
    /// the binary data arrays.
    pub fn get_load_extended_spectrum_data(&self) -> bool {
        self.inner
            .as_ref()
            .unwrap()
            .get_load_extended_spectrum_data()
    }

    /// Set whether or not to load extended spectrum signal information for the spectrum.
    ///
    /// The loaded data isn't incorporated into a peak list, instead access them under
    /// the binary data arrays.
    pub fn set_load_extended_spectrum_data(&mut self, load_extended_spectrum_data: bool) {
        self.inner
            .as_mut()
            .unwrap()
            .set_load_extended_spectrum_data(load_extended_spectrum_data)
    }

    pub fn get_index(&self) -> &crate::io::OffsetIndex {
        self.inner.as_ref().unwrap().get_index()
    }

    pub fn as_stream(&mut self) -> impl SpectrumStream<C, D, MultiLayerSpectrum<C, D>> + '_ {
        Box::pin(stream::unfold(self, |reader| async {
            let spec = reader.read_next();
            spec.await.map(|val| (val, reader))
        }))
    }

    pub async fn read_next(&mut self) -> Option<MultiLayerSpectrum<C, D>> {
        let mut inner = self.inner.take().unwrap();
        let (inner, spec) = tokio::task::spawn_blocking(move || {
            let spec = inner.read_next_spectrum();
            (inner, spec)
        })
        .await
        .unwrap();
        self.inner = Some(inner);
        spec
    }

    pub async fn get_spectrum_by_id(&mut self, id: &str) -> Option<MultiLayerSpectrum<C, D>> {
        let mut inner = self.inner.take().unwrap();
        let id = id.to_string();
        let (inner, spec) = tokio::task::spawn_blocking(move || {
            let spec = inner.get_spectrum_by_id(&id);
            (inner, spec)
        })
        .await
        .unwrap();
        self.inner = Some(inner);
        spec
    }

    pub async fn get_spectrum_by_index(
        &mut self,
        index: usize,
    ) -> Option<MultiLayerSpectrum<C, D>> {
        let mut inner = self.inner.take().unwrap();
        let (inner, spec) = tokio::task::spawn_blocking(move || {
            let spec = inner.get_spectrum_by_index(index);
            (inner, spec)
        })
        .await
        .unwrap();
        self.inner = Some(inner);
        spec
    }

    pub async fn get_spectrum_by_time(&mut self, time: f64) -> Option<MultiLayerSpectrum<C, D>> {
        let mut inner = self.inner.take().unwrap();
        let (inner, spec) = tokio::task::spawn_blocking(move || {
            let spec = inner.get_spectrum_by_time(time);
            (inner, spec)
        })
        .await
        .unwrap();
        self.inner = Some(inner);
        spec
    }
}

pub type ThermoRawReader = ThermoRawReaderType<CentroidPeak, DeconvolutedPeak>;

impl<
    C: CentroidLike + From<CentroidPeak> + Send + Sync + 'static,
    D: DeconvolutedCentroidLike + Send + Sync + 'static,
> AsyncRandomAccessSpectrumIterator<C, D, MultiLayerSpectrum<C, D>> for ThermoRawReaderType<C, D>
{
    async fn start_from_id(&mut self, id: &str) -> Result<&mut Self, SpectrumAccessError> {
        let mut inner = self.inner.take().unwrap();
        let id = id.to_string();
        let (inner, spec) = tokio::task::spawn_blocking(move || {
            let spec = inner.start_from_id(&id);
            let res = spec.err();
            (inner, res)
        })
        .await
        .unwrap();
        if let Some(e) = spec {
            return Err(e);
        }
        self.inner = Some(inner);
        Ok(self)
    }

    async fn start_from_index(&mut self, index: usize) -> Result<&mut Self, SpectrumAccessError> {
        let mut inner = self.inner.take().unwrap();
        let (inner, spec) = tokio::task::spawn_blocking(move || {
            let spec = inner.start_from_index(index);
            let res = spec.err();
            (inner, res)
        })
        .await
        .unwrap();
        if let Some(e) = spec {
            return Err(e);
        }
        self.inner = Some(inner);
        Ok(self)
    }

    async fn start_from_time(&mut self, time: f64) -> Result<&mut Self, SpectrumAccessError> {
        let mut inner = self.inner.take().unwrap();
        let (inner, spec) = tokio::task::spawn_blocking(move || {
            let spec = inner.start_from_time(time);
            let res = spec.err();
            (inner, res)
        })
        .await
        .unwrap();
        if let Some(e) = spec {
            return Err(e);
        }
        self.inner = Some(inner);
        Ok(self)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[tokio::test(flavor = "multi_thread", worker_threads = 4)]
    async fn test_read() -> io::Result<()> {
        let mut reader: ThermoRawReaderType<CentroidPeak, DeconvolutedPeak> =
            ThermoRawReaderType::open_path("./test/data/small.RAW").await?;

        let n = reader.len();
        let mut ms1_counter = 0;
        let mut msn_counter = 0;
        while let Some(spec) = reader.read_next().await {
            if spec.ms_level() > 1 {
                msn_counter += 1;
            } else {
                ms1_counter += 1;
            }
        }

        assert_eq!(n, ms1_counter + msn_counter);
        Ok(())
    }
}
