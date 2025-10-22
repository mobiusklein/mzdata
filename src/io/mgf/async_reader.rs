use core::str;
use std::io::SeekFrom;
use std::marker::PhantomData;

use std::collections::HashMap;

use futures::stream::{self};
use tokio::io::{self, AsyncBufReadExt, AsyncSeekExt};

use log::{error, warn};

use mzpeaks::{CentroidLike, CentroidPeak, DeconvolutedCentroidLike, DeconvolutedPeak};

use super::{
    super::{offset_index::OffsetIndex, utils::DetailLevel},
    reader::{MGFLineParsing, SpectrumBuilder},
    MGFError, MGFParserState, MGFIndexing, DefaultTitleIndexing
};

#[cfg(feature = "async")]
use crate::io::traits::AsyncMZFileReader;

use crate::{
    io::traits::{AsyncRandomAccessSpectrumIterator, SpectrumStream},
    meta::{
        DataProcessing, FileDescription, InstrumentConfiguration, MSDataFileMetadata,
        MassSpectrometryRun, Sample, Software,
    },
    prelude::{SpectrumAccessError, SpectrumLike},
};

use crate::params::{ControlledVocabulary, Param, ParamDescribed};

use crate::spectrum::spectrum_types::MultiLayerSpectrum;

use crate::io::traits::AsyncSpectrumSource;

/// An MGF (Mascot Generic Format) file parser that supports iteration and random access.
/// The parser produces [`Spectrum`](crate::spectrum::Spectrum) instances. These may be
/// converted directly into [`CentroidSpectrum`](crate::spectrum::CentroidSpectrum)
/// instances which better represent the nature of this preprocessed data type.
pub struct MGFReaderType<
    R: io::AsyncRead,
    C: CentroidLike + From<CentroidPeak> = CentroidPeak,
    D: DeconvolutedCentroidLike + From<DeconvolutedPeak> = DeconvolutedPeak,
> {
    pub handle: io::BufReader<R>,
    pub state: MGFParserState,
    pub offset: usize,
    pub error: Option<MGFError>,
    index: OffsetIndex,
    file_description: FileDescription,
    instrument_configurations: HashMap<u32, InstrumentConfiguration>,
    softwares: Vec<Software>,
    samples: Vec<Sample>,
    data_processings: Vec<DataProcessing>,
    run: MassSpectrometryRun,
    pub detail_level: DetailLevel,
    centroid_type: PhantomData<C>,
    deconvoluted_type: PhantomData<D>,
    read_counter: usize,
    indexer: Box<dyn MGFIndexing>,
}

#[cfg(feature = "async")]
impl<
        C: CentroidLike + From<CentroidPeak> + Send + Sync,
        D: DeconvolutedCentroidLike + From<DeconvolutedPeak> + Send + Sync,
    > AsyncMZFileReader<C, D, MultiLayerSpectrum<C, D>> for MGFReaderType<tokio::fs::File, C, D>
{
    async fn construct_index_from_stream(&mut self) -> u64 {
        self.build_index().await
    }

    async fn open_file(
        source: tokio::fs::File,
    ) -> std::io::Result<Self> {
        Ok(Self::new(source).await)
    }
}

impl<R: io::AsyncRead, C: CentroidLike + From<CentroidPeak>, D: DeconvolutedCentroidLike + From<DeconvolutedPeak>> MGFLineParsing<C, D>
    for MGFReaderType<R, C, D>
{
    fn state(&self) -> &MGFParserState {
        &self.state
    }

    fn state_mut(&mut self) -> &mut MGFParserState {
        &mut self.state
    }

    fn error_mut(&mut self) -> &mut Option<MGFError> {
        &mut self.error
    }

    fn indexer(&self) -> &Box<dyn MGFIndexing> {
        &self.indexer
    }

    fn set_indexer(&mut self, indexer: Box<dyn MGFIndexing>) {
        self.indexer = indexer;
        self.index = OffsetIndex::new("spectrum".into());
    }
}

impl<R: io::AsyncRead + Unpin, C: CentroidLike + From<CentroidPeak>, D: DeconvolutedCentroidLike + From<DeconvolutedPeak>>
    MGFReaderType<R, C, D>
{
    async fn read_line(&mut self, buffer: &mut String) -> io::Result<usize> {
        self.handle.read_line(buffer).await
    }

    /// Read the next spectrum from the file, if there is one.
    pub async fn read_next(&mut self) -> Option<MultiLayerSpectrum<C, D>> {
        let mut builder = SpectrumBuilder::<C, D>::default();
        self._parse_into(&mut builder)
            .await
            .inspect_err(|e| error!("An error occurred while reading MGF spectrum: {e}"))
            .ok()
            .and_then(|(_, started_spectrum)| {
                let mut spec: Option<MultiLayerSpectrum<C, D>> =
                    (started_spectrum && !builder.is_empty()).then(|| builder.into());
                if let Some(spec) = spec.as_mut() {
                    spec.description_mut().index = self.read_counter;
                    self.read_counter += 1;
                }
                spec
            })
    }

    /// Read the next spectrum's contents directly into the passed [`SpectrumBuilder`].
    async fn _parse_into(
        &mut self,
        builder: &mut SpectrumBuilder<C, D>,
    ) -> Result<(usize, bool), MGFError> {
        let mut buffer = String::new();
        let mut work = true;
        let mut offset: usize = 0;
        let mut had_begin_ions = false;

        while work {
            buffer.clear();
            let b = match self.read_line(&mut buffer).await {
                Ok(b) => b,
                Err(err) => {
                    self.state = MGFParserState::Error;
                    return Err(MGFError::IOError(err));
                }
            };

            // Count how many bytes we've read from the source
            offset += b;
            if b == 0 {
                self.state = MGFParserState::Done;
                break;
            }

            let line = buffer.trim();

            // Skip empty lines
            if line.is_empty() {
                continue;
            }

            work = match self.state {
                MGFParserState::Start | MGFParserState::FileHeader => self.handle_start(line),
                MGFParserState::ScanHeaders => {
                    had_begin_ions = true;
                    self.handle_scan_header(line, builder)
                }
                MGFParserState::Peaks => self.handle_peak(line, builder),
                MGFParserState::Between => self.handle_between(line),
                MGFParserState::Done => false,
                MGFParserState::Error => {
                    return Err(self.error.take().unwrap_or(MGFError::NoError));
                }
            };

            if self.state == MGFParserState::Error {
                return Err(self.error.take().unwrap_or(MGFError::NoError));
            }
        }
        Ok((offset, had_begin_ions))
    }

    pub async fn read_into(
        &mut self,
        spectrum: &mut MultiLayerSpectrum<C, D>,
    ) -> Result<usize, MGFError> {
        let mut accumulator = SpectrumBuilder::default();
        match self._parse_into(&mut accumulator).await {
            Ok((sz, started_spectrum)) => {
                if !started_spectrum {
                    Err(MGFError::IOError(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "EOF found before spectrum started",
                    )))
                } else {
                    accumulator.into_spectrum(spectrum);
                    spectrum.description_mut().index = self.read_counter;
                    self.read_counter += 1;
                    Ok(sz)
                }
            }
            Err(err) => Err(err),
        }
    }

    fn default_file_description() -> FileDescription {
        let mut fd = FileDescription::default();
        let mut term = Param::new();
        term.name = "MSn spectrum".to_owned();
        term.accession = Some(1000580);
        term.controlled_vocabulary = Some(ControlledVocabulary::MS);
        fd.add_param(term);
        fd
    }

    /// Create a new, unindexed MGF parser
    pub async fn new(file: R) -> MGFReaderType<R, C, D> {
        let handle = io::BufReader::with_capacity(500, file);
        MGFReaderType {
            handle,
            state: MGFParserState::Start,
            offset: 0,
            error: None,
            index: OffsetIndex::new("spectrum".to_owned()),
            centroid_type: PhantomData,
            deconvoluted_type: PhantomData,
            instrument_configurations: HashMap::new(),
            data_processings: Vec::new(),
            softwares: Vec::new(),
            samples: Vec::new(),
            file_description: Self::default_file_description(),
            detail_level: DetailLevel::Full,
            run: MassSpectrometryRun::default(),
            read_counter: 0,
            indexer: Box::new(DefaultTitleIndexing::default()),
        }
    }
}

impl<
        R: io::AsyncRead + io::AsyncSeek + Unpin,
        C: CentroidLike + From<CentroidPeak>,
        D: DeconvolutedCentroidLike + From<DeconvolutedPeak>,
    > MGFReaderType<R, C, D>
{
    pub fn as_stream(&mut self) -> impl SpectrumStream<C, D, MultiLayerSpectrum<C, D>> + '_ {
        Box::pin(stream::unfold(self, |reader| async {
            let spec = reader.read_next();
            spec.await.map(|val| (val, reader))
        }))
    }

    /// Construct a new MGFReaderType and build an offset index
    /// using [`Self::build_index`]
    pub async fn new_indexed(file: R) -> MGFReaderType<R, C, D> {
        let mut reader = Self::new(file).await;
        reader.build_index().await;
        reader
    }

    pub async fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        self.handle.seek(pos).await
    }

    /// Builds an offset index to each `BEGIN IONS` line
    /// by doing a fast pre-scan of the text file.
    pub async fn build_index(&mut self) -> u64 {
        let mut offset: u64 = 0;
        let mut last_start: u64 = 0;

        let mut found_start = false;

        let start = self
            .handle
            .stream_position()
            .await
            .expect("Failed to save restore location");
        self.seek(SeekFrom::Start(0))
            .await
            .expect("Failed to reset stream to beginning");

        let mut buffer: Vec<u8> = Vec::new();

        self.index.clear();

        loop {
            buffer.clear();
            let b = match self.handle.read_until(b'\n', &mut buffer).await {
                Ok(b) => b,
                Err(err) => {
                    panic!("Error while reading file: {}", err);
                }
            };
            if b == 0 {
                break;
            }
            if buffer.starts_with(b"BEGIN IONS") {
                found_start = true;
                last_start = offset;
            } else if found_start {
                let string_buffer = String::from_utf8_lossy(&buffer);
                let indexer = self.indexer();
                if let Some((key, value)) = string_buffer.split_once('=') {
                    if indexer.is_index_key(key) {
                        self.index.insert(indexer.handle_key(key, value), last_start);
                        found_start = false;
                        last_start = 0;
                    }
                }
            }
            offset += b as u64;
        }
        self.seek(SeekFrom::Start(start))
            .await
            .expect("Failed to restore location");
        self.index.init = true;
        if self.index.is_empty() {
            warn!("An index was built but no entries were found")
        }
        offset
    }
}

/// The MGF format does not contain any consistent metadata, but additional
/// information can be included after creation.
impl<R: io::AsyncRead + Unpin, C: CentroidLike + From<CentroidPeak>, D: DeconvolutedCentroidLike + From<DeconvolutedPeak>>
    MSDataFileMetadata for MGFReaderType<R, C, D>
{
    crate::impl_metadata_trait!();

    fn spectrum_count_hint(&self) -> Option<u64> {
        if self.index.init {
            Some(self.index.len() as u64)
        } else {
            None
        }
    }

    fn run_description(&self) -> Option<&MassSpectrometryRun> {
        Some(&self.run)
    }

    fn run_description_mut(&mut self) -> Option<&mut MassSpectrometryRun> {
        Some(&mut self.run)
    }
}

pub type MGFReader<R> = MGFReaderType<R, CentroidPeak, DeconvolutedPeak>;

impl<
        R: io::AsyncRead + io::AsyncSeek + io::AsyncSeekExt + Unpin,
        C: CentroidLike + From<CentroidPeak> + Send + Sync,
        D: DeconvolutedCentroidLike + From<DeconvolutedPeak> + Send + Sync,
    > MGFReaderType<R, C, D>
{
    /// Helper method to support seeking to an ID
    fn _offset_of_id(&self, id: &str) -> Option<u64> {
        self.get_index().get(id)
    }

    /// Helper method to support seeking to an index
    fn _offset_of_index(&self, index: usize) -> Option<u64> {
        self.get_index()
            .get_index(index)
            .map(|(_id, offset)| offset)
    }

    /// Helper method to support seeking to a specific time.
    /// Considerably more complex than seeking by ID or index.
    async fn _offset_of_time(&mut self, time: f64) -> Option<u64> {
        match self.get_spectrum_by_time(time).await {
            Some(scan) => self._offset_of_index(scan.index()),
            None => None,
        }
    }

    /// Read the length of the spectrum offset index
    pub fn len(&self) -> usize {
        self.index.len()
    }

    pub fn is_empty(&self) -> bool {
        self.index.is_empty()
    }

    /// Retrieve a spectrum by its scan start time
    /// Considerably more complex than seeking by ID or index, this involves
    /// a binary search over the spectrum index and assumes that spectra are stored
    /// in chronological order.
    pub async fn get_spectrum_by_time(&mut self, time: f64) -> Option<MultiLayerSpectrum<C, D>> {
        let n = self.len();
        let mut lo: usize = 0;
        let mut hi: usize = n;

        let mut best_error: f64 = f64::INFINITY;
        let mut best_match: Option<MultiLayerSpectrum<_, _>> = None;

        if lo == hi {
            return None;
        }
        while hi != lo {
            let mid = (hi + lo) / 2;
            let scan = self.get_spectrum_by_index(mid).await?;
            let scan_time = scan.start_time();
            let err = (scan_time - time).abs();

            if err < best_error {
                best_error = err;
                best_match = Some(scan);
            }
            if hi.saturating_sub(1) == lo {
                return best_match;
            } else if scan_time > time {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        best_match
    }

    /// Retrieve a spectrum by it's native ID
    pub async fn get_spectrum_by_id(&mut self, id: &str) -> Option<MultiLayerSpectrum<C, D>> {
        let offset = self.index.get(id)?;
        let index = self.index.index_of(id)?;
        let start = self
            .handle
            .stream_position()
            .await
            .expect("Failed to save checkpoint");
        self.seek(SeekFrom::Start(offset))
            .await
            .expect("Failed to move seek to offset");
        let result = self.read_next().await;
        self.seek(SeekFrom::Start(start))
            .await
            .expect("Failed to restore offset");
        result.map(|mut scan| {
            scan.description.index = index;
            scan
        })
    }

    /// Retrieve a spectrum by it's integer index
    pub async fn get_spectrum_by_index(
        &mut self,
        index: usize,
    ) -> Option<MultiLayerSpectrum<C, D>> {
        let (_id, byte_offset) = self.index.get_index(index)?;
        let start = self
            .handle
            .stream_position()
            .await
            .expect("Failed to save checkpoint");
        self.seek(SeekFrom::Start(byte_offset)).await.ok()?;
        let result = self.read_next().await;
        self.seek(SeekFrom::Start(start))
            .await
            .expect("Failed to restore offset");
        result.map(|mut scan| {
            scan.description.index = index;
            scan
        })
    }

    /// Return the data stream to the beginning
    pub async fn reset(&mut self) {
        self.handle
            .seek(SeekFrom::Start(0))
            .await
            .expect("Failed to reset file stream");
    }

    pub fn get_index(&self) -> &OffsetIndex {
        if !self.index.init {
            warn!("Attempting to use an uninitialized offset index on MzMLReaderType")
        }
        &self.index
    }

    pub fn set_index(&mut self, index: OffsetIndex) {
        self.index = index
    }
}

impl<
        R: io::AsyncRead + io::AsyncSeek + io::AsyncSeekExt + Unpin + Send,
        C: CentroidLike + From<CentroidPeak> + Send + Sync,
        D: DeconvolutedCentroidLike + From<DeconvolutedPeak> + Send + Sync,
    > AsyncSpectrumSource<C, D, MultiLayerSpectrum<C, D>> for MGFReaderType<R, C, D>
{
    async fn reset(&mut self) {
        self.reset().await;
    }

    async fn get_spectrum_by_id(&mut self, id: &str) -> Option<MultiLayerSpectrum<C, D>> {
        self.get_spectrum_by_id(id).await
    }

    async fn get_spectrum_by_index(&mut self, index: usize) -> Option<MultiLayerSpectrum<C, D>> {
        self.get_spectrum_by_index(index).await
    }

    fn detail_level(&self) -> &DetailLevel {
        &self.detail_level
    }

    fn set_detail_level(&mut self, detail_level: DetailLevel) {
        self.detail_level = detail_level;
    }

    fn get_index(&self) -> &OffsetIndex {
        &self.index
    }

    fn set_index(&mut self, index: OffsetIndex) {
        self.index = index;
    }

    async fn read_next(&mut self) -> Option<MultiLayerSpectrum<C, D>> {
        self.read_next().await
    }
}

impl<
        R: io::AsyncRead + io::AsyncSeek + io::AsyncSeekExt + Unpin + Send,
        C: CentroidLike + From<CentroidPeak> + Send + Sync,
        D: DeconvolutedCentroidLike + From<DeconvolutedPeak> + Send + Sync,
    > AsyncRandomAccessSpectrumIterator<C, D, MultiLayerSpectrum<C, D>> for MGFReaderType<R, C, D> {

    async fn start_from_id(&mut self, id: &str) -> Result<&mut Self, SpectrumAccessError> {
        let idx = match self._offset_of_id(id) {
            Some(i) => i,
            None => return Err(crate::io::SpectrumAccessError::SpectrumIdNotFound(id.to_string())),
        };

        self.handle.seek(SeekFrom::Start(idx)).await.map_err(|e| SpectrumAccessError::IOError(Some(e)))?;
        Ok(self)
    }

    async fn start_from_index(&mut self, index: usize) -> Result<&mut Self, SpectrumAccessError> {
        let idx = match self._offset_of_index(index) {
            Some(i) => i,
            None => return Err(crate::io::SpectrumAccessError::SpectrumIndexNotFound(index)),
        };

        self.handle.seek(SeekFrom::Start(idx)).await.map_err(|e| SpectrumAccessError::IOError(Some(e)))?;
        Ok(self)
    }

    async fn start_from_time(&mut self, time: f64) -> Result<&mut Self, SpectrumAccessError> {
        let idx = match self._offset_of_time(time).await {
            Some(i) => i,
            None => return Err(crate::io::SpectrumAccessError::SpectrumNotFound),
        };

        self.handle.seek(SeekFrom::Start(idx)).await.map_err(|e| SpectrumAccessError::IOError(Some(e)))?;
        Ok(self)
    }
}
