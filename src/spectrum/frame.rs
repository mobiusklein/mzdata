#![allow(unused)]
use std::{borrow::Cow, convert::TryFrom, marker::PhantomData, mem};

use mzpeaks::{
    coordinate::{IonMobility, Mass, MZ},
    feature::{ChargedFeature, Feature, FeatureLike},
    feature_map::FeatureMap,
};

use super::{
    bindata::{
        ArrayRetrievalError, ArraysAvailable, BinaryArrayMap3D, BinaryCompressionType,
        BinaryDataArrayType, BuildArrayMap3DFrom, BuildFromArrayMap3D,
    },
    Acquisition, ArrayType, BinaryArrayMap, CentroidPeakAdapting, DeconvolutedPeakAdapting,
    Precursor, ScanPolarity, SignalContinuity, SpectrumConversionError, SpectrumDescription,
};
use super::{scan_properties::SCAN_TITLE, MultiLayerSpectrum};
use crate::{
    io::IonMobilityFrameGrouping,
    params::{ParamDescribed, ParamList},
};
use crate::{prelude::*, RawSpectrum};

/// Represent an owned representation of one the kinds of feature data that a [`IonMobilityFrameLike`](crate::spectrum::IonMobilityFrameLike) instance
/// might otherwise carry.
///
/// # See also
/// [`RefFeatureDataLevel`] for borrowed version.
#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum FeatureDataLevel<
    C: FeatureLike<MZ, IonMobility> = Feature<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge = ChargedFeature<Mass, IonMobility>,
> {
    Missing,
    RawData(BinaryArrayMap3D),
    Centroid(FeatureMap<MZ, IonMobility, C>),
    Deconvoluted(FeatureMap<Mass, IonMobility, D>),
}

/// An variant for dispatching to different strategies of computing
/// common statistics of different levels of feature data.
#[derive(Debug)]
pub enum RefFeatureDataLevel<
    'a,
    C: FeatureLike<MZ, IonMobility> = Feature<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge = ChargedFeature<Mass, IonMobility>,
> {
    Missing,
    RawData(&'a BinaryArrayMap3D),
    Centroid(&'a FeatureMap<MZ, IonMobility, C>),
    Deconvoluted(&'a FeatureMap<Mass, IonMobility, D>),
}

/// The set of descriptive metadata that give context for how an ion mobility frame was acquired
/// within a particular run. This forms the basis for a large portion of the [`IonMobilityFrameLike`]
/// trait.
///
/// This is the equivalent of the [`SpectrumDescription`] type.
#[derive(Debug, Default, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct IonMobilityFrameDescription {
    /// The spectrum's native identifier
    pub id: String,

    /// The ordinal sequence number for the spectrum
    pub index: usize,

    /// The degree of exponentiation of the spectrum, e.g MS1, MS2, MS3, etc
    pub ms_level: u8,

    /// The spectrum is in positive or negative mode
    pub polarity: ScanPolarity,

    /// The spectrum's main representation is as a peak list or a continuous
    /// profile
    pub signal_continuity: SignalContinuity,

    /// A set of controlled or uncontrolled descriptors of the spectrum not already
    /// covered by fields
    pub params: ParamList,

    /// A description of how the spectrum was acquired including time, scan windows, and more
    pub acquisition: Acquisition,

    /// The parent ion or ions and their isolation and activation description
    pub precursor: Option<Precursor>,
}

impl_param_described!(IonMobilityFrameDescription);

impl IonMobilityFrameDescription {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        id: String,
        index: usize,
        ms_level: u8,
        polarity: ScanPolarity,
        signal_continuity: SignalContinuity,
        params: ParamList,
        acquisition: Acquisition,
        precursor: Option<Precursor>,
    ) -> Self {
        Self {
            id,
            index,
            ms_level,
            polarity,
            signal_continuity,
            params,
            acquisition,
            precursor,
        }
    }

    pub fn title(&self) -> Option<Cow<'_, str>> {
        self.get_param_by_curie(&SCAN_TITLE).map(|p| p.as_str())
    }
}

impl From<SpectrumDescription> for IonMobilityFrameDescription {
    fn from(value: SpectrumDescription) -> Self {
        Self::new(
            value.id,
            value.index,
            value.ms_level,
            value.polarity,
            value.signal_continuity,
            value.params,
            value.acquisition,
            value.precursor,
        )
    }
}

impl From<IonMobilityFrameDescription> for SpectrumDescription {
    fn from(value: IonMobilityFrameDescription) -> Self {
        Self::new(
            value.id,
            value.index,
            value.ms_level,
            value.polarity,
            value.signal_continuity,
            value.params,
            value.acquisition,
            value.precursor,
        )
    }
}

/// A trait for providing a uniform delegated access to spectrum metadata
pub trait IonMobilityFrameLike<
    C: FeatureLike<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge,
>: ParamDescribed
{
    /// The method to access the spectrum description itself, which supplies
    /// the data for most other methods on this trait.
    fn description(&self) -> &IonMobilityFrameDescription;

    /// The method to access the spectrum descript itself, mutably.
    fn description_mut(&mut self) -> &mut IonMobilityFrameDescription;

    /// Access the acquisition information for this spectrum.
    #[inline]
    fn acquisition(&self) -> &Acquisition {
        &self.description().acquisition
    }

    /// Access the precursor information, if it exists.
    #[inline]
    fn precursor(&self) -> Option<&Precursor> {
        let desc = self.description();
        if let Some(precursor) = &desc.precursor {
            Some(precursor)
        } else {
            None
        }
    }

    /// Iterate over all precursors of the spectrum
    fn precursor_iter(&self) -> impl Iterator<Item = &Precursor> {
        let desc = self.description();
        desc.precursor.iter()
    }

    /// Mutably access the precursor information, if it exists
    fn precursor_mut(&mut self) -> Option<&mut Precursor> {
        let desc = self.description_mut();
        if let Some(precursor) = desc.precursor.as_mut() {
            Some(precursor)
        } else {
            None
        }
    }

    /// Iterate over all precursors of the spectrum mutably
    fn precursor_iter_mut(&mut self) -> impl Iterator<Item = &mut Precursor> {
        let desc = self.description_mut();
        desc.precursor.iter_mut()
    }

    /// A shortcut method to retrieve the scan start time of a spectrum
    #[inline]
    fn start_time(&self) -> f64 {
        let acq = self.acquisition();
        match acq.scans.first() {
            Some(evt) => evt.start_time,
            None => 0.0,
        }
    }

    /// Access the MS exponentiation level
    #[inline]
    fn ms_level(&self) -> u8 {
        self.description().ms_level
    }

    /// Access the native ID string for the spectrum
    #[inline]
    fn id(&self) -> &str {
        &self.description().id
    }

    /// Access the index of the spectrum in the source file
    #[inline]
    fn index(&self) -> usize {
        self.description().index
    }

    /// Access a description of how raw the signal is, whether a
    /// profile spectrum is available or only centroids are present.
    #[inline]
    fn signal_continuity(&self) -> SignalContinuity {
        self.description().signal_continuity
    }

    #[inline]
    fn polarity(&self) -> ScanPolarity {
        self.description().polarity
    }

    fn raw_arrays(&'_ self) -> Option<&'_ BinaryArrayMap3D>;

    fn features(&self) -> RefFeatureDataLevel<C, D>;

    fn into_features_and_parts(self) -> (FeatureDataLevel<C, D>, IonMobilityFrameDescription);
}

/// An ion mobility dimension-equivalent of [`MultiLayerSpectrum`].
///
/// An ion mobility frame is just like a mass spectrum, except the signal
/// is spread across multiple ion mobility "points". This means that instead
/// producing "points" in the mass coordinate system, this produces traces-over-time
/// or "features".
///
/// An ion mobility frames can usually be inter-converted with spectra.
#[derive(Debug, Default, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MultiLayerIonMobilityFrame<
    C: FeatureLike<MZ, IonMobility> = Feature<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge = ChargedFeature<Mass, IonMobility>,
> {
    /// The raw data arrays, split by ion mobility point
    pub arrays: Option<BinaryArrayMap3D>,
    /// The (potentially absent) m/z coordinate feature map
    pub features: Option<FeatureMap<MZ, IonMobility, C>>,
    /// The (potentially absent) neutral mass coordinate feature map
    pub deconvoluted_features: Option<FeatureMap<Mass, IonMobility, D>>,
    pub description: IonMobilityFrameDescription,
}

impl<
        C: FeatureLike<MZ, IonMobility> + BuildFromArrayMap3D,
        D: FeatureLike<Mass, IonMobility> + KnownCharge + BuildFromArrayMap3D,
    > MultiLayerIonMobilityFrame<C, D>
{
    pub fn try_build_features(
        &mut self,
    ) -> Result<RefFeatureDataLevel<'_, C, D>, SpectrumConversionError> {
        if matches!(self.signal_continuity(), SignalContinuity::Centroid) {
            if let Some(arrays) = self.arrays.as_ref() {
                if let ArraysAvailable::Ok = D::has_arrays_3d_for(arrays) {
                    let peaks = FeatureMap::new(D::try_from_arrays_3d(arrays)?);
                    self.deconvoluted_features = Some(peaks);
                    return Ok(self.features());
                }

                if let ArraysAvailable::Ok = C::has_arrays_3d_for(arrays) {
                    let peaks = FeatureMap::new(C::try_from_arrays_3d(arrays)?);
                    self.features = Some(peaks);
                    return Ok(self.features());
                }
                return Ok(RefFeatureDataLevel::Missing);
            } else {
                return Ok(self.features());
            }
        } else {
            return Ok(self.features());
        }
    }
}

impl<C: FeatureLike<MZ, IonMobility>, D: FeatureLike<Mass, IonMobility> + KnownCharge>
    ParamDescribed for MultiLayerIonMobilityFrame<C, D>
{
    fn params(&self) -> &[crate::Param] {
        &self.description().params
    }

    fn params_mut(&mut self) -> &mut ParamList {
        &mut self.description_mut().params
    }
}

impl<C: FeatureLike<MZ, IonMobility>, D: FeatureLike<Mass, IonMobility> + KnownCharge>
    MultiLayerIonMobilityFrame<C, D>
{
    pub fn new(
        arrays: Option<BinaryArrayMap3D>,
        features: Option<FeatureMap<MZ, IonMobility, C>>,
        deconvoluted_features: Option<FeatureMap<Mass, IonMobility, D>>,
        description: IonMobilityFrameDescription,
    ) -> Self {
        Self {
            arrays,
            features,
            deconvoluted_features,
            description,
        }
    }
}

impl<C: FeatureLike<MZ, IonMobility>, D: FeatureLike<Mass, IonMobility> + KnownCharge>
    IonMobilityFrameLike<C, D> for MultiLayerIonMobilityFrame<C, D>
{
    fn description(&self) -> &IonMobilityFrameDescription {
        &self.description
    }

    fn description_mut(&mut self) -> &mut IonMobilityFrameDescription {
        &mut self.description
    }

    fn raw_arrays(&'_ self) -> Option<&'_ BinaryArrayMap3D> {
        self.arrays.as_ref()
    }

    fn features(&self) -> RefFeatureDataLevel<C, D> {
        let state = if let Some(d) = self.deconvoluted_features.as_ref() {
            RefFeatureDataLevel::Deconvoluted(d)
        } else if let Some(c) = self.features.as_ref() {
            RefFeatureDataLevel::Centroid(c)
        } else if let Some(arrays) = self.arrays.as_ref() {
            RefFeatureDataLevel::RawData(arrays)
        } else {
            RefFeatureDataLevel::Missing
        };
        state
    }

    fn into_features_and_parts(self) -> (FeatureDataLevel<C, D>, IonMobilityFrameDescription) {
        let state = if let Some(d) = self.deconvoluted_features {
            FeatureDataLevel::Deconvoluted(d)
        } else if let Some(c) = self.features {
            FeatureDataLevel::Centroid(c)
        } else if let Some(arrays) = self.arrays {
            FeatureDataLevel::RawData(arrays)
        } else {
            FeatureDataLevel::Missing
        };

        (state, self.description)
    }
}

impl<CFeat: FeatureLike<MZ, IonMobility>, DFeat: FeatureLike<Mass, IonMobility> + KnownCharge>
    TryFrom<RawSpectrum> for MultiLayerIonMobilityFrame<CFeat, DFeat>
{
    type Error = ArrayRetrievalError;

    fn try_from(value: RawSpectrum) -> Result<Self, Self::Error> {
        let arrays = BinaryArrayMap3D::try_from(value.arrays)?;
        let descr = IonMobilityFrameDescription::from(value.description);

        Ok(MultiLayerIonMobilityFrame::new(
            Some(arrays),
            None,
            None,
            descr,
        ))
    }
}

impl<
        CFeat: FeatureLike<MZ, IonMobility>,
        DFeat: FeatureLike<Mass, IonMobility> + KnownCharge,
        CPeak: CentroidPeakAdapting,
        DPeak: DeconvolutedPeakAdapting,
    > TryFrom<MultiLayerSpectrum<CPeak, DPeak>> for MultiLayerIonMobilityFrame<CFeat, DFeat>
{
    type Error = ArrayRetrievalError;

    fn try_from(value: MultiLayerSpectrum<CPeak, DPeak>) -> Result<Self, Self::Error> {
        let arrays = match value.arrays {
            Some(arrays) => BinaryArrayMap3D::try_from(arrays)?,
            None => return Err(ArrayRetrievalError::NotFound(ArrayType::IonMobilityArray)),
        };

        let descr = IonMobilityFrameDescription::from(value.description);

        Ok(MultiLayerIonMobilityFrame::new(
            Some(arrays),
            None,
            None,
            descr,
        ))
    }
}

impl<
        C: FeatureLike<MZ, IonMobility> + BuildArrayMap3DFrom,
        D: FeatureLike<Mass, IonMobility> + KnownCharge + BuildArrayMap3DFrom,
    > From<MultiLayerIonMobilityFrame<C, D>> for RawSpectrum
{
    fn from(value: MultiLayerIonMobilityFrame<C, D>) -> Self {
        let arrays = if let Some(d) = value.deconvoluted_features {
            BuildArrayMap3DFrom::as_arrays_3d(&d[..]).unstack().unwrap()
        } else if let Some(c) = value.features {
            BuildArrayMap3DFrom::as_arrays_3d(&c[..]).unstack().unwrap()
        } else if let Some(arrays) = value.arrays {
            arrays.unstack().unwrap()
        } else {
            BinaryArrayMap::default()
        };

        let descr = value.description.into();

        RawSpectrum::new(descr, arrays)
    }
}

impl<
        CFeat: FeatureLike<MZ, IonMobility> + BuildArrayMap3DFrom,
        DFeat: FeatureLike<Mass, IonMobility> + KnownCharge + BuildArrayMap3DFrom,
        CPeak: CentroidPeakAdapting + BuildFromArrayMap,
        DPeak: DeconvolutedPeakAdapting + BuildFromArrayMap,
    > From<MultiLayerIonMobilityFrame<CFeat, DFeat>> for MultiLayerSpectrum<CPeak, DPeak>
{
    fn from(value: MultiLayerIonMobilityFrame<CFeat, DFeat>) -> Self {
        let raw: RawSpectrum = value.into();
        raw.into()
    }
}

#[cfg(feature = "mzsignal")]
mod mzsignal_impl {
    use super::*;

    use mzpeaks::{feature::Feature, peak_set::PeakSetVec, CentroidPeak, Tolerance};
    use mzsignal::{
        feature_mapping::{FeatureExtracterType, FeatureGraphBuilder, MapState, PeakMapState},
        peak_picker::PeakPicker,
        FittedPeak, PeakFitType,
    };

    impl<
            C: FeatureLike<MZ, IonMobility> + FeatureLikeMut<MZ, IonMobility> + Default + Clone,
            D: FeatureLike<Mass, IonMobility> + KnownCharge,
        > MultiLayerIonMobilityFrame<C, D>
    {
        pub fn extract_features_with<
            E: MapState<CentroidPeak, MZ, IonMobility, FeatureType = C> + Default,
        >(
            &mut self,
            error_tolerance: Tolerance,
            min_length: usize,
            maximum_gap_size: f64,
            peak_picker: Option<PeakPicker>,
        ) -> Result<(), ArrayRetrievalError> {
            if let Some(arrays) = self.arrays.as_ref() {
                let features = match self.signal_continuity() {
                    SignalContinuity::Unknown => {
                        panic!("Cannot extract feature, unknown signal continuities")
                    }
                    SignalContinuity::Centroid => {
                        let mut extractor: FeatureExtracterType<E, _, _, IonMobility> = arrays
                            .iter()
                            .map(
                                |(time, submap)| -> Result<
                                    (f64, PeakSetVec<CentroidPeak, MZ>),
                                    ArrayRetrievalError,
                                > {
                                    let mzs = submap.mzs()?;
                                    let intens = submap.intensities()?;
                                    let peaks = mzs
                                        .iter()
                                        .zip(intens.iter())
                                        .map(|(mz, inten)| CentroidPeak::new(*mz, *inten, 0))
                                        .collect();
                                    Ok((time as f64, peaks))
                                },
                            )
                            .flatten()
                            .collect();
                        extractor.extract_features(error_tolerance, min_length, maximum_gap_size)
                    }
                    SignalContinuity::Profile => {
                        let peak_picker = peak_picker.unwrap_or_else(|| {
                            PeakPicker::new(1.0, 1.0, 1.0, PeakFitType::Quadratic)
                        });
                        let mut extractor: FeatureExtracterType<E, _, _, IonMobility> = arrays
                            .iter()
                            .map(
                                |(time, submap)| -> Result<
                                    (f64, PeakSetVec<CentroidPeak, MZ>),
                                    ArrayRetrievalError,
                                > {
                                    let mzs = submap.mzs()?;
                                    let intens = submap.intensities()?;
                                    let mut peaks = Vec::new();
                                    peak_picker
                                        .discover_peaks(&mzs, &intens, &mut peaks)
                                        .unwrap();
                                    let peaks =
                                        peaks.into_iter().map(|p| p.as_centroid()).collect();
                                    Ok((time as f64, peaks))
                                },
                            )
                            .flatten()
                            .collect();
                        extractor.extract_features(error_tolerance, min_length, maximum_gap_size)
                    }
                };
                self.features = Some(features);
            }
            Ok(())
        }
    }

    impl<D: FeatureLike<Mass, IonMobility> + KnownCharge>
        MultiLayerIonMobilityFrame<Feature<MZ, IonMobility>, D>
    {
        pub fn extract_features_simple(
            &mut self,
            error_tolerance: Tolerance,
            min_length: usize,
            maximum_gap_size: f64,
            peak_picker: Option<PeakPicker>,
        ) -> Result<(), ArrayRetrievalError> {
            self.extract_features_with::<PeakMapState<CentroidPeak, MZ>>(
                error_tolerance,
                min_length,
                maximum_gap_size,
                peak_picker,
            )
        }
    }
}

#[cfg(test)]
mod test {
    use std::{fs, io};

    use super::*;

    use crate::io::Generic3DIonMobilityFrameSource;
    use crate::prelude::*;

    macro_rules! assert_is_close {
        ($t1:expr, $t2:expr, $tol:expr, $label:literal) => {
            assert!(
                ($t1 - $t2).abs() < $tol,
                "Observed {} {}, expected {}, difference {}",
                $label,
                $t1,
                $t2,
                $t1 - $t2,
            );
        };
        ($t1:expr, $t2:expr, $tol:expr, $label:literal, $obj:ident) => {
            assert!(
                ($t1 - $t2).abs() < $tol,
                "Observed {} {}, expected {}, difference {} from {:?}",
                $label,
                $t1,
                $t2,
                $t1 - $t2,
                $obj
            );
        };
    }

    #[test]
    fn test_loader_conversion() -> io::Result<()> {
        let fh = fs::File::open("./test/data/20200204_BU_8B8egg_1ug_uL_7charges_60_min_Slot2-11_1_244.mzML.gz")?;
        let handle = crate::io::compression::RestartableGzDecoder::new(io::BufReader::new(fh));
        let mzml_reader = crate::MzMLReader::new_indexed(handle);
        let mut frame_reader = mzml_reader.try_into_frame_source::<Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>>().unwrap();
        let frame = frame_reader.get_frame_by_id("merged=42926 frame=9728 scanStart=1 scanEnd=705").unwrap();
        assert_eq!(frame.ms_level(), 1);
        Ok(())
    }

    #[test]
    fn test_reader_wrapper_iter() -> io::Result<()> {
        let group = crate::mz_read!("./test/data/20200204_BU_8B8egg_1ug_uL_7charges_60_min_Slot2-11_1_244.mzML.gz".as_ref(), reader => {
            let mut wrapper: Generic3DIonMobilityFrameSource<_, _, _, Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>> = Generic3DIonMobilityFrameSource::new(reader);
            let mut group_iter = wrapper.into_groups();
            group_iter.start_from_id("merged=42926 frame=9728 scanStart=1 scanEnd=705")?;
            group_iter.next().unwrap()
        })?;

        assert_eq!(group.lowest_ms_level().unwrap(), 1);
        assert_eq!(group.highest_ms_level().unwrap(), 2);

        Ok(())
    }

    #[test]
    fn test_reader_wrapper_extract() -> io::Result<()> {
        crate::mz_read!("./test/data/20200204_BU_8B8egg_1ug_uL_7charges_60_min_Slot2-11_1_244.mzML.gz".as_ref(), reader => {
            let mut wrapper: Generic3DIonMobilityFrameSource<_, _, _, Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>> = Generic3DIonMobilityFrameSource::new(reader);
            let mut frame = wrapper.get_frame_by_id("merged=42926 frame=9728 scanStart=1 scanEnd=705").unwrap();
            assert_eq!(frame.ms_level(), 1);
            assert_eq!(frame.polarity(), ScanPolarity::Positive);

            #[cfg(feature = "mzsignal")]
            {
                frame.extract_features_simple(Tolerance::PPM(15.0), 2, 0.1, None)?;

                let fm = frame.features.as_ref().unwrap();
                let query = 1456.95;
                let hits = fm.all_features_for(query, Tolerance::PPM(15.0));
                assert_eq!(hits.len(), 2);

                let feat = &hits[0];

                if false {
                    let mut writer = io::BufWriter::new(std::fs::File::create("tmp/peaks_over_time.txt")?);
                    let arrays = frame.arrays.as_ref().unwrap();
                    for (im, arrs) in arrays.ion_mobility_dimension.iter().zip(arrays.arrays.iter()) {
                        let mzs = arrs.mzs()?;
                        let intensities = arrs.intensities()?;
                        for (mz, int) in mzs.iter().zip(intensities.iter()) {
                            writeln!(writer, "{mz}\t{int}\t{im}")?;
                        }
                    }
                }

                if false {
                    let mut writer = io::BufWriter::new(std::fs::File::create("tmp/features_graph.txt")?);
                    writer.write_all(b"feature_id\tmz\trt\tintensity\n")?;
                    for (i, f) in fm.iter().enumerate() {
                        for (mz, rt, inten) in f.iter() {
                            writer.write_all(format!("{i}\t{mz}\t{rt}\t{inten}\n").as_bytes())?;
                        }
                    }
                }

                eprintln!("{:?}\t{:?}\t{:?}\t{}", feat.start_time(), feat.end_time(), feat.apex_time(), feat.len());

                assert_is_close!(feat.start_time().unwrap(), 0.947379971788502, 1e-3, "start_time");
                assert_is_close!(feat.end_time().unwrap(), 1.2638564827665548, 1e-3, "end_time");
                assert_is_close!(feat.apex_time().unwrap(), 1.212666, 1e-3, "apex_time");
                assert_eq!(feat.len(), 98);
            }
        });
        Ok(())
    }
}
