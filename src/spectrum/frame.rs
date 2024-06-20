#![allow(unused)]
use std::{borrow::Cow, convert::TryFrom};

use mzpeaks::{
    coordinate::{IonMobility, Mass, MZ},
    feature::{ChargedFeature, Feature, FeatureLike},
    feature_map::FeatureMap,
};

use super::{
    bindata::{ArrayRetrievalError, BinaryArrayMap3D, BinaryCompressionType, BinaryDataArrayType},
    Acquisition, ArrayType, CentroidPeakAdapting, DeconvolutedPeakAdapting, Precursor,
    ScanPolarity, SignalContinuity, SpectrumDescription,
};
use super::{scan_properties::SCAN_TITLE, MultiLayerSpectrum};
use crate::params::{ParamDescribed, ParamList};
use crate::prelude::*;

#[derive(Debug, Default, Clone)]
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

/// A trait for providing a uniform delegated access to spectrum metadata
pub trait IonMobilityFrameLike<
    C: FeatureLike<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge,
>
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

    #[inline]
    fn params(&self) -> &ParamList {
        &self.description().params
    }
}

#[derive(Debug, Default, Clone)]
pub struct MultiLayerIonMobilityFrame<
    C: FeatureLike<MZ, IonMobility> = Feature<MZ, IonMobility>,
    D: FeatureLike<Mass, IonMobility> + KnownCharge = ChargedFeature<Mass, IonMobility>,
> {
    pub arrays: Option<BinaryArrayMap3D>,
    pub features: Option<FeatureMap<MZ, IonMobility, C>>,
    pub deconvoluted_features: Option<FeatureMap<Mass, IonMobility, D>>,
    description: IonMobilityFrameDescription,
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
            M: FeatureGraphBuilder<MZ, IonMobility, C>,
            E: MapState<CentroidPeak, MZ, IonMobility, FeatureType = C, FeatureMergerType = M>
                + Default,
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
            self.extract_features_with::<
            <PeakMapState<CentroidPeak, MZ> as MapState<
                CentroidPeak,
                MZ,
                IonMobility,
            >>::FeatureMergerType,
            PeakMapState<CentroidPeak, MZ>
            >(
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
    use std::io;

    use super::*;

    use crate::io::Generic3DIonMobilityFrameSource;
    use crate::prelude::*;

    #[test]
    fn test_reader_wrapper() -> io::Result<()> {
        crate::mz_read!("./test/data/20200204_BU_8B8egg_1ug_uL_7charges_60_min_Slot2-11_1_244.mzML.gz".as_ref(), reader => {
            let mut wrapper: Generic3DIonMobilityFrameSource<_, _, _, Feature<MZ, IonMobility>, ChargedFeature<Mass, IonMobility>> = Generic3DIonMobilityFrameSource::new(reader);
            let mut frame = wrapper.get_frame_by_id("merged=42926 frame=9728 scanStart=1 scanEnd=705").unwrap();
            assert_eq!(frame.ms_level(), 1);
            assert_eq!(frame.polarity(), ScanPolarity::Positive);

            #[cfg(feature = "mzsignal")]
            {
                frame.extract_features_simple(Tolerance::PPM(10.0), 3, 0.1, None)?;

                let fm = frame.features.as_ref().unwrap();
                let query = 1456.95;
                let hits = fm.all_features_for(query, Tolerance::PPM(10.0));
                assert_eq!(hits.len(), 1);

                let feat = &hits[0];

                assert!((feat.start_time().unwrap() - 1.173088550567627).abs() < 1e-3);
                assert!((feat.end_time().unwrap() - 1.2521668672561646).abs() < 1e-3);
                assert!((feat.apex_time().unwrap() - 1.2126655578613281).abs() < 1e-3);
                assert_eq!(feat.len(), 49);
            }
        });
        Ok(())
    }
}
