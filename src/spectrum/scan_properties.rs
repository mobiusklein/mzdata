use mzdata_spectrum::{CentroidPeakAdapting, DeconvolutedPeakAdapting, SpectrumLike};
pub use mzdata_spectrum::{
    Acquisition, Activation, AsPrecursorCollection, ChromatogramType, ChromatogramDescription,
    IonProperties, IsolationWindow, IsolationWindowState, Precursor, PrecursorSelection, Product,
    ScanPolarity, ScanCombination, ScanEvent, ScanWindow, SelectedIon, SignalContinuity, IonMobilityMeasure,
    SpectrumDescription
};
use crate::io::traits::SpectrumSource;


pub trait PrecursorRetrieval : PrecursorSelection {
    /// Given a [`SpectrumSource`] object, look up the precursor scan in it.
    /// This is useful when examining the area *around* where the precursor
    /// ion was or to obtain a snapshot of the retention time when the spectrum
    /// was scheduled.
    fn precursor_spectrum<C, D, S, R>(&self, source: &mut R) -> Option<S>
    where
        C: CentroidPeakAdapting,
        D: DeconvolutedPeakAdapting,
        S: SpectrumLike<C, D>,
        R: SpectrumSource<C, D, S>,
    {
        match self.precursor_id().as_ref() {
            Some(id) => source.get_spectrum_by_id(id),
            None => None,
        }
    }

    /// Given a [`SpectrumSource`] object, look up the product scan in it.
    /// This is rarely needed unless you have manually separated [`Precursor`]
    /// objects from their spectra.
    fn product_spectrum<C, D, S, R>(&self, source: &mut R) -> Option<S>
    where
        C: CentroidPeakAdapting,
        D: DeconvolutedPeakAdapting,
        S: SpectrumLike<C, D>,
        R: SpectrumSource<C, D, S>,
    {
        match self.product_id().as_ref() {
            Some(id) => source.get_spectrum_by_id(id),
            None => None,
        }
    }
}

impl<T: PrecursorSelection + ?Sized> PrecursorRetrieval for T {}


#[cfg(test)]
mod contains_tests {
    use super::*;

    // Regression tests for the isolation/scan-window `contains` membership check. A prior bug used
    // `lower <= point && upper <= point`, which wrongly excluded interior points and wrongly included
    // points above the window; `contains` must be an inclusive `[lower, upper]` test.

    #[test]
    fn isolation_window_contains_is_inclusive_range() {
        let window = IsolationWindow::around(810.789, 1.0); // [809.789, 811.789]
        assert!(window.contains(window.target), "target must be inside");
        assert!(window.contains(window.lower_bound), "lower bound is inclusive");
        assert!(window.contains(window.upper_bound), "upper bound is inclusive");
        assert!(window.contains(810.0_f64), "interior point");
        assert!(!window.contains(809.0_f64), "below the lower bound");
        assert!(!window.contains(812.0_f64), "above the upper bound");
    }

    #[test]
    fn scan_window_contains_is_inclusive_range() {
        let window = ScanWindow::new(200.0, 2000.0);
        assert!(window.contains(1100.0_f64), "interior point");
        assert!(window.contains(200.0_f64), "lower bound is inclusive");
        assert!(window.contains(2000.0_f64), "upper bound is inclusive");
        assert!(!window.contains(199.0_f64), "below the lower bound");
        assert!(!window.contains(2500.0_f64), "above the upper bound (the old bug returned true here)");
    }
}
