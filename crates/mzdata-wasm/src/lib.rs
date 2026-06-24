use log;
use wasm_bindgen::prelude::*;
use wasm_logger;

mod binds;
mod mem_reader;
mod mem_writer;
mod utils;
mod blobio;
// mod webio;
// mod worker_reader;
pub mod asyncio;

pub use asyncio::{test_reader, WebReaderAsyncRead};
pub use blobio::test_reader_blob;
pub use binds::*;
pub use mem_reader::{MemWebIMMZReader as WebIMMZReader, MemWebMZReader as WebMZReader};

pub fn set_panic_hook() {
    // When the `console_error_panic_hook` feature is enabled, we can call the
    // `set_panic_hook` function at least once during initialization, and then
    // we will get better error messages if our code ever panics.
    //
    // For more details see
    // https://github.com/rustwasm/console_error_panic_hook#readme
    console_error_panic_hook::set_once();
}

#[wasm_bindgen(start)]
fn start() {
    set_panic_hook();
    wasm_logger::init(wasm_logger::Config::new(log::Level::Debug));
}

/*
#[wasm_bindgen]
pub fn read_file(event: &Event) {
    event.prevent_default();
    event.stop_propagation();

    let source = event
        .target()
        .unwrap()
        .dyn_into::<HtmlInputElement>()
        .unwrap();
    let files = source.files().unwrap();
    if files.length() > 0 {
        let file = files.item(0).unwrap();
        let file_reader = FileReader::new().unwrap();
        file_reader.read_as_array_buffer(&file).unwrap();
        let onload = Closure::wrap(Box::new(move |event: Event| {
            log::info!("Reading File");
            let file_reader: FileReader = event.target().unwrap().dyn_into().unwrap();
            let js_buff = file_reader.result().unwrap();
            let js_bytes = Uint8Array::new(&js_buff);
            let mut rs_buff: Vec<u8> = vec![0; js_bytes.length() as usize];
            js_bytes.copy_to(&mut rs_buff);

            let rs_buff = io::Cursor::new(rs_buff);
            let mut mz_reader = MzMLReader::new_indexed(rs_buff);
            let mut hexnac_intensity = Vec::new();
            let queries = vec![("HexNAc", 203.0793)];

            log::info!("Reader with {} scans", mz_reader.len());
            let (ms1_count, ms2_count) = mz_reader.iter().enumerate().fold(
                (0usize, 0usize),
                |(ms1_count, ms2_count), (i, mut scan)| {
                    if scan.ms_level() > 1 {
                        match scan.signal_continuity() {
                            mzdata::spectrum::SignalContinuity::Unknown => {
                                panic!("Don't know how to handle spectrum")
                            }
                            mzdata::spectrum::SignalContinuity::Centroid => {
                                // We have a pre-deconvoluted neutral mass peak list
                                let matched_peaks =
                                    if let Ok(peaks) = scan.try_build_deconvoluted_centroids() {
                                        let hits = queries
                                            .iter()
                                            .map(|(q, mass)| {
                                                (
                                                    q,
                                                    peaks
                                                        .all_peaks_for(*mass, Tolerance::PPM(10.0))
                                                        .to_vec(),
                                                )
                                            })
                                            .collect::<Vec<_>>();
                                        hits
                                    }
                                    // We have a centroided peak list but still in the m/z dimension
                                    else if let Ok(peaks) = scan.try_build_centroids() {
                                        if i % 100 == 0 {
                                            log::info!("Deconvolving centroid spectrum {}", i);
                                        }
                                        let deconv_peaks: MassPeakSetType<_> =
                                        mzdeisotope::deconvolute_peaks(
                                            peaks.clone(),
                                            mzdeisotope::isotopic_model::IsotopicModels::Glycopeptide,
                                            Tolerance::PPM(20.0),
                                            (1, 8),
                                            mzdeisotope::scorer::MSDeconvScorer::new(0.02),
                                            mzdeisotope::scorer::MaximizingFitFilter::new(5.0),
                                            1,
                                            mzdeisotope::isotopic_model::IsotopicPatternParams::default(
                                            ),
                                            true,
                                        )
                                        .unwrap()
                                        .into_iter()
                                        .map(|p| p.as_centroid())
                                        .collect();

                                        let hits = queries
                                            .iter()
                                            .map(|(q, mass)| {
                                                (
                                                    q,
                                                    deconv_peaks
                                                        .all_peaks_for(*mass, Tolerance::PPM(10.0))
                                                        .to_vec(),
                                                )
                                            })
                                            .collect::<Vec<_>>();
                                        scan.deconvoluted_peaks = Some(deconv_peaks);
                                        hits
                                    } else {
                                        Vec::new()
                                    };
                                hexnac_intensity.push(matched_peaks);
                            }
                            mzdata::spectrum::SignalContinuity::Profile => {
                                if i % 100 == 0 {
                                    log::info!("Processing profile spectrum {}", i);
                                }
                                scan.pick_peaks(1.0)
                                    .unwrap();
                                let deconv_peaks: MassPeakSetType<_> =
                                    mzdeisotope::deconvolute_peaks(
                                        scan.peaks.as_ref().unwrap().clone(),
                                        mzdeisotope::isotopic_model::IsotopicModels::Glycopeptide,
                                        Tolerance::PPM(20.0),
                                        (1, 8),
                                        mzdeisotope::scorer::MSDeconvScorer::new(0.02),
                                        mzdeisotope::scorer::MaximizingFitFilter::new(5.0),
                                        1,
                                        mzdeisotope::isotopic_model::IsotopicPatternParams::default(),
                                        true,
                                    )
                                    .unwrap()
                                    .into_iter()
                                    .map(|p| p.as_centroid())
                                    .collect();
                                let hits = queries
                                    .iter()
                                    .map(|(q, mass)| {
                                        (
                                            q,
                                            deconv_peaks
                                                .all_peaks_for(*mass, Tolerance::PPM(10.0))
                                                .to_vec(),
                                        )
                                    })
                                    .collect::<Vec<_>>();
                                hexnac_intensity.push(hits);
                                scan.deconvoluted_peaks = Some(deconv_peaks);
                            }
                        }
                    }

                    if i % 1000 == 0 {
                        log::info!("{} -> {:?}", scan.id(), scan.signal_continuity());
                    };
                    if scan.ms_level() == 1 {
                        (ms1_count + 1, ms2_count)
                    } else {
                        (ms1_count, ms2_count + 1)
                    }
                },
            );
            log::info!(
                "Reader with {} MS1 scans and {} MSn scans",
                ms1_count,
                ms2_count
            );

            // Do something with `hexnac_intensity` here
            let total_intensity: f32 = hexnac_intensity.iter().map(|peak_matches| {
                peak_matches.iter().map(|(_query, peaks)| {
                    peaks.iter().map(|p| p.intensity()).sum::<f32>()
                }).sum::<f32>()
            }).sum();
            log::info!("Matched a total intensity of {total_intensity} HexNAc oxonium ions");

            log::info!("Done reading");
        }) as Box<dyn FnMut(_)>);

        file_reader.set_onload(Some(onload.as_ref().unchecked_ref()));
        onload.forget();
    }

}
 */
