use mzdata::prelude::*;
use mzdata::spectrum::SignalContinuity;
use mzpeaks::Tolerance;
use std::io;

fn main() -> io::Result<()> {
    mzdata::mz_read!("./test/data/small.mzML".as_ref(), reader =>
    {
            let mut ms1_count = 0;
            let mut msn_count = 0;
            for spectrum in reader {
                if spectrum.ms_level()==1 {
                        ms1_count+=1;
                }else {
                        msn_count+=1;
                    }
                println!("Scan {} => BP {}",spectrum.id(),spectrum.peaks().base_peak().mz);
                if spectrum.signal_continuity()<SignalContinuity::Profile {
                    let peak_picked = spectrum.into_centroid().unwrap();
                        println!("Matches for 579.155: {:?}",peak_picked.peaks.all_peaks_for(579.155,Tolerance::Da(0.02)));
                    }
                }
                println!("MS1 Count: {}\nMSn Count: {}",ms1_count,msn_count);
                assert_eq!(ms1_count,14);
                assert_eq!(msn_count,34);
            }
        )
}
