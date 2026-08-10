/*! Metadata describing mass spectrometry data files and their contents.
 */

pub use mzdata_meta::*;


#[cfg(test)]
mod test {
    use std::io;
    use crate::prelude::*;
    use crate::params::*;
    use crate::params::ParamCow;

    use super::*;

    #[test]
    fn test_source_file() -> io::Result<()> {
        let sf = SourceFile::from_path("test/data/small.mzML")?;
        assert_eq!(
            MassSpectrometerFileFormatTerm::MzML.to_param(),
            sf.file_format.clone().unwrap()
        );
        Ok(())
    }

    #[test]
    fn test_parser() {
        let ident = NativeSpectrumIdentifierFormatTerm::ThermoNativeIDFormat
            .parse("controllerType=0 controllerNumber=1 scan=25788")
            .unwrap();
        let scan_number = ident.name("scan").unwrap().as_str();
        assert_eq!(scan_number, "25788");
    }

    #[test]
    fn test_meta() {
        assert_eq!(
            NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat.flags(),
            r"frame=(?<frame>\d+) scan=(?<scan>\d+)"
        );
        assert_eq!(
            NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat.parents(),
            [NativeSpectrumIdentifierFormatTerm::NativeSpectrumIdentifierFormat]
        );
        let param: ParamCow<'static> =
            (NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat).into();
        assert_eq!(
            param.curie().unwrap(),
            CURIE::new(ControlledVocabulary::MS, 1002818)
        );
        let param: ParamCow<'static> =
            (&NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat).into();
        assert_eq!(
            param.curie().unwrap(),
            CURIE::new(ControlledVocabulary::MS, 1002818)
        );

        assert_eq!(
            NativeSpectrumIdentifierFormatTerm::from_param(&param),
            Some(NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat)
        );

        let param: Param = (NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat).into();
        assert_eq!(
            param.curie().unwrap(),
            CURIE::new(ControlledVocabulary::MS, 1002818)
        );
        let param: Param = (&NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat).into();
        assert_eq!(
            param.curie().unwrap(),
            CURIE::new(ControlledVocabulary::MS, 1002818)
        );

        assert_eq!(
            NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat.accession(),
            1002818
        );
        assert_eq!(
            NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat.controlled_vocabulary(),
            ControlledVocabulary::MS
        );
        assert_eq!(
            NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat.name(),
            "Bruker TDF nativeID format"
        );

        assert_eq!(
            NativeSpectrumIdentifierFormatTerm::from_name("Bruker TDF nativeID format"),
            Some(NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat)
        );

        assert_eq!(
            NativeSpectrumIdentifierFormatTerm::from_curie(&CURIE::new(
                ControlledVocabulary::MS,
                1002818
            )),
            Some(NativeSpectrumIdentifierFormatTerm::BrukerTDFNativeIDFormat)
        )
    }

    #[test]
    fn test_format() {
        let fmt = NativeSpectrumIdentifierFormatTerm::ThermoNativeIDFormat.format([
            ValueRef::Int(0),
            ValueRef::Int(1),
            ValueRef::Int(25788),
        ]);
        assert_eq!(fmt, "controllerType=0 controllerNumber=1 scan=25788");

        let fmt = NativeSpectrumIdentifierFormatTerm::ThermoNativeIDFormat
            .build()
            .format([ValueRef::Int(0), ValueRef::Int(1)])
            .unwrap_err();
        assert_eq!(
            fmt,
            NativeIDFormatError::IncorrectArgumentNumber {
                term: NativeSpectrumIdentifierFormatTerm::ThermoNativeIDFormat,
                expected: 3,
                received: 2
            }
        );
    }
}
