use quick_xml::Reader;
use quick_xml::events::{Event, BytesStart, BytesText, BytesEnd};

use std::io;

use crate::spectrum::params::{Param, ParamList};
use crate::spectrum::scan_properties::*;
use crate::spectrum::signal::{DataArray, BinaryArrayMap, BinaryCompressionType, BinaryDataArrayType, ArrayType};

pub type Bytes = Vec<u8>;

#[derive(Debug, Clone, Copy)]
pub enum MzMLParserState {
    Start,
    Spectrum,
    SpectrumList,
    BinaryDataArrayList,
    BinaryDataArray,
    Binary,
    ScanList,
    Scan,
    ScanWindowList,
    ScanWindow,
    PrecursorList,
    Precursor,
    IsolationWindow,
    SelectedIonList,
    SelectedIon,
    Activation,
    ParserError,
    SpectrumDone,
}


const BUFFER_SIZE: usize = 10000;


#[derive(Default)]
pub struct MzMLSpectrumBuilder {
    pub params: ParamList,
    pub acquisition: AcquisitionInformation,
    pub precursor: Precursor,
    pub arrays: BinaryArrayMap,
    pub current_array: DataArray,
    pub index: usize,
    pub scan_id: String,
}

pub type ParserResult = Result<MzMLParserState, String>;

impl MzMLSpectrumBuilder {
    pub fn new() -> MzMLSpectrumBuilder {
        MzMLSpectrumBuilder {..Default::default()}
    }

    pub fn reset(&mut self) {
        self.params.clear();
        self.acquisition = AcquisitionInformation::default();
        self.arrays.clear();
        self.current_array.clear();
        self.scan_id.clear();
        self.precursor = Precursor::default();
        self.index = 0;
    }

    fn handle_param<B: io::BufRead>(&self, event: &BytesStart, reader: &Reader<B>) -> Param {
        let mut param = Param::new();
        for attr_parsed in event.attributes() {
            match attr_parsed {
                Ok(attr) => {
                    match attr.key {
                        b"name" => {
                            param.name = attr.unescape_and_decode_value(reader).expect("Error decoding name");

                        },
                        b"value" => {
                            param.value = attr.unescape_and_decode_value(reader).expect("Error decoding value");
                        },
                        b"cvRef" => {
                            param.controlled_vocabulary = Some(attr.unescape_and_decode_value(reader).expect("Error decoding CV Ref"));
                        },
                        b"accession" => {
                            param.accession = attr.unescape_and_decode_value(reader).expect("Error decoding accession");
                        },
                        b"unitName" => {
                            param.unit_info = Some(attr.unescape_and_decode_value(reader).expect("Error decoding unit name"));
                        },
                        b"unitAccession" => {

                        },
                        b"unitCvRef" => {

                        },
                        _ => {}
                    }
                },
                Err(msg) => {
                    panic!("{} during parameter parsing.", self.handle_xml_error(msg));
                }
            }
        }
        return param;
    }

    pub fn handle_xml_error(&self, _error: quick_xml::Error) -> String {
        return String::from("Unknown Error");
    }

    pub fn fill_isolation_window(&mut self, param: Param) {
        let window = &mut self.precursor.isolation_window;
        match param.name.as_ref() {
            "isolation window target m/z" => {
                window.target = param.coerce().expect("Failed to parse isolation window target");
                window.flags = match window.flags {
                    IsolationWindowState::Unknown => {
                        IsolationWindowState::Complete
                    },
                    IsolationWindowState::Explicit => {
                        IsolationWindowState::Complete
                    },
                    IsolationWindowState::Offset => {
                        window.lower_bound = window.target - window.lower_bound;
                        window.upper_bound = window.target + window.upper_bound;
                        IsolationWindowState::Complete
                    },
                    IsolationWindowState::Complete => {
                        IsolationWindowState::Complete
                    }
                };
            },
            "isolation window lower offset" => {
                let lower_bound: f64 = param.coerce().expect("Failed to parse isolation window limit");
                match window.flags {
                    IsolationWindowState::Unknown => {
                        window.flags = IsolationWindowState::Offset;
                        window.lower_bound = lower_bound;
                    },
                    IsolationWindowState::Complete => {
                        window.lower_bound = window.target - lower_bound;
                    },
                    _ => {}
                }
            },
            "isolation window upper offset" => {
                let upper_bound: f64 = param.coerce().expect("Failed to parse isolation window limit");
                match window.flags {
                    IsolationWindowState::Unknown => {
                        window.flags = IsolationWindowState::Offset;
                        window.upper_bound = upper_bound;
                    },
                    IsolationWindowState::Complete => {
                        window.upper_bound = window.target + upper_bound;
                    },
                    _ => {}
                }
            },
            "isolation window lower limit" => {
                let lower_bound: f64 = param.coerce().expect("Failed to parse isolation window limit");
                match window.flags {
                    IsolationWindowState::Unknown => {
                        window.flags = IsolationWindowState::Explicit;
                        window.lower_bound = lower_bound;
                    },
                    _ => {}
                }
            },
            "isolation window upper limit" => {
                let upper_bound: f64 = param.coerce().expect("Failed to parse isolation window limit");
                match window.flags {
                    IsolationWindowState::Unknown => {
                        window.flags = IsolationWindowState::Explicit;
                        window.upper_bound = upper_bound;
                    },
                    _ => {}
                }
            },
            &_ => {

            }
        }
    }

    pub fn fill_scan_window(&mut self, param: Param) {
        let window: &mut ScanWindow = self.acquisition.scans.last_mut().unwrap().scan_windows.last_mut().unwrap();
        match param.name.as_ref() {
            "scan window lower limit" => {
                window.lower_bound = param.coerce().expect("Failed to parse scan window limit");
            },
            "scan window upper limit" => {
                window.upper_bound = param.coerce().expect("Failed to parse scan window limit");
            },
            &_ => {

            }
        }
    }

    pub fn fill_selected_ion(&mut self, param: Param) {
        match param.name.as_ref() {
            "selected ion m/z" => {
                self.precursor.ion.mz = param.coerce().expect("Failed to parse ion m/z");
            },
            "peak intensity" => {
                self.precursor.ion.intensity = param.coerce().expect("Failed to parse peak intensity");
            },
            "charge state" => {
                self.precursor.ion.charge = param.coerce().expect("Failed to parse ion charge");
            },
            &_ => {
                self.precursor.ion.params.push(param);
            }
        };
    }

    pub fn fill_binary_data_array(&mut self, param: Param) {
        match param.name.as_ref() {
            // Compression types
            "zlib compression" => {
                self.current_array.compression = BinaryCompressionType::Zlib;
            },
            "no compression" => {
                self.current_array.compression = BinaryCompressionType::NoCompression;
            },

            // Array data types
            "64-bit float" => {
                self.current_array.dtype = BinaryDataArrayType::Float64;
            },
            "32-bit float" => {
                self.current_array.dtype = BinaryDataArrayType::Float32;
            },
            "64-bit integer" => {
                self.current_array.dtype = BinaryDataArrayType::Int64;
            },
            "32-bit integer" => {
                self.current_array.dtype = BinaryDataArrayType::Int32;
            },
            "null-terminated ASCII string" => {
                self.current_array.dtype = BinaryDataArrayType::ASCII;
            },

            // Array types
            "m/z array" => {
                self.current_array.name = ArrayType::MZArray
            },
            "intensity array" => {
                self.current_array.name = ArrayType::IntensityArray
            },
            "charge array" => {
                self.current_array.name = ArrayType::ChargeArray
            },
            "non-standard data array" => {
                self.current_array.name = ArrayType::NonStandardDataArray {name: param.value};
            },
            "mean ion mobility array" | "mean drift time array" | "mean inverse reduced ion mobility array" => {
                self.current_array.name = ArrayType::MeanIonMobilityArray
            },
            "ion mobility array" | "drift time array" | "inverse reduced ion mobility array" => {
                self.current_array.name = ArrayType::IonMobilityArray
            },
            "deconvoluted ion mobility array" | "deconvoluted drift time array" | "deconvoluted inverse reduced ion mobility array" => {
                self.current_array.name = ArrayType::DeconvolutedIonMobilityArray
            },

            &_ => {
                self.current_array.params.push(param);
            }
        }
    }

    pub fn start_element<B: io::BufRead>(&mut self, event: &BytesStart, state: MzMLParserState, reader: &Reader<B>) -> ParserResult {
        let elt_name = event.name();
        match elt_name {
            b"spectrum" => {
                for attr_parsed in event.attributes() {
                    match attr_parsed {
                        Ok(attr) => {
                            match attr.key {
                                b"id" => {
                                    self.scan_id = attr.unescape_and_decode_value(reader).expect("Error decoding id");
                                },
                                b"index" => {
                                    self.index = (&String::from_utf8_lossy(&attr.value)).parse::<usize>().expect("Failed to parse index");
                                },
                                _ => {}
                            }
                        },
                        Err(msg) => {
                            return Err(self.handle_xml_error(msg));
                        }
                    }
                }

                return Ok(MzMLParserState::Spectrum);

            }
            b"spectrumList" => {
                return Ok(MzMLParserState::SpectrumList);
            }
            b"cvParam" | b"userParam" => {
                let param = self.handle_param(event, reader);
                match state {
                    MzMLParserState::Spectrum => {
                        self.params.push(param);
                    },
                    MzMLParserState::ScanList => {
                        self.acquisition.params.push(param)
                    }
                    MzMLParserState::Scan => {
                        self.acquisition.scans.last_mut().unwrap().params.push(param)
                    },
                    MzMLParserState::ScanWindowList => {
                        self.acquisition.scans.last_mut().unwrap().params.push(param)
                    },
                    MzMLParserState::ScanWindow => {
                        self.fill_scan_window(param);
                    },
                    MzMLParserState::IsolationWindow => {
                        self.fill_isolation_window(param);
                    },
                    MzMLParserState::SelectedIon | MzMLParserState::SelectedIonList => {
                        self.fill_selected_ion(param);
                    },
                    MzMLParserState::Activation => {
                        match param.name.as_ref() {
                            "collision energy" | "activation energy" => {
                                self.precursor.activation.energy = param.coerce().expect("Failed to parse collision energy");
                            },
                            &_ => {
                                self.precursor.activation.params.push(param);
                            }
                        }
                    },
                    MzMLParserState::BinaryDataArrayList => {},
                    MzMLParserState::BinaryDataArray => {
                        self.fill_binary_data_array(param);
                    },
                    MzMLParserState::Precursor | MzMLParserState::PrecursorList => {
                        self.precursor.params.push(param);
                    },
                    _ => {}

                };
                return Ok(state);
            },
            b"scanList" => {
                return Ok(MzMLParserState::ScanList);
            },
            b"scan" => {
                return Ok(MzMLParserState::Scan);
            },
            b"scanWindow" => {
                return Ok(MzMLParserState::ScanWindow);
            },
            b"scanWindowList" => {
                return Ok(MzMLParserState::ScanWindowList);
            },
            b"precursorList" => {
                return Ok(MzMLParserState::PrecursorList);
            },
            b"precursor" => {
                return Ok(MzMLParserState::Precursor);
            },
            b"isolationWindow" => {
                return Ok(MzMLParserState::IsolationWindow);
            },
            b"selectedIonList" => {
                return Ok(MzMLParserState::SelectedIonList);
            },
            b"selectedIon" => {
                return Ok(MzMLParserState::SelectedIon);
            },
            b"activation" => {
                return Ok(MzMLParserState::Activation);
            },
            b"binaryDataArrayList" => {
                return Ok(MzMLParserState::BinaryDataArrayList);
            },
            b"binaryDataArray" => {
                return Ok(MzMLParserState::BinaryDataArray);
            },
            b"binary" => {
                return Ok(MzMLParserState::Binary);
            }
            _ => {}
        };
        return Ok(state);
    }

    pub fn end_element(&mut self, event: &BytesEnd, state: MzMLParserState) -> ParserResult {
        let elt_name = event.name();
        match elt_name {
            b"spectrum" => {
                return Ok(MzMLParserState::SpectrumDone)
            }
            b"scanList" => {
                return Ok(MzMLParserState::Spectrum)
            },
            b"scan" => {
                return Ok(MzMLParserState::ScanList)
            },
            b"scanWindow" => {
                return Ok(MzMLParserState::ScanWindowList)
            },
            b"scanWindowList" => {
                return Ok(MzMLParserState::Scan)
            },
            b"precursorList" => {
                return Ok(MzMLParserState::Spectrum)
            },
            b"precursor" => {
                return Ok(MzMLParserState::PrecursorList)
            },
            b"isolationWindow" => {
                return Ok(MzMLParserState::Precursor)
            },
            b"selectedIonList" => {
                return Ok(MzMLParserState::Precursor)
            },
            b"selectedIon" => {
                return Ok(MzMLParserState::SelectedIonList)
            },
            b"activation" => {
                return Ok(MzMLParserState::Precursor)
            },
            b"binaryDataArrayList" => {
                return Ok(MzMLParserState::Spectrum);
            },
            b"binaryDataArray" => {
                let array = self.current_array.clone();
                self.arrays.add(array);
                return Ok(MzMLParserState::BinaryDataArrayList)
            },
            b"binary" => {
                return Ok(MzMLParserState::BinaryDataArray)
            }
            _ => {

            }
        };
        return Ok(state);
    }

    pub fn text(&mut self, event: &BytesText, state: MzMLParserState) -> ParserResult {
        match state {
            MzMLParserState::Binary => {
                let bin = event.unescaped().expect("Failed to unescape binary data array content");
                self.current_array.data = Bytes::from(&*bin);
            },
            _ => {}
        }
        return Ok(state);
    }
}


pub struct MzMLReader<R: io::Read> {
    pub state: MzMLParserState,
    pub handle: io::BufReader<R>,
    pub error: String,
}


impl<R: io::Read> MzMLReader<R> {
    pub fn new(file: R) -> MzMLReader<R> {
        let handle = io::BufReader::with_capacity(BUFFER_SIZE, file);
        return MzMLReader {
            handle: handle,
            state: MzMLParserState::Start,
            error: String::from(""),
        }
    }

    fn _parse_into(mut self, buffer: &mut Bytes, accumulator: &mut MzMLSpectrumBuilder) -> Self {
        let mut reader = Reader::from_reader(self.handle);

        reader.trim_text(true);
        loop {
            match reader.read_event(buffer) {
                Ok(Event::Start(ref e)) => {
                    match accumulator.start_element(e, self.state, &reader) {
                        Ok(state) => {self.state = state;},
                        Err(message) => {
                            self.state = MzMLParserState::ParserError;
                            self.error = message;
                        }
                    }
                },
                Ok(Event::End(ref e)) => {
                    match accumulator.end_element(e, self.state) {
                        Ok(state) => {self.state = state;},
                        Err(message) => {
                            self.state = MzMLParserState::ParserError;
                            self.error = message;
                        }
                    }
                },
                Ok(Event::Text(ref e)) => {
                    match accumulator.text(e, self.state) {
                        Ok(state) => {self.state = state;},
                        Err(message) => {
                            self.state = MzMLParserState::ParserError;
                            self.error = message;
                        }
                    }
                },
                Ok(Event::Eof) => {
                    break;
                },
                Err(err) => panic!("Error at position {}: {:?}", reader.buffer_position(), err),
                _ => {}
            }
            buffer.clear();
        }
        self.handle = reader.into_underlying_reader();
        return self;
    }

    pub fn get(mut self) {
        let mut buffer = Bytes::new();
        let mut accumulator = MzMLSpectrumBuilder::new();
        self = self._parse_into(&mut buffer, &mut accumulator);
        println!("{}", self.error);
    }


}
