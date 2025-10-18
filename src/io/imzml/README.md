# imzML Support in mzdata

This module provides support for reading imzML (Imaging Mass Spectrometry Markup Language) files in the mzdata library.

## Overview

imzML is a data format for mass spectrometry imaging (MSI) data that extends the mzML format. Instead of storing binary data as base64-encoded text within the XML, imzML stores the binary data in a separate `.ibd` (imaging binary data) file. This allows for more efficient storage and access of large imaging datasets.

## File Structure

An imzML dataset consists of two files:
- `.imzML` file: XML metadata based on the mzML schema with imaging-specific controlled vocabulary terms
- `.ibd` file: Binary data file containing the actual spectral data

## Data Modes

imzML supports two data storage modes:

1. **Continuous Mode**: All spectra share the same m/z values, stored once at the beginning of the `.ibd` file
2. **Processed Mode**: Each spectrum has its own m/z and intensity arrays

## Usage

```rust
use mzdata::io::ImzMLReader;

// Open an imzML file (will automatically look for the .ibd file)
let mut reader = ImzMLReader::open_path("sample.imzML")?;

// Iterate through spectra
for spectrum in reader {
    println!("Spectrum ID: {}", spectrum.id());
    println!("Number of peaks: {}", spectrum.peaks().len());
}

// Or access spectra by ID/index
let spectrum = reader.get_spectrum_by_id("spectrum=1")?;
```

## Implementation Details

- The `ImzMLReader` wraps the existing `MzMLReader` and extends it with `.ibd` file support
- Format detection is automatic based on the presence of imaging-specific CV terms
- Currently provides basic functionality with room for enhancement (e.g., actual IBD data reading)

## Controlled Vocabulary

imzML uses additional controlled vocabulary terms defined in `imagingMS.obo` for imaging-specific metadata such as:
- `IMS:1000101`: external data
- `IMS:1000102`: external offset  
- `IMS:1000103`: external array length
- `IMS:1000104`: external encoded length

## References

- [imzML Specification](https://www.ms-imaging.org/imzml/)
- [Imaging MS Ontology](https://www.ms-imaging.org/controlled-vocabulary/)
- [Data Structure Documentation](https://www.ms-imaging.org/imzml/data-structure/)
