import { write_spectra_to_mgf, Spectrum, _write_spectra_to_mzml } from "mzdata-wasm"
import { MZReader } from "./mem_reader"


export function writeMGF(spectra: Spectrum[]) {
    return write_spectra_to_mgf(spectra.map((s) => s.copy()))
}


export function writeMzML(spectra: Spectrum[], reader: MZReader) {
    return _write_spectra_to_mzml(
        spectra.map((s) => s.copy()),
        reader.reader
    )
}
