import * as wasm from "mzdata-wasm";

export class IMMZReader {
  reader: wasm.WebIMMZReader;

  static async fromByteArray(bytes: Uint8Array) {
    return new IMMZReader(await wasm.WebIMMZReader.byte_array(bytes));
  }

  static async fromBuffer(buffer: ArrayBuffer) {
    return new IMMZReader(await wasm.WebIMMZReader.buffer(buffer));
  }

  static async fromBlob(blob: Blob) {
    return new IMMZReader(await wasm.WebIMMZReader.blob(blob));
  }

  static async fromFile(file: File) {
    return new IMMZReader(await wasm.WebIMMZReader.file(file));
  }

  static __wrap(reader: wasm.WebIMMZReader)  {
    return new IMMZReader(reader)
  }

  private constructor(reader: wasm.WebIMMZReader) {
    this.reader = reader;
  }
}

export class MZReader {
  reader: wasm.WebMZReader;

  static async fromByteArray(bytes: Uint8Array) {
    return new MZReader(await wasm.WebMZReader.byte_array(bytes));
  }

  static async fromBuffer(buffer: ArrayBuffer) {
    return new MZReader(await wasm.WebMZReader.buffer(buffer));
  }

  static async fromBlob(blob: Blob) {
    return new MZReader(await wasm.WebMZReader.blob(blob));
  }

  static async fromFile(file: File) {
    return new MZReader(await wasm.WebMZReader.file(file));
  }

  private constructor(reader: wasm.WebMZReader) {
    this.reader = reader;
  }

  fileFormat() {
    this.reader.file_format;
  }

  setDataLoading(value: boolean) {
    this.reader.set_data_loading(value);
    return this;
  }

  setPeakPicking(value: boolean) {
    this.reader.set_peak_picking(value);
    return this;
  }

  get length() {
    return this.reader.length;
  }

  at(index: number) {
    return this.getSpectrumByIndex(index);
  }

  *[Symbol.iterator]() {
    const n = this.reader.length;
    for (let i = 0; i < n; i++) {
      yield this.getSpectrumByIndex(i);
    }
  }

  iter() {
    return this[Symbol.iterator]();
  }

  getSpectrumByIndex(index: number) {
    return this.reader.get_spectrum_by_index(index);
  }

  getSpectrumById(id: string) {
    return this.reader.get_spectrum_by_id(id);
  }

  getSpectrumByTime(time: number) {
    return this.reader.get_spectrum_by_time(time);
  }

  hasIonMobility() {
    return this.reader.has_ion_mobility_dimension()
  }

  asIonMobilityReader() {
    return this.reader.as_ion_mobility_reader().then((reader) => {
      reader === undefined ? null : IMMZReader.__wrap(reader);
    })
  }
}