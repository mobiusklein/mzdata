"""File-level metadata: file description, instrument configs, software, samples, processing."""

from __future__ import annotations


def test_file_description(reader) -> None:
    description = reader.file_description()
    assert isinstance(description.source_files(), list)
    assert isinstance(description.has_msn_spectra, bool)
    for source in description.source_files():
        assert isinstance(source.name, str)
        assert isinstance(source.location, str)


def test_instrument_configurations(reader) -> None:
    configs = reader.instrument_configurations()
    assert isinstance(configs, list)
    for config in configs:
        assert isinstance(config.id, int)
        for component in config.components():
            assert component.order >= 0
            for accessor in (component.mass_analyzer, component.detector, component.ionization_type):
                assert accessor is None or isinstance(accessor, str)


def test_softwares(reader) -> None:
    softwares = reader.softwares()
    assert isinstance(softwares, list)
    for software in softwares:
        assert isinstance(software.id, str)
        assert isinstance(software.version, str)


def test_samples(reader) -> None:
    for sample in reader.samples():
        assert isinstance(sample.id, str)
        assert sample.name is None or isinstance(sample.name, str)


def test_run_description(reader) -> None:
    run = reader.run_description()
    if run is not None:
        assert run.id is None or isinstance(run.id, str)
        assert run.default_instrument_id is None or isinstance(run.default_instrument_id, int)


def test_data_processings(reader) -> None:
    for processing in reader.data_processings():
        assert isinstance(processing.id, str)
        for method in processing.methods():
            assert isinstance(method.order, int)
            assert isinstance(method.software_reference, str)
