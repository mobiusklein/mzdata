set dotenv-load := true

test-units:
    cargo nextest run --lib --features nalgebra,parallelism,mzsignal,zlib-ng-compat,thermo,async,numpress

test-coverage:
    cargo llvm-cov --lib --tests nextest --features nalgebra,parallelism,mzsignal,zlib-ng-compat,thermo,mzmlb,async,numpress --html

alias t := test-units

test-units-more:
    cargo nextest run --lib --features nalgebra,parallelism,mzsignal,zlib-ng-compat,thermo,async,numpress

quick-docs:
    cargo doc --no-deps -p mzdata

docs:
    cargo doc --no-deps --features nalgebra,parallelism,mzsignal,mzmlb,zlib-ng-compat,thermo,async,proxi,bruker_tdf -p mzdata -p mzsignal -p mzpeaks

install-mzdata:
    cargo install --path . --features nalgebra,parallelism,mzsignal,mzmlb,zlib-ng-compat,hdf5_static

update-cv:
    curl --insecure \
     --location \
     https://github.com/HUPO-PSI/psi-ms-CV/releases/latest/download/psi-ms.obo | gzip -c > cv/psi-ms.obo.gz

    gzip -d -c cv/psi-ms.obo.gz | head -n 5

update-cv-terms:
    cog -c -r -U src/meta/software.rs src/meta/instrument.rs src/meta/file_description.rs src/io/mzml/writer.rs src/meta/activation.rs

changelog version:
    #!/usr/bin/env python

    import subprocess
    import re

    new_content = subprocess.check_output(['git', 'cliff', '-s', 'all', '-u', '-t', '{{version}}'], stderr=subprocess.DEVNULL).decode()

    new_version = "{{version}}"

    buffer = open('CHANGELOG.md').read()

    buffer = buffer.replace("## ", f"{new_content}## ", 1).splitlines()

    offset = buffer.index("<!-- Versions -->") + 1
    line_to_patch = buffer[offset + 1]
    previous_version = re.search(r"(v\d+\.\d+\.\d+[^\.]*)", line_to_patch).group(1)
    buffer[offset] = re.sub(r"v\d+\.\d+\.\d+[^\.]*", new_version, line_to_patch)

    version_link_template = buffer[offset + 2]
    version_link_template = re.sub(
        r"\d+\.\d+\.\d+[^\.]*(?=\])", new_version[1:], version_link_template
    )
    version_link_template = version_link_template.rsplit("/", 1)[0] + f"/{previous_version}...{new_version}"
    buffer[offset + 1] = version_link_template

    buffer.insert(offset, '')
    buffer = '\n'.join(buffer)
    open('CHANGELOG.md', 'wt').write(buffer)
    print(buffer)


release tag: (patch-version) (changelog tag)
    #!/usr/bin/env bash

    git add crates
    git add CHANGELOG.md
    git commit -m "chore: update changelog"
    git tag {{tag}}

    cargo publish
    cd crates/mzdata-spectra && cargo publish


patch-version:
    #!/usr/bin/env python
    import sys
    import re


    ref_toml = "Cargo.toml"
    target_toml = "crates/mzdata-spectra/Cargo.toml"

    pattern = re.compile(r"^version\s*=\s*\"(.+?)\"")
    dep_pattern = re.compile(r"version\s*=\s*\"(.+?)\"")

    version = None

    with open(ref_toml) as fh:
        for line in fh:
            if match := pattern.match(line):
                version = match.string
                break

    if not version:
        raise ValueError("Version not found in reference")


    buffer = []
    with open(target_toml) as fh:
        for line in fh:
            if pattern.match(line):
                line = version
            if line.startswith("mzdata"):
                line = dep_pattern.sub(version.strip(), line)

            buffer.append(line.strip())

    with open(target_toml, 'w') as fh:
        fh.write('\n'.join(buffer))
