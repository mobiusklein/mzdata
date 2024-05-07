set dotenv-load := true

test-units:
    cargo test --lib --features nalgebra,parallelism,mzsignal,mzmlb,zlib-ng-compat,thermo

alias t := test-units

test-units-more:
    cargo test --lib --features nalgebra,parallelism,mzsignal,mzmlb,zlib-ng-compat,thermo,async

docs:
    cargo doc --no-deps --features nalgebra,parallelism,mzsignal,mzmlb,zlib-ng-compat,thermo,async -p mzdata

install-mzdata:
    cargo install --path . --features nalgebra,parallelism,mzsignal,mzmlb,zlib-ng-compat,hdf5_static

update-cv:
    curl --insecure \
     --location \
     https://github.com/HUPO-PSI/psi-ms-CV/releases/latest/download/psi-ms.obo | gzip -c > cv/psi-ms.obo.gz

update-cv-terms:
    cog -c -r -U src/meta/software.rs src/meta/instrument.rs

changelog version:
    #!/usr/bin/env python

    import subprocess
    import re

    new_content = subprocess.check_output(['git', 'cliff', '-s', 'all', '-u', '-t', '{{version}}'], stderr=subprocess.DEVNULL).decode()

    new_version = "{{version}}"

    buffer = open('CHANGELOG.md').read()

    buffer = buffer.replace("## ", f"{new_content}## ", 1).splitlines()

    offset = buffer.index("<!-- Versions -->") + 1
    line_to_patch = buffer[offset]
    buffer[offset] = re.sub(r"v\d+\.\d+\.\d+[^\.]*", new_version, line_to_patch)
    version_link_template = buffer[offset + 1]
    version_link_template = re.sub(
        r"\d+\.\d+\.\d+[^\.]*(?=\])", new_version[1:], version_link_template
    )
    version_link_template = re.sub(
        r"v\d+\.\d+\.\d+[^\.]*", new_version, version_link_template
    )

    buffer.insert(offset + 1, version_link_template)

    buffer = '\n'.join(buffer)
    open('CHANGELOG.md', 'wt').write(buffer)

