set dotenv-load := true

test-units:
    cargo nextest run --workspace --lib --features nalgebra,parallelism,mzsignal,zlib-ng-compat,thermo,async,numpress,bruker_tdf,imzml,cv

test-coverage:
    cargo llvm-cov --lib --workspace --tests nextest --features nalgebra,parallelism,mzsignal,zlib-ng-compat,thermo,async,numpress,imzml,cv --html

alias t := test-units

test-units-more:
    cargo nextest run --lib --features nalgebra,parallelism,mzsignal,zlib-ng-compat,thermo,async,numpress,bruker_tdf,imzml,cv

quick-docs:
    cargo doc --no-deps -p mzdata

docs:
    cargo doc --no-deps --features nalgebra,parallelism,mzsignal,zlib-ng-compat,thermo,async,proxi,bruker_tdf,imzml,cv -p mzdata -p mzsignal -p mzpeaks

install-mzdata:
    cargo install --path . --features nalgebra,parallelism,mzsignal,mzmlb,zlib-ng-compat,hdf5_static

update-cv:
    curl --insecure \
     --location \
     https://github.com/HUPO-PSI/psi-ms-CV/releases/latest/download/psi-ms.obo | gzip -c > cv/psi-ms.obo.gz

    gzip -d -c cv/psi-ms.obo.gz | head -n 5
    cargo r -p cvbuild cv/psi-ms.obo.gz "crates/mzdata-param/src/ms.dat"

update-cv-build:
    cargo r -p cvbuild cv/psi-ms.obo.gz "crates/mzdata-param/src/ms.dat"

update-cv-terms:
    cog -c -r -U \
        src/io/mzml/writer.rs \
        crates/mzdata-meta/src/software.rs \
        crates/mzdata-meta/src/instrument.rs \
        crates/mzdata-meta/src/file_description.rs \
        crates/mzdata-meta/src/activation.rs

pytest:
    maturin develop -m "./crates/pymzdata/Cargo.toml"
    pytest -v -s ./crates/pymzdata/

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
    set -e

    git add crates
    git add CHANGELOG.md
    git commit -m "chore: update changelog"
    git tag {{tag}}

    cargo publish --allow-dirty -p mzdata-param -p mzdata-bindata -p mzdata-meta -p mzdata-spectrum  -p mzdata
    cd crates/mzdata-spectra && cargo publish --allow-dirty


patch-version tag="-":
    #!/usr/bin/env python
    import sys
    import re


    ref_toml = "Cargo.toml"
    pattern = re.compile(r"^version\s*=\s*\"(.+?)\"")
    dep_pattern = re.compile(r"version\s*=\s*\"(.+?)\"")

    version = None
    if "{{tag}}" == "-":
        with open(ref_toml) as fh:
            for line in fh:
                if match := pattern.match(line):
                    version = match.string
                    break
    else:
        tag = "{{tag}}".replace("v", "")
        version = f"version = \"{tag}\""
        with open(ref_toml) as fh:
            lines = fh.readlines()
        buf = []
        for line in lines:
            if match := pattern.match(line):
                line = version + '\n'
            buf.append(line)
        with open(ref_toml, 'w') as fh:
            fh.writelines(buf)

    if version is None:
        sys.stderr.write("Failed to resolve reference mzdata version")
        if not version:
            raise ValueError("Version not found in reference")
    else:
        sys.stderr.write(f"mzdata version = {version}")

    with open(ref_toml) as fh:
        lines = fh.readlines()
    in_workspace = False
    buf = []
    for line in lines:
        if not in_workspace and line.startswith("[workspace.dependencies"):
            in_workspace = True
        elif in_workspace and line.startswith("mzdata-"):
            line_before = line
            line = dep_pattern.sub(version.strip(), line)
            mat_before = dep_pattern.search(line_before)
            mat_after = dep_pattern.search(line)
            print(line_before, "=>", line)
            print(mat_before, "=>", mat_after)
        buf.append(line)
    with open(ref_toml, 'w') as fh:
        fh.writelines(buf)


    target_toml_files = [
        "./crates/mzdata-spectra/Cargo.toml",
        "./crates/pymzdata/Cargo.toml",
        "./crates/mzdata-wasm/Cargo.toml"
    ]
    for target_toml in target_toml_files:
        print(f"Updating {target_toml}")
        buffer = []
        with open(target_toml) as fh:
            for line in fh:
                if pattern.match(line):
                    line = version
                if line.startswith("mzdata"):
                    line_before = line
                    line = dep_pattern.sub(version.strip(), line)
                    mat_before = dep_pattern.search(line_before)
                    mat_after = dep_pattern.search(line)
                    print(mat_before, "=>", mat_after)

                buffer.append(line.strip())

        with open(target_toml, 'w') as fh:
            fh.write('\n'.join(buffer))
