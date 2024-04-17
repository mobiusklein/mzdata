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