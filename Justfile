set dotenv-load := true

test-units:
    cargo test --lib --features nalgebra,parallelism,mzsignal,mzmlb,zlib-ng-compat

test-units-more:
    cargo test --lib --features nalgebra,parallelism,mzsignal,mzmlb,zlib-ng-compat,async

docs:
    cargo doc --no-deps --features nalgebra,parallelism,mzsignal,mzmlb,zlib-ng-compat,async -p mzdata

install-mzdata:
    cargo install --path . --features nalgebra,parallelism,mzsignal,mzmlb,zlib-ng-compat,hdf5_static