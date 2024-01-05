test-units:
    cargo test --lib --features nalgebra,parallelism,mzsignal,mzmlb,zlib-ng-compat

docs:
    cargo doc --no-deps --features nalgebra,parallelism,mzsignal,mzmlb,zlib-ng-compat
