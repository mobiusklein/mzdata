export RUST_LOG=info

RAW=$1
OUTFILE=$2

cargo b --example mzconvert --features mzmlb

msconvert --mzML $RAW -o- | cargo r --example mzconvert --features mzmlb - $OUTFILE 2> log.txt