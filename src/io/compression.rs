use std::path;


pub fn is_gzipped(header: &[u8]) -> bool {
    header.starts_with(b"\x1f\x8b")
}


pub fn is_gzipped_extension(path: path::PathBuf) -> (bool, path::PathBuf) {
    if let Some(ext) = path.extension() {
        if ext.to_ascii_lowercase() == "gz" {
            (true, path.with_extension(""))
        } else {
            (false, path)
        }
    } else {
        (false, path)
    }
}

