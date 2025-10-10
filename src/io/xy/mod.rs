mod reader;

pub use reader::*;

pub fn is_xy(buf: &[u8]) -> bool {
    if let Ok(string) = String::from_utf8(buf.to_vec()) {
        let trimmed = string.trim();
        match trimmed.split_once(' ').or(trimmed.split_once('\t')) {
            None => false,
            Some((mz, intensity)) => {
                mz.trim().parse::<f64>().is_ok() && intensity.trim().parse::<f32>().is_ok()
            }
        }
    } else {
        false
    }
}
