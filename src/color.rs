//! Colors that can be used for spatial color quantization.
//!
//! Enable the `palette_color` feature for Lab color space quantization.

#[cfg(feature = "palette_color")]
mod lab;
mod rgb;

pub use crate::color::rgb::Rgb;
