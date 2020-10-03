//! A library and command line tool for spatial color quantization.
//!
//! # Overview
//!
//! Rust port of Derrick Coetzee's [`scolorq`][scolorq], based on the 1998 paper
//! "On spatial quantization of color images" by Jan Puzicha, Markus Held, Jens
//! Ketterer, Joachim M. Buhmann, & Dieter Fellner. *Spatial quantization* is
//! defined as simultaneously performing halftoning (dithering) and color
//! quantization (limiting the colors in an image). For more information, visit
//! [the original implementation's website][scolorq].
//!
//! The algorithm is excellent for retaining image detail and minimizing visual
//! distortions for color palettes in the neighborhood of 4, 8, or 16 colors,
//! especially as the image size is reduced. It combines limiting the color palette
//! and dithering the image into a simultaneous process as opposed to sequentially
//! limiting the colors then dithering. Colors are chosen based on their context in
//! the image, hence the "spatial" aspect of spatial color quantization. The
//! colors are selected based on their neighbors to mix as an average illusory
//! color in the human eye.
//!
//! To use as a library, add the following to your `Cargo.toml`; add the
//! `palette_color` feature to enable Lab color quantization. See the
//! [README.md][readme] for image examples, output, and usage.
//!
//! ```toml
//! [dependencies.rscolorq]
//! version = "0.1"
//! default-features = false
//! ```
//!
//! [readme]: https://github.com/okaneco/rscolorq/blob/master/README.md
//! [scolorq]: http://people.eecs.berkeley.edu/~dcoetzee/downloads/scolorq/
//!
//! ## Usage
//!
//! The following example shows the mapping of an image buffer in Rgb from
//! `[u8; 3]` to `[f64; 3]`, performing the color quantization, then filling
//! a buffer with `u8` to be saved as an image.
//!
//! [Matrix2d]: struct.Matrix2d.html
//! ```
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! # let width = 2;
//! # let height = 1;
//! # let palette_size = 4;
//! # let img = vec![[0, 0, 0], [255, 255, 255]];
//! use rscolorq::{color::Rgb, spatial_color_quant, Matrix2d, Params};
//!
//! // Create the output buffer and quantized palette index buffer
//! let mut imgbuf = Vec::with_capacity(width * height * 3);
//! let mut quantized_image = Matrix2d::new(width, height);
//!
//! // Build the quantization parameters, verify if accepting user input
//! let mut conditions = Params::new();
//! conditions.palette_size(palette_size);
//! conditions.verify_parameters()?;
//!
//! // Convert the input image buffer from Rgb<u8> to Rgb<f64>
//! let image = Matrix2d::from_vec(
//!     img.iter()
//!         .map(|&c| Rgb {
//!             red: c[0] as f64 / 255.0,
//!             green: c[1] as f64 / 255.0,
//!             blue: c[2] as f64 / 255.0,
//!         })
//!         .collect(),
//!     width,
//!     height,
//! );
//!
//! let mut palette = Vec::with_capacity(palette_size as usize);
//!
//! spatial_color_quant(&image, &mut quantized_image, &mut palette, &conditions)?;
//!
//! // Convert the Rgb<f64> palette to Rgb<u8>
//! let palette = palette
//!     .iter()
//!     .map(|&c| {
//!         let color = 255.0 * c;
//!         [
//!             color.red.round() as u8,
//!             color.green.round() as u8,
//!             color.blue.round() as u8,
//!         ]
//!     })
//!     .collect::<Vec<[u8; 3]>>();
//!
//! // Create the final image by color lookup from the palette
//! quantized_image.iter().for_each(|&c| {
//!     imgbuf.extend_from_slice(&*palette.get(c as usize).unwrap());
//! });
//!
//! # Ok(())
//! # }
//! ```
//!
//! ## Features
//! - use RGB or Lab color space for calculations
//! - can dither based on fixed color palette
//! - seedable RNG for reproducible results
//!
//! ## Limitations
//!
//! ### It's "slow"
//! - Larger images or images with smooth transitions/gradients will take longer.
//! Higher palette sizes will take longer.
//! - The algorithm is suited towards retaining detail with smaller color palettes.
//! You can still use it on larger images but be aware it's not close to real-time
//! unless the image is small.
//!
//! ### Filter size 1x1
//! - Doesn't produce an image resembling the input, nor does the original.
//!
//! ### Filter size 5x5
//! - Doesn't always converge.
//! - I'm unsure if this is an error in this implementation or a problem with the
//! random number generator being used. The original implementation may take a while
//! but eventually completes with filter size 5.
//! - Any help on this would be appreciated.
#![warn(missing_docs, rust_2018_idioms)]

pub mod color;
mod error;
mod matrix;
mod quant;

pub use error::QuantError;
pub use matrix::{Matrix2d, Matrix3d};
pub use quant::{spatial_color_quant, FilterSize, Params};

/// A trait required to calculate the spatial quantization on a color type.
pub trait SpatialQuant: Sized {
    /// Calculates the filter weights matrix based on the
    /// [`FilterSize`](enum.FilterSize.html).
    fn calculate_filter_weights(dithering_level: f64, filter_size: FilterSize) -> Matrix2d<Self>;
    /// Calculates the difference between two colors.
    fn color_difference(&self, other: &Self) -> f64;
    /// Returns the minimal threshold before a color is considered to be
    /// different.
    fn difference_threshold() -> f64;
    /// Multiplies the components of a color directly with the corresponding
    /// components of `other`.
    fn direct_product(&self, other: &Self) -> Self;
    /// Calculates the dot product.
    fn dot_product(&self, other: &Self) -> f64;
    /// Calculate the squared magnitude of `self` and `other`.
    fn norm_squared(&self) -> f64;
    /// Create a random color.
    fn random(rng: &mut impl rand::Rng) -> Self;
    /// Update the color palette.
    fn refine_palette(
        s: &mut Matrix2d<Self>,
        coarse_variables: &Matrix3d<f64>,
        a: &Matrix2d<Self>,
        palette: &mut Vec<Self>,
    ) -> Result<(), QuantError>;
    /// Either one or the maximal intensity of the color type.
    ///
    /// Note: This is included for thoroughness. It's only used for the `One`
    /// filter size, and the `One` filter does not produce a usable result.
    fn one() -> Self;
}

/// A trait for calculating the inverse of a matrix.
pub trait MatrixComponent {
    /// Returns the identity of the component.
    fn identity() -> Self;
    /// Returns the reciprocal of the component.
    fn inverse(&self) -> Self;
}

impl MatrixComponent for f32 {
    fn identity() -> Self {
        1.0
    }

    fn inverse(&self) -> Self {
        self.recip()
    }
}

impl MatrixComponent for f64 {
    fn identity() -> Self {
        1.0
    }

    fn inverse(&self) -> Self {
        self.recip()
    }
}
