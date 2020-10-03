use crate::{quant::FilterSize, Matrix2d, SpatialQuant};

/// Red, green, and blue color with component range from 0.0 to 1.0.
#[derive(Clone, Copy, Debug, Default)]
pub struct Rgb {
    /// Red color channel.
    pub red: f64,
    /// Green color channel.
    pub green: f64,
    /// Blue color channel.
    pub blue: f64,
}

impl Rgb {
    /// Create a new `Rgb` struct.
    pub fn new(red: f64, green: f64, blue: f64) -> Self {
        Self { red, green, blue }
    }

    /// Create an `Rgb` struct from a slice.
    pub fn from_slice(slice: &[f64; 3]) -> Self {
        Self {
            red: slice[0],
            green: slice[1],
            blue: slice[2],
        }
    }
}

impl core::ops::Add for Rgb {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            red: self.red + other.red,
            green: self.green + other.green,
            blue: self.blue + other.blue,
        }
    }
}

impl core::ops::AddAssign for Rgb {
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            red: self.red + other.red,
            green: self.green + other.green,
            blue: self.blue + other.blue,
        };
    }
}

impl core::ops::DivAssign<f64> for Rgb {
    fn div_assign(&mut self, other: f64) {
        *self = Self {
            red: self.red / other,
            green: self.green / other,
            blue: self.blue / other,
        };
    }
}

impl core::ops::Mul for Rgb {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Self {
            red: self.red * other.red,
            green: self.green * other.green,
            blue: self.blue * other.blue,
        }
    }
}

impl core::ops::MulAssign for Rgb {
    fn mul_assign(&mut self, other: Self) {
        *self = Self {
            red: self.red * other.red,
            green: self.green * other.green,
            blue: self.blue * other.blue,
        };
    }
}

impl core::ops::Mul<f64> for Rgb {
    type Output = Self;

    fn mul(self, other: f64) -> Self {
        Self {
            red: self.red * other,
            green: self.green * other,
            blue: self.blue * other,
        }
    }
}

impl core::ops::MulAssign<f64> for Rgb {
    fn mul_assign(&mut self, other: f64) {
        *self = Self {
            red: self.red * other,
            green: self.green * other,
            blue: self.blue * other,
        };
    }
}

impl core::ops::Mul<Rgb> for f64 {
    type Output = Rgb;

    fn mul(self, other: Rgb) -> Self::Output {
        Rgb {
            red: self * other.red,
            green: self * other.green,
            blue: self * other.blue,
        }
    }
}

impl core::ops::Sub for Rgb {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            red: self.red - other.red,
            green: self.green - other.green,
            blue: self.blue - other.blue,
        }
    }
}

impl core::ops::SubAssign for Rgb {
    fn sub_assign(&mut self, other: Self) {
        *self = Self {
            red: self.red - other.red,
            green: self.green - other.green,
            blue: self.blue - other.blue,
        };
    }
}

impl SpatialQuant for Rgb {
    fn calculate_filter_weights(dithering_level: f64, filter_size: FilterSize) -> Matrix2d<Self> {
        let mut sum = 0.0;
        match filter_size {
            FilterSize::One => Matrix2d::from_vec(vec![Self::one()], 1, 1),
            FilterSize::Three => {
                let mut filter_weights = Matrix2d::new(3, 3);
                for i in 0i8..3i8 {
                    for j in 0i8..3i8 {
                        let val = (-f64::from((i - 1) * (i - 1) + (j - 1) * (j - 1)).sqrt()
                            / dithering_level.powi(2))
                        .exp();
                        sum += val;
                        *filter_weights.get_mut(i as usize, j as usize).unwrap() =
                            Self::new(val, val, val);
                    }
                }
                filter_weights.iter_mut().for_each(|fw| *fw /= sum);
                filter_weights
            }
            FilterSize::Five => {
                let mut filter_weights = Matrix2d::new(5, 5);
                for i in 0i8..5i8 {
                    for j in 0i8..5i8 {
                        let val = (-f64::from((i - 2) * (i - 2) + (j - 2) * (j - 2)).sqrt()
                            / dithering_level.powi(2))
                        .exp();
                        sum += val;
                        *filter_weights.get_mut(i as usize, j as usize).unwrap() =
                            Self::new(val, val, val);
                    }
                }
                filter_weights.iter_mut().for_each(|fw| *fw /= sum);
                filter_weights
            }
        }
    }

    fn color_difference(&self, other: &Self) -> f64 {
        (*self - *other).norm_squared()
    }

    fn difference_threshold() -> f64 {
        1.0 / (255.0 * 255.0)
    }

    fn direct_product(&self, other: &Self) -> Self {
        Self {
            red: self.red * other.red,
            green: self.green * other.green,
            blue: self.blue * other.blue,
        }
    }

    fn dot_product(&self, other: &Self) -> f64 {
        self.red * other.red + self.green * other.green + self.blue * other.blue
    }

    fn one() -> Self {
        Self {
            red: 1.0,
            green: 1.0,
            blue: 1.0,
        }
    }

    fn norm_squared(&self) -> f64 {
        self.red * self.red + self.green * self.green + self.blue * self.blue
    }

    fn random(rng: &mut impl rand::Rng) -> Self {
        Self {
            red: rng.gen_range(0.0, 1.0),
            green: rng.gen_range(0.0, 1.0),
            blue: rng.gen_range(0.0, 1.0),
        }
    }

    fn refine_palette(
        s: &mut crate::Matrix2d<Self>,
        coarse_variables: &crate::Matrix3d<f64>,
        a: &crate::Matrix2d<Self>,
        palette: &mut Vec<Self>,
    ) -> Result<(), crate::QuantError> {
        // Reflect s above the diagonal
        for v in 0..s.width() {
            for alpha in 0..v {
                *s.get_mut(v, alpha).unwrap() = *s.get(alpha, v).unwrap();
            }
        }

        let mut r: Vec<Self> = (0..palette.len()).map(|_| Self::default()).collect();
        for (v, r_item) in r.iter_mut().enumerate().take(palette.len()) {
            for i_y in 0..coarse_variables.height() {
                for i_x in 0..coarse_variables.width() {
                    *r_item +=
                        *coarse_variables.get(i_x, i_y, v).unwrap() * *a.get(i_x, i_y).unwrap();
                }
            }
        }

        // Iterate through each of the component channels, clamp to valid range
        // Red
        let s_k = Matrix2d::from_vec(s.iter().map(|a| a.red).collect(), s.width(), s.height());
        let r_k = r.iter().map(|a| a.red).collect::<Vec<_>>();
        let palette_channel = (((s_k * 2.0).matrix_inverse()?) * -1.0) * r_k;
        for (v, pal) in palette.iter_mut().zip(&palette_channel) {
            v.red = pal.min(1.0).max(0.0);
        }
        // Green
        let s_k = Matrix2d::from_vec(s.iter().map(|a| a.green).collect(), s.width(), s.height());
        let r_k = r.iter().map(|a| a.green).collect::<Vec<_>>();
        let palette_channel = (((s_k * 2.0).matrix_inverse()?) * -1.0) * r_k;
        for (v, pal) in palette.iter_mut().zip(&palette_channel) {
            v.green = pal.min(1.0).max(0.0);
        }
        // Blue
        let s_k = Matrix2d::from_vec(s.iter().map(|a| a.blue).collect(), s.width(), s.height());
        let r_k = r.iter().map(|a| a.blue).collect::<Vec<_>>();
        let palette_channel = (((s_k * 2.0).matrix_inverse()?) * -1.0) * r_k;
        for (v, pal) in palette.iter_mut().zip(&palette_channel) {
            v.blue = pal.min(1.0).max(0.0);
        }

        Ok(())
    }
}
