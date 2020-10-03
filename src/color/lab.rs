use core::marker::PhantomData;

use crate::{quant::FilterSize, Matrix2d, QuantError, SpatialQuant};
use palette::{white_point::D65, Lab};

impl SpatialQuant for Lab<D65, f64> {
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
                            Lab::new(val, val, val);
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
        // Using the 1976 delta E*_ab for just noticeable difference,
        // 2.3 seemed too high to produce good results
        1.0
    }

    fn direct_product(&self, other: &Self) -> Self {
        Self {
            l: self.l * other.l,
            a: self.a * other.a,
            b: self.b * other.b,
            white_point: PhantomData,
        }
    }

    fn dot_product(&self, other: &Self) -> f64 {
        self.l * other.l + self.a * other.a + self.b * other.b
    }

    fn one() -> Self {
        Self::new(1.0, 1.0, 1.0)
    }

    fn norm_squared(&self) -> f64 {
        self.l * self.l + self.a * self.a + self.b * self.b
    }

    fn random(rng: &mut impl rand::Rng) -> Self {
        Self {
            l: rng.gen_range(0.0, 100.0),
            a: rng.gen_range(-128.0, 127.0),
            b: rng.gen_range(-128.0, 127.0),
            white_point: PhantomData,
        }
    }

    fn refine_palette(
        s: &mut Matrix2d<Self>,
        coarse_variables: &crate::Matrix3d<f64>,
        a: &Matrix2d<Self>,
        palette: &mut Vec<Self>,
    ) -> Result<(), QuantError> {
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
                        *a.get(i_x, i_y).unwrap() * *coarse_variables.get(i_x, i_y, v).unwrap();
                }
            }
        }

        // Iterate through each of the component channels, clamp to valid range
        // Lightness
        let s_k = Matrix2d::from_vec(s.iter().map(|a| a.l).collect(), s.width(), s.height());
        let r_k = r.iter().map(|a| a.l).collect::<Vec<_>>();
        let palette_channel = (((s_k * 2.0).matrix_inverse()?) * -1.0) * r_k;
        for (v, pal) in palette.iter_mut().zip(&palette_channel) {
            v.l = pal.min(100.0).max(0.0);
        }
        // a star
        let s_k = Matrix2d::from_vec(s.iter().map(|a| a.a).collect(), s.width(), s.height());
        let r_k = r.iter().map(|a| a.a).collect::<Vec<_>>();
        let palette_channel = (((s_k * 2.0).matrix_inverse()?) * -1.0) * r_k;
        for (v, pal) in palette.iter_mut().zip(&palette_channel) {
            v.a = pal.min(127.0).max(-128.0);
        }
        // b star
        let s_k = Matrix2d::from_vec(s.iter().map(|a| a.b).collect(), s.width(), s.height());
        let r_k = r.iter().map(|a| a.b).collect::<Vec<_>>();
        let palette_channel = (((s_k * 2.0).matrix_inverse()?) * -1.0) * r_k;
        for (v, pal) in palette.iter_mut().zip(&palette_channel) {
            v.b = pal.min(127.0).max(-128.0);
        }

        Ok(())
    }
}
