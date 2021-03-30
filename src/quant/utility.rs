use rand::seq::SliceRandom;

use crate::matrix::{Matrix2d, Matrix3d};
use crate::{QuantError, SpatialQuant};

pub fn compute_max_coarse_level(mut width: usize, mut height: usize) -> usize {
    const MAX_PIXELS: usize = 4000;
    let mut result = 0;

    while width * height > MAX_PIXELS {
        width >>= 1;
        height >>= 1;
        result += 1;
    }

    result
}

/// Fill a matrix with random f64 from 0 to 1.
pub fn fill_random(arr: &mut Matrix3d<f64>, rng: &mut impl rand::Rng) {
    arr.iter_mut().for_each(|a| {
        *a = rng.gen_range(0.0..=1.0);
    });
}

// Generate random permutation of indices in 2d matrix.
pub fn random_permutation_2d(
    width: usize,
    height: usize,
    result: &mut std::collections::VecDeque<(i32, i32)>,
    rng: &mut impl rand::Rng,
) {
    let mut perm1d = (0..width as i32 * height as i32).collect::<Vec<_>>();
    perm1d.shuffle(rng);

    while let Some(val) = perm1d.pop() {
        result.push_back((val % width as i32, val / width as i32));
    }
}

pub fn compute_b_array<T>(
    filter_weights: &Matrix2d<T>,
    b: &mut Matrix2d<T>,
) -> Result<(), QuantError>
where
    T: Clone + Default + Copy + core::ops::AddAssign<T> + SpatialQuant,
{
    // Assume pixel i is located in the center of b, vary j throughout b
    let radius_width = (filter_weights.width() as isize - 1) / 2;
    let radius_height = (filter_weights.height() as isize - 1) / 2;
    let offset_x = (b.width() as isize - 1) / 2 - radius_width;
    let offset_y = (b.height() as isize - 1) / 2 - radius_height;

    for (j_y, b_chunk) in b.rows_mut().enumerate() {
        for (j_x, b_entry) in b_chunk.iter_mut().enumerate() {
            for (k_y, fw_chunk) in filter_weights.rows().enumerate() {
                for (k_x, fw_entry) in fw_chunk.iter().enumerate() {
                    if k_x as isize + offset_x >= j_x as isize - radius_width
                    && k_x as isize + offset_x <= j_x as isize + radius_width
                    // Original has the following as radius_width which seems
                    // to be an error however the filter width and height are equal
                    && k_y as isize + offset_y >= j_y as isize - radius_height
                    && k_y as isize + offset_y <= j_y as isize + radius_height
                    {
                        *b_entry += fw_entry.direct_product(
                            filter_weights
                                .get(
                                    (k_x as isize + offset_x - j_x as isize + radius_width)
                                        as usize,
                                    (k_y as isize + offset_y - j_y as isize + radius_height)
                                        as usize,
                                )
                                .ok_or("Could not compute b array")?,
                        );
                    }
                }
            }
        }
    }

    Ok(())
}

pub fn b_value<T>(b: &Matrix2d<T>, i_x: isize, i_y: isize, j_x: isize, j_y: isize) -> T
where
    T: Clone + Default + Copy,
{
    let radius_width = b.width().saturating_sub(1) / 2;
    let radius_height = b.height().saturating_sub(1) / 2;
    let k_x = j_x - i_x + radius_width as isize;
    let k_y = j_y - i_y + radius_height as isize;

    *b.get(k_x as usize, k_y as usize).unwrap_or(&T::default())
}

pub fn compute_a_image<T>(
    image: &Matrix2d<T>,
    b: &Matrix2d<T>,
    a: &mut Matrix2d<T>,
) -> Result<(), QuantError>
where
    T: Clone + Default + Copy + core::ops::AddAssign<T> + core::ops::MulAssign<f64> + SpatialQuant,
{
    let radius_width = b.width().saturating_sub(1) / 2;
    let radius_height = b.height().saturating_sub(1) / 2;
    let a_width = a.width();
    let a_height = a.height();

    for (i_y, a_chunk) in a.rows_mut().enumerate() {
        for (i_x, a_entry) in a_chunk.iter_mut().enumerate() {
            for j_y in i_y.saturating_sub(radius_height)..=(i_y + radius_height) {
                if j_y >= a_height {
                    break;
                }

                for j_x in (i_x.saturating_sub(radius_width))..=(i_x + radius_width) {
                    if j_x >= a_width {
                        break;
                    }

                    if let Some(inner) = image.get(j_x, j_y) {
                        *a_entry +=
                            b_value(b, i_x as isize, i_y as isize, j_x as isize, j_y as isize)
                                .direct_product(inner);
                    }
                }
            }
            *a_entry *= -2.0;
        }
    }

    Ok(())
}

pub fn sum_coarsen<T>(fine: &Matrix2d<T>, coarse: &mut Matrix2d<T>) -> Result<(), QuantError>
where
    T: Clone + Default + Copy + core::ops::AddAssign<T> + core::ops::Mul<f64, Output = T>,
{
    for (y, coarse_chunk) in coarse.rows_mut().enumerate() {
        for (x, coarse_entry) in coarse_chunk.iter_mut().enumerate() {
            // coarse width and height should be at most half the fine matrix
            let mut val = *fine
                .get(x * 2, y * 2)
                .ok_or("Could not access fine in sum coarsen")?;
            // let mut divisor = 1.0f64;

            if let Some(v) = fine.get(x * 2 + 1, y * 2) {
                val += *v;
                // divisor += 1.0;
            }
            if let Some(v) = fine.get(x * 2, y * 2 + 1) {
                val += *v;
                // divisor += 1.0;
            }
            if let Some(v) = fine.get(x * 2 + 1, y * 2 + 1) {
                val += *v;
                // divisor += 1.0;
            }

            *coarse_entry = val;
            // NOTE: Original code has divisor term commented out
            // *coarse.get_mut(x, y).ok_or("") = val * (divisor.recip());
        }
    }
    Ok(())
}

pub fn best_match_color<T>(
    vars: &Matrix3d<f64>,
    i_x: usize,
    i_y: usize,
    palette: &[T],
) -> Result<usize, QuantError> {
    let mut max_v = 0;
    let mut max_weight = *vars
        .get(i_x, i_y, 0)
        .ok_or("Could not compute best match color")?;
    for v in 1..palette.len() {
        if let Some(var) = vars.get(i_x, i_y, v) {
            if *var > max_weight {
                max_v = v;
                max_weight = *var;
            }
        }
    }

    Ok(max_v)
}

pub fn zoom_double(small: &Matrix3d<f64>, big: &mut Matrix3d<f64>) -> Result<(), QuantError> {
    // Scale the weights array by mixing the pixels under each fine pixel,
    // weighted by area. The fine pixel is assumed to be 1.2 x 1.2 pixels.

    for y in 0..(big.height() / 2 * 2) {
        for x in 0..(big.width() / 2 * 2) {
            let left = ((x as f64 - 0.1) / 2.0).max(0.0);
            let right = (small.width() as f64 - 0.001).min((x as f64 + 1.1) / 2.0);
            let top = ((y as f64 - 0.1) / 2.0).max(0.0);
            let bottom = (small.height() as f64 - 0.001).min((y as f64 + 1.1) / 2.0);
            let x_left = left.floor() as usize;
            let x_right = right.floor() as usize;
            let y_top = top.floor() as usize;
            let y_bottom = bottom.floor() as usize;
            let area = (right - left) * (bottom - top);

            let top_weight = (right - left) * (top.ceil() - top) / area;
            let bottom_weight = (right - left) * (bottom - bottom.floor()) / area;
            let left_weight = (bottom - top) * (left.ceil() - left) / area;
            let right_weight = (bottom - top) * (right - right.floor()) / area;
            let top_left_weight = (left.ceil() - left) * (top.ceil() - top) / area;
            let top_right_weight = (right - right.floor()) * (top.ceil() - top) / area;
            let bottom_left_weight = (left.ceil() - left) * (bottom - bottom.floor()) / area;
            let bottom_right_weight = (right - right.floor()) * (bottom - bottom.floor()) / area;

            for z in 0..big.depth() {
                if x_left == x_right && y_top == y_bottom {
                    *big.get_mut(x, y, z).ok_or("Could not access big matrix")? =
                        *small.get(x_left, y_top, z).ok_or("Could not zoom double")?;
                } else if x_left == x_right {
                    *big.get_mut(x, y, z).ok_or("Could not access big matrix")? = top_weight
                        * *small.get(x_left, y_top, z).ok_or("Could not zoom double")?
                        + bottom_weight
                            * *small
                                .get(x_left, y_bottom, z)
                                .ok_or("Could not zoom double")?;
                } else if y_top == y_bottom {
                    *big.get_mut(x, y, z).ok_or("Could not access big matrix")? = left_weight
                        * *small.get(x_left, y_top, z).ok_or("Could not zoom double")?
                        + right_weight
                            * *small
                                .get(x_right, y_top, z)
                                .ok_or("Could not zoom double")?;
                } else {
                    *big.get_mut(x, y, z).ok_or("Could not access big matrix")? = top_left_weight
                        * *small.get(x_left, y_top, z).ok_or("Could not zoom double")?
                        + top_right_weight
                            * *small
                                .get(x_right, y_top, z)
                                .ok_or("Could not zoom double")?
                        + bottom_left_weight
                            * *small
                                .get(x_left, y_bottom, z)
                                .ok_or("Could not zoom double")?
                        + bottom_right_weight
                            * *small
                                .get(x_right, y_bottom, z)
                                .ok_or("Could not zoom double")?;
                }
            }
        }
    }

    Ok(())
}

pub fn compute_initial_s<T>(
    s: &mut Matrix2d<T>,
    coarse_variables: &Matrix3d<f64>,
    b: &Matrix2d<T>,
) -> Result<(), QuantError>
where
    T: Clone + Default + Copy + core::ops::AddAssign<T> + core::ops::Mul<f64, Output = T>,
{
    let palette_size = s.width();
    let coarse_width = coarse_variables.width();
    let coarse_height = coarse_variables.height();
    let center_x = b.width().saturating_sub(1) / 2;
    let center_y = b.height().saturating_sub(1) / 2;
    let center_b = b_value(&b, 0, 0, 0, 0);
    for v in 0..palette_size {
        for alpha in v..palette_size {
            *s.get_mut(v, alpha)
                .ok_or("Could not access s to compute initial s")? = T::default();
        }
    }

    for i_y in 0..coarse_height {
        for i_x in 0..coarse_width {
            let max_j_x = coarse_width.min(i_x.saturating_sub(center_x) + b.width());
            let max_j_y = coarse_height.min(i_y.saturating_sub(center_y) + b.height());

            for j_y in i_y.saturating_sub(center_y)..max_j_y {
                for j_x in i_x.saturating_sub(center_x)..max_j_x {
                    if i_x == j_x && i_y == j_y {
                        continue;
                    }
                    let b_ij = b_value(b, i_x as isize, i_y as isize, j_x as isize, j_y as isize);
                    for v in 0..palette_size {
                        for alpha in v..palette_size {
                            let mult = coarse_variables
                                .get(i_x as usize, i_y as usize, v)
                                .ok_or("Could not compute initial s")?
                                * coarse_variables.get(j_x, j_y, alpha).ok_or(
                                    "Could not access coarse variables to compute initial s",
                                )?;
                            *s.get_mut(v, alpha)
                                .ok_or("Could not access s to compute initial s")? += b_ij * mult;
                        }
                    }
                }
            }
            for v in 0..palette_size {
                *s.get_mut(v, v)
                    .ok_or("Could not access s to compute initial s")? += center_b
                    * *coarse_variables
                        .get(i_x as usize, i_y as usize, v)
                        .ok_or("Could not access coarse variables to compute initial s")?;
            }
        }
    }

    Ok(())
}

pub fn update_s<T>(
    s: &mut Matrix2d<T>,
    coarse_variables: &Matrix3d<f64>,
    b: &Matrix2d<T>,
    j_x: usize,
    j_y: usize,
    alpha: usize,
    delta: f64,
) -> Result<(), QuantError>
where
    T: Clone + Default + Copy + core::ops::AddAssign<T> + core::ops::Mul<f64, Output = T>,
{
    let palette_size = s.width();
    let center_x = b.width().saturating_sub(1) / 2;
    let center_y = b.height().saturating_sub(1) / 2;
    let max_i_x = coarse_variables.width().min(j_x + center_x + 1);
    let max_i_y = coarse_variables.height().min(j_y + center_y + 1);

    for i_y in j_y.saturating_sub(center_y)..max_i_y {
        for i_x in j_x.saturating_sub(center_x)..max_i_x {
            let delta_b_ij =
                b_value(b, i_x as isize, i_y as isize, j_x as isize, j_y as isize) * delta;
            if i_x == j_x && i_y == j_y {
                continue;
            }

            for v in 0..=alpha {
                let mult = *coarse_variables
                    .get(i_x, i_y, v)
                    .ok_or("Could not access coarse variables to update s")?;
                *s.get_mut(v, alpha)
                    .ok_or("Could not access s to update s")? += delta_b_ij * mult;
            }
            for v in alpha..palette_size {
                let mult = *coarse_variables
                    .get(i_x, i_y, v)
                    .ok_or("Could not access coarse variables to update s")?;
                *s.get_mut(alpha, v)
                    .ok_or("Could not access s to update s")? += delta_b_ij * mult;
            }
        }
    }
    *s.get_mut(alpha, alpha)
        .ok_or("Could not access s to update s")? += b_value(b, 0, 0, 0, 0) * delta;

    Ok(())
}

pub fn compute_initial_j_palette_sum<T>(
    j_palette_sum: &mut Matrix2d<T>,
    coarse_variables: &Matrix3d<f64>,
    palette: &[T],
) -> Result<(), QuantError>
where
    T: Clone + Default + Copy + core::ops::AddAssign<T> + core::ops::Mul<f64, Output = T>,
{
    for (j_y, j_palette_chunk) in j_palette_sum.rows_mut().enumerate() {
        for (j_x, j_palette_entry) in j_palette_chunk.iter_mut().enumerate() {
            let mut palette_sum = T::default();
            for (alpha, &item) in palette.iter().enumerate() {
                palette_sum += item
                    * *coarse_variables
                        .get(j_x, j_y, alpha)
                        .ok_or("Could not access coarse variables in j palette sum")?;
            }

            *j_palette_entry = palette_sum;
        }
    }

    Ok(())
}
