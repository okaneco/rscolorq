mod utility;

use rand::SeedableRng;

use crate::matrix::{Matrix2d, Matrix3d};
use crate::SpatialQuant;
use utility::*;

/// Square filter size for the quantization.
///
/// Note that `Three` is the default. `Five` can take a long time to complete.
/// `One` does not produce an image resembling the original input but is
/// included for completion's sake.
#[derive(Clone, Copy, Debug)]
pub enum FilterSize {
    /// 1x1 pixel filter. This filter does not produce a usable resulting image.
    One,
    /// 3x3 pixel filter.
    Three,
    /// 5x5 pixel filter, may be quite slow.
    Five,
}

impl Default for FilterSize {
    fn default() -> Self {
        FilterSize::Three
    }
}

/// Input parameter struct for spatial color quantization and simulated
/// annealing. The parameters can be validated with
/// [`verify_parameters`][verify] before quantization.
///
/// [verify]: struct.Params.html#method.verify_parameters
///
/// If a palette is supplied to the struct, the dithering will be performed with
/// fixed colors and will override the `palette_size` field.
#[derive(Clone, Debug)]
pub struct Params<T: SpatialQuant> {
    initial_temp: f64,
    final_temp: f64,
    palette_size: u8,
    iters_per_level: usize,
    repeats_per_temp: usize,
    seed: Option<u64>,
    filter_size: FilterSize,
    dithering_level: f64,
    palette: Option<Vec<T>>,
}

impl<T: SpatialQuant> Params<T> {
    /// Crate a new input parameter struct.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the initial temperature. This must be greater than 0 and greater
    /// than the final temperature.
    pub fn initial_temp(&mut self, initial_temp: f64) -> &mut Self {
        self.initial_temp = initial_temp;
        self
    }

    /// Set the final temperature. This value must be greater than 0 and less
    /// than the initial temperature.
    pub fn final_temp(&mut self, final_temp: f64) -> &mut Self {
        self.final_temp = final_temp;
        self
    }

    /// Set the desired palette size. This parameter must be greater than 2.
    pub fn palette_size(&mut self, palette_size: u8) -> &mut Self {
        self.palette_size = palette_size;
        self
    }

    /// Set the iterations per level. This parameter must be greater than 0.
    pub fn iters_per_level(&mut self, iters_per_level: usize) -> &mut Self {
        self.iters_per_level = iters_per_level;
        self
    }

    /// Set the amount of repeats per temperature. This parameter must be greater
    /// than 0.
    pub fn repeats_per_temp(&mut self, repeats_per_temp: usize) -> &mut Self {
        self.repeats_per_temp = repeats_per_temp;
        self
    }

    /// Set the random number generator seed.
    pub fn seed(&mut self, seed: Option<u64>) -> &mut Self {
        self.seed = seed;
        self
    }

    /// Set the filter size for quantization. The default is 3.
    pub fn filter_size(&mut self, filter_size: FilterSize) -> &mut Self {
        self.filter_size = filter_size;
        self
    }

    /// Set the dithering level manually, default is 0.8. The range should be
    /// between roughly 0.5 and 1.5, higher numbers will result in higher
    /// frequency noise.
    ///
    /// Dithering level must be greater than 0.
    pub fn dithering_level(&mut self, dithering_level: f64) -> &mut Self {
        self.dithering_level = dithering_level;
        self
    }

    /// Calculate and set the dithering level automatically based on image size
    /// and palette size.
    ///
    /// Dithering level must be greater than 0.
    pub fn dithering_level_auto(
        &mut self,
        width: u32,
        height: u32,
        palette_size: usize,
    ) -> &mut Self {
        self.dithering_level =
            0.09 * ((width * height) as f64).ln() - 0.04 * (palette_size as f64).ln() + 0.001;
        self
    }

    /// Supply the palette colors to use for dithering. Length must be at least
    /// two.
    pub fn palette(&mut self, palette: Vec<T>) -> Result<&mut Self, crate::QuantError> {
        if palette.len() < 2 {
            return Err(crate::QuantError::Parameter(String::from(
                "Palette size must be at least 2.",
            )));
        }
        self.palette = Some(palette);
        Ok(self)
    }

    /// Returns an error if any of the parameters are incorrectly set.
    pub fn verify_parameters(&self) -> Result<(), crate::QuantError> {
        if self.palette.is_none() && self.palette_size < 2 {
            return Err(crate::QuantError::Parameter(String::from(
                "Palette size must be at least 2.",
            )));
        }
        if self.iters_per_level < 1 {
            return Err(crate::QuantError::Parameter(String::from(
                "iters_per_level must be greater than 0.",
            )));
        }
        if self.repeats_per_temp < 1 {
            return Err(crate::QuantError::Parameter(String::from(
                "repeats_per_temp must be greater than 0.",
            )));
        }
        if self.dithering_level <= 0.0 {
            return Err(crate::QuantError::Parameter(String::from(
                "Dithering level must be greater than 0.",
            )));
        }
        if self.initial_temp <= self.final_temp {
            return Err(crate::QuantError::Parameter(String::from(
                "Initial temperature must be greater than final temperature.",
            )));
        }
        if self.initial_temp <= 0.0 || self.final_temp <= 0.0 {
            return Err(crate::QuantError::Parameter(String::from(
                "Temperatures must be greater than 0.",
            )));
        }

        Ok(())
    }
}

impl<T: SpatialQuant> Default for Params<T> {
    fn default() -> Self {
        Self {
            initial_temp: 1.0,
            final_temp: 0.001,
            palette_size: 4,
            iters_per_level: 3,
            repeats_per_temp: 1,
            seed: None,
            filter_size: FilterSize::default(),
            dithering_level: 0.8,
            palette: None,
        }
    }
}

/// Perform the spatial color quantization. The mutated input parameters are
/// a 2-dimensional matrix with the quantized color palette indices and the
/// quantized color palette.
pub fn spatial_color_quant<T>(
    image: &Matrix2d<T>,
    quantized_image: &mut Matrix2d<u8>,
    palette: &mut Vec<T>,
    conditions: &Params<T>,
) -> Result<(), crate::QuantError>
where
    T: Clone
        + Default
        + Copy
        + core::ops::Add<T, Output = T>
        + core::ops::AddAssign<T>
        + core::ops::Mul<f64, Output = T>
        + core::ops::Mul<T, Output = T>
        + core::ops::MulAssign<f64>
        + core::ops::Sub<T, Output = T>
        + SpatialQuant,
{
    let mut rng = conditions.seed.map_or_else(
        rand_pcg::Pcg64Mcg::from_entropy,
        rand_pcg::Pcg64Mcg::seed_from_u64,
    );

    // Initialize palette
    if conditions.palette.is_none() {
        for _ in 0..conditions.palette_size {
            palette.push(T::random(&mut rng));
        }
    } else {
        conditions
            .palette
            .clone()
            .ok_or("Could not access conditions palette")?
            .iter()
            .for_each(|&a| palette.push(a));
    }

    let max_coarse_level = compute_max_coarse_level(image.width(), image.height());
    let mut coarse_variables = Matrix3d::<f64>::new(
        image.width() >> max_coarse_level,
        image.height() >> max_coarse_level,
        palette.len(),
    );
    fill_random(&mut coarse_variables, &mut rng);

    let filter_weights =
        T::calculate_filter_weights(conditions.dithering_level, conditions.filter_size);
    let mut temperature = conditions.initial_temp;

    // Compute a_i, b_{ij} according to (11)
    let extended_neighborhood_width = filter_weights.width() * 2 - 1;
    let extended_neighborhood_height = filter_weights.height() * 2 - 1;
    let mut b0 = Matrix2d::new(extended_neighborhood_width, extended_neighborhood_height);
    compute_b_array(&filter_weights, &mut b0)?;

    let mut a0 = Matrix2d::new(image.width(), image.height());
    compute_a_image(&image, &b0, &mut a0)?;

    // Compute a_I^l, b_{IJ}^l according to (18)
    let mut a_vec = Vec::new();
    let mut b_vec = Vec::new();
    a_vec.push(a0);
    b_vec.push(b0);

    for coarse_level in 1..=max_coarse_level {
        let radius_width = filter_weights.width().saturating_sub(1) / 2;
        let radius_height = filter_weights.height().saturating_sub(1) / 2;
        let mut bi = Matrix2d::new(
            (3).max(
                b_vec
                    .last()
                    .ok_or("Could not get last in b_vec")?
                    .width()
                    .saturating_sub(2),
            ),
            (3).max(
                b_vec
                    .last()
                    .ok_or("Could not get last in b_vec")?
                    .height()
                    .saturating_sub(2),
            ),
        );

        for big_j_y in 0..bi.height() {
            for big_j_x in 0..bi.width() {
                for i_y in (radius_height * 2)..(radius_height * 2 + 2) {
                    for i_x in (radius_width * 2)..(radius_width * 2 + 2) {
                        for j_y in (big_j_y * 2)..(big_j_y * 2 + 2) {
                            for j_x in (big_j_x * 2)..(big_j_x * 2 + 2) {
                                *bi.get_mut(big_j_x, big_j_y).ok_or("Could not access bi")? +=
                                    b_value(
                                        b_vec
                                            .last()
                                            .ok_or("Could not calculate b_value with b_vec")?,
                                        i_x as isize,
                                        i_y as isize,
                                        j_x as isize,
                                        j_y as isize,
                                    );
                            }
                        }
                    }
                }
            }
        }

        b_vec.push(bi);

        let mut ai = Matrix2d::new(
            image.width() >> coarse_level,
            image.height() >> coarse_level,
        );
        sum_coarsen(
            a_vec
                .last()
                .ok_or("Could not access a_vec for sum_coarsen")?,
            &mut ai,
        )?;
        a_vec.push(ai);
    }

    // Multiscale annealing
    let mut coarse_level = max_coarse_level as isize;
    let temp_multiplier = (conditions.final_temp / conditions.initial_temp)
        .powf(1.0 / (3.0f64).max(max_coarse_level as f64 * conditions.iters_per_level as f64));
    let mut iters_at_current_level = 0;
    let mut skip_palette_maintenance = false;

    let mut s = Matrix2d::new(palette.len(), palette.len());
    compute_initial_s(
        &mut s,
        &coarse_variables,
        b_vec
            .get(coarse_level as usize)
            .ok_or("Could not access b_vec for computing s")?,
    )?;
    let mut j_palette_sum = Matrix2d::new(coarse_variables.width(), coarse_variables.height());
    compute_initial_j_palette_sum(&mut j_palette_sum, &coarse_variables, &palette)?;

    while coarse_level >= 0 || temperature > conditions.final_temp {
        let a = a_vec
            .get(coarse_level as usize)
            .ok_or("Could not access a_vec")?;
        let b = b_vec
            .get(coarse_level as usize)
            .ok_or("Could not access b_vec")?;
        let middle_b = b_value(b, 0, 0, 0, 0);
        let center_x = (b.width().saturating_sub(1) / 2) as i32;
        let center_y = (b.height().saturating_sub(1) / 2) as i32;

        // let mut steps = 0u32;
        for _ in 0..conditions.repeats_per_temp {
            let mut visit_queue = std::collections::VecDeque::new();
            random_permutation_2d(
                coarse_variables.width(),
                coarse_variables.height(),
                &mut visit_queue,
                &mut rng,
            );

            // Compute 2*sum(j in extended neighborhood of i, j != i) b_ij
            while !visit_queue.is_empty() {
                // Revisit everything if 10% above initial size
                if visit_queue.len()
                    > coarse_variables.width() * coarse_variables.height() * 11 / 10
                {
                    visit_queue.clear();
                    random_permutation_2d(
                        coarse_variables.width(),
                        coarse_variables.height(),
                        &mut visit_queue,
                        &mut rng,
                    );
                }

                let (i_x, i_y) = visit_queue
                    .pop_front()
                    .ok_or("Could not pop from visit queue")?;

                // Compute (25)
                let mut p_i = T::default();
                for y in 0..b.height() {
                    for x in 0..b.width() {
                        let j_x = x as i32 - center_x + i_x;
                        let j_y = y as i32 - center_y + i_y;
                        if i_x == j_x && i_y == j_y {
                            continue;
                        }
                        if j_x < 0
                            || j_y < 0
                            || j_x as usize >= coarse_variables.width()
                            || j_y as usize >= coarse_variables.height()
                        {
                            continue;
                        }
                        let b_ij =
                            b_value(&b, i_x as isize, i_y as isize, j_x as isize, j_y as isize);
                        let j_pal = *j_palette_sum
                            .get(j_x as usize, j_y as usize)
                            .ok_or("Could not access j_palette_sum")?;
                        p_i += b_ij * j_pal;
                    }
                }
                p_i *= 2.0;
                if let Some(a) = a.get(i_x as usize, i_y as usize) {
                    p_i += *a;
                }

                let mut meanfield_logs = Vec::new();
                let mut meanfields = Vec::new();
                let mut max_meanfield_log = core::f64::NEG_INFINITY;
                let mut meanfield_sum = 0.0;
                for v in palette.iter() {
                    // Update m_{pi(i)v}^I according to (23)
                    meanfield_logs
                        .push(-(v.dot_product(&(p_i + middle_b.direct_product(&v)))) / temperature);
                    if let Some(meanfield_logs_last) = meanfield_logs.last() {
                        if *meanfield_logs_last > max_meanfield_log {
                            max_meanfield_log = *meanfield_logs_last;
                        }
                    }
                }
                for v in meanfield_logs.iter().take(palette.len()) {
                    meanfields.push((v - max_meanfield_log + 100.0).exp());
                    meanfield_sum += meanfields
                        .last()
                        .ok_or("Could not access last in meanfields")?;
                }

                if !meanfield_sum.is_normal() {
                    return Err(crate::QuantError::Quantization(String::from(
                        "Meanfield sum underflowed",
                    )));
                }

                let old_max_v =
                    best_match_color(&coarse_variables, i_x as usize, i_y as usize, palette)?;
                for (v, &palette_item) in palette.iter().enumerate() {
                    let mut new_val =
                        meanfields.get(v).ok_or("Could not access meanfields")? / meanfield_sum;
                    // Prevent S from becoming singular
                    if new_val <= 0.0 {
                        new_val = 1e-10;
                    }
                    if new_val >= 1.0 {
                        new_val = 1.0 - 1e-10;
                    }
                    let delta_m_iv = new_val
                        - coarse_variables
                            .get(i_x as usize, i_y as usize, v)
                            .ok_or("Could not access coarse_variables for delta v")?;
                    *coarse_variables
                        .get_mut(i_x as usize, i_y as usize, v)
                        .ok_or("Could not access coarse_variables for delta v")? = new_val;
                    *j_palette_sum
                        .get_mut(i_x as usize, i_y as usize)
                        .ok_or("Could not access j_palette_sum for delta v")? +=
                        palette_item * delta_m_iv;
                    if delta_m_iv.abs() > 0.001 && !skip_palette_maintenance {
                        update_s(
                            &mut s,
                            &coarse_variables,
                            &b,
                            i_x as usize,
                            i_y as usize,
                            v,
                            delta_m_iv,
                        )?;
                    }
                }

                let max_v =
                    best_match_color(&coarse_variables, i_x as usize, i_y as usize, palette)?;
                // Color difference delta, determines if pixel has changed
                if palette
                    .get(max_v)
                    .ok_or("Could not access palette for max_v")?
                    .color_difference(
                        palette
                            .get(old_max_v)
                            .ok_or("Could not access palette for old_max_v")?,
                    )
                    >= T::difference_threshold()
                {
                    for y in (1).min(center_y - 1)
                        ..(b.height().saturating_sub(1) as i32).max(center_y + 1)
                    {
                        for x in (1).min(center_x - 1)
                            ..(b.width().saturating_sub(1) as i32).max(center_x + 1)
                        {
                            let j_x = x - center_x + i_x as i32;
                            let j_y = y - center_y + i_y as i32;
                            if j_x < 0
                                || j_y < 0
                                || j_x >= coarse_variables.width() as i32
                                || j_y >= coarse_variables.height() as i32
                            {
                                continue;
                            }
                            visit_queue.push_back((j_x, j_y));
                        }
                    }
                }

                // steps += 1;
                // if steps % 10_000 == 0 {
                //     print!("{},", visit_queue.len());
                //     std::io::Write::flush(&mut std::io::stdout()).unwrap();
                // }
            }

            if skip_palette_maintenance {
                compute_initial_s(
                    &mut s,
                    &coarse_variables,
                    b_vec
                        .get(coarse_level as usize)
                        .ok_or("Could not access b_vec in skip maintenance")?,
                )?;
            }
            if conditions.palette.is_none() {
                T::refine_palette(&mut s, &coarse_variables, &a, palette)?;
            }
            compute_initial_j_palette_sum(&mut j_palette_sum, &coarse_variables, &palette)?;
        }

        iters_at_current_level += 1;
        skip_palette_maintenance = false;
        if (temperature <= conditions.final_temp || coarse_level > 0)
            && iters_at_current_level >= conditions.iters_per_level
        {
            coarse_level -= 1;
            if coarse_level < 0 {
                break;
            };
            let mut new_coarse_variables = Matrix3d::new(
                image.width() >> coarse_level,
                image.height() >> coarse_level,
                palette.len(),
            );
            zoom_double(&coarse_variables, &mut new_coarse_variables)?;
            coarse_variables = new_coarse_variables;
            iters_at_current_level = 0;
            j_palette_sum = Matrix2d::new(coarse_variables.width(), coarse_variables.height());
            compute_initial_j_palette_sum(&mut j_palette_sum, &coarse_variables, &palette)?;
            skip_palette_maintenance = true;
        }

        if temperature > conditions.final_temp {
            temperature *= temp_multiplier;
        }
    }

    while coarse_level > 0 {
        coarse_level -= 1;
        let mut new_coarse_variables = Matrix3d::new(
            image.width() >> coarse_level,
            image.height() >> coarse_level,
            palette.len(),
        );
        zoom_double(&coarse_variables, &mut new_coarse_variables)?;
        coarse_variables = new_coarse_variables;
    }

    for i_y in 0..image.height() {
        for i_x in 0..image.width() {
            *quantized_image
                .get_mut(i_x, i_y)
                .ok_or("Could not access quantized image")? =
                best_match_color(&coarse_variables, i_x, i_y, palette)? as u8;
        }
    }

    Ok(())
}
