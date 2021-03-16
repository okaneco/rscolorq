use palette::{white_point::D65, Lab, Pixel, Srgb};
use rscolorq::{color::Rgb, spatial_color_quant, FilterSize, Matrix2d, Params};
use structopt::StructOpt;

mod args;
mod utils;

use utils::{parse_color, print_colors, print_colors_lab};
use utils::{save_image, save_palette, save_palette_lab};

fn main() {
    if let Err(e) = try_main() {
        eprintln!("{}", e);
        std::process::exit(1);
    }
}

fn try_main() -> Result<(), Box<dyn std::error::Error>> {
    let opt = args::Opt::from_args();

    if opt.output.is_none() && opt.palette_output.is_none() && !opt.print && !opt.no_file {
        return Err("No output specified.".into());
    }

    let img = image::open(opt.input)?.to_rgb8();
    let (width, height) = (img.dimensions().0, img.dimensions().1);

    let dithering_level = if opt.dither <= 0.0 {
        eprintln!("Dithering level must be greater than 0, setting to 0.8");
        0.8
    } else {
        opt.dither
    };

    let filter_size = match opt.filter {
        1 => FilterSize::One,
        3 => FilterSize::Three,
        5 => FilterSize::Five,
        _ => {
            eprintln!("Filter must be 3 or 5; setting to 3");
            FilterSize::Three
        }
    };

    let mut imgbuf = Vec::with_capacity(width as usize * height as usize * 3);
    let mut quantized_image = Matrix2d::new(width as usize, height as usize);

    if !opt.lab {
        // Set quantization parameters
        let mut conditions = Params::new();
        conditions.initial_temp(opt.start);
        conditions.final_temp(opt.end);
        conditions.iters_per_level(opt.iters);
        conditions.repeats_per_temp(opt.repeats);
        conditions.seed(opt.seed);
        conditions.filter_size(filter_size);
        conditions.verify_parameters()?;

        // Determine user supplied colors, or set palette size
        let mut palette;
        if !opt.colors.is_empty() {
            let mut color_vec = Vec::with_capacity(opt.colors.len());
            for c in &opt.colors {
                let color = parse_color(c.trim_start_matches('#'))?;
                color_vec.push(Rgb {
                    red: f64::from(color[0]) / 255.0,
                    green: f64::from(color[1]) / 255.0,
                    blue: f64::from(color[2]) / 255.0,
                });
            }
            palette = Vec::with_capacity(color_vec.len());

            // set dithering level if automatic
            if opt.auto {
                conditions.dithering_level_auto(width, height, color_vec.len());
            } else {
                conditions.dithering_level(dithering_level);
            }
            conditions.palette(color_vec)?;
        } else {
            if opt.auto {
                conditions.dithering_level_auto(width, height, opt.n as usize);
            } else {
                conditions.dithering_level(dithering_level);
            }
            conditions.palette_size(opt.n);
            palette = Vec::with_capacity(opt.n as usize);
        }

        // Convert image to Rgb<f64>
        let image = Matrix2d::from_vec(
            img.pixels()
                .map(|&c| Rgb {
                    red: f64::from(c.0[0]) / 255.0,
                    green: f64::from(c.0[1]) / 255.0,
                    blue: f64::from(c.0[2]) / 255.0,
                })
                .collect(),
            width as usize,
            height as usize,
        );

        spatial_color_quant(&image, &mut quantized_image, &mut palette, &conditions)?;

        // Convert the Rgb<f64> palette to Rgb<u8>
        let palette = palette
            .iter()
            .map(|&c| {
                let color = 255.0 * c;
                [
                    color.red.round() as u8,
                    color.green.round() as u8,
                    color.blue.round() as u8,
                ]
            })
            .collect::<Vec<[u8; 3]>>();

        // Print the colors to console
        if opt.print {
            print_colors(&palette)?;
        }

        // Save a palette image if the palette_output is Some
        if let Some(title) = opt.palette_output {
            save_palette(&palette, opt.height, opt.width, &title)?;
        }

        // Create the final image by color lookup from the palette
        for &c in quantized_image.iter() {
            let color = palette
                .get(c as usize)
                .ok_or("Could not retrieve color from palette")?;
            imgbuf.extend_from_slice(color);
        }
    } else {
        // Set quantization parameters
        let mut conditions = Params::new();
        conditions.initial_temp(opt.start);
        conditions.final_temp(opt.end);
        conditions.iters_per_level(opt.iters);
        conditions.repeats_per_temp(opt.repeats);
        conditions.seed(opt.seed);
        conditions.filter_size(filter_size);
        conditions.dithering_level(dithering_level);
        conditions.verify_parameters()?;

        // Determine user supplied colors, or set palette size
        let mut palette;
        if !opt.colors.is_empty() {
            let mut color_vec = Vec::with_capacity(opt.colors.len());
            for c in &opt.colors {
                let color = parse_color(c.trim_start_matches('#'))?;
                color_vec.push(Srgb::new(color[0], color[1], color[2]).into_format().into());
            }
            palette = Vec::with_capacity(color_vec.len());

            // set dithering level if automatic
            if opt.auto {
                conditions.dithering_level_auto(width, height, color_vec.len());
            } else {
                conditions.dithering_level(dithering_level);
            }
            conditions.palette(color_vec)?;
        } else {
            if opt.auto {
                conditions.dithering_level_auto(width, height, opt.n as usize);
            } else {
                conditions.dithering_level(dithering_level);
            }
            conditions.palette_size(opt.n);
            palette = Vec::with_capacity(opt.n as usize);
        }

        // Convert image to Lab<f64>
        let img_vec = img.into_raw();
        let image = Matrix2d::from_vec(
            Srgb::from_raw_slice(&img_vec)
                .iter()
                .map(|x| x.into_format().into())
                .collect::<Vec<Lab<D65, f64>>>(),
            width as usize,
            height as usize,
        );

        spatial_color_quant(&image, &mut quantized_image, &mut palette, &conditions)?;

        // Convert the Lab<f64> palette to Srgb<u8>
        let palette = palette
            .iter()
            .map(|&c| Srgb::from(c).into_format())
            .collect::<Vec<Srgb<u8>>>();

        // Print the colors to console
        if opt.print {
            print_colors_lab(&palette)?;
        }

        // Save a palette image if the palette_output is Some
        if let Some(title) = opt.palette_output {
            save_palette_lab(&palette, opt.height, opt.width, &title)?;
        }

        // Create the final image by color lookup from the palette
        for &c in quantized_image.iter() {
            let color = *palette
                .get(c as usize)
                .ok_or("Could not retrieve color from palette")?;
            imgbuf.extend_from_slice(Srgb::into_raw_slice(&[color]));
        }
    }

    if !opt.no_file {
        if let Some(title) = opt.output {
            save_image(&title, &imgbuf, width, height)?;
        }
    }

    Ok(())
}
