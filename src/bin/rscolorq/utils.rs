use std::fmt::Write;

/// Parse hex string to Rgb color.
pub fn parse_color(c: &str) -> Result<[u8; 3], Box<dyn std::error::Error>> {
    let red = u8::from_str_radix(
        match &c.get(0..2) {
            Some(x) => x,
            None => {
                eprintln!("Invalid color: {}", c);
                return Err("Could not parse hex color length, must be 6 characters.".into());
            }
        },
        16,
    )?;
    let green = u8::from_str_radix(
        match &c.get(2..4) {
            Some(x) => x,
            None => {
                eprintln!("Invalid color: {}", c);
                return Err("Could not parse hex color length, must be 6 characters.".into());
            }
        },
        16,
    )?;
    let blue = u8::from_str_radix(
        match &c.get(4..6) {
            Some(x) => x,
            None => {
                eprintln!("Invalid color: {}", c);
                return Err("Could not parse hex color length, must be 6 characters.".into());
            }
        },
        16,
    )?;
    Ok([red, green, blue])
}

/// Prints the palette colors in hexadecimal format.
pub fn print_colors(colors: &[[u8; 3]]) -> Result<(), Box<dyn std::error::Error>> {
    let mut col = String::new();
    if let Some((last, elements)) = colors.split_last() {
        for elem in elements {
            write!(&mut col, "{:02x}{:02x}{:02x},", elem[0], elem[1], elem[2])?;
        }
        writeln!(&mut col, "{:02x}{:02x}{:02x}", last[0], last[1], last[2])?;
    }
    print!("{}", col);

    Ok(())
}

/// Prints the palette colors in hexadecimal format from Srgb.
pub fn print_colors_lab(colors: &[palette::Srgb<u8>]) -> Result<(), Box<dyn std::error::Error>> {
    let mut col = String::new();
    if let Some((last, elements)) = colors.split_last() {
        for elem in elements {
            write!(&mut col, "{:x},", elem)?;
        }
        writeln!(&mut col, "{:x}", last)?;
    }
    print!("{}", col);

    Ok(())
}

/// Saves image buffer to file.
pub fn save_image(
    output: &std::path::PathBuf,
    imgbuf: &[u8],
    width: u32,
    height: u32,
) -> Result<(), Box<dyn std::error::Error>> {
    let w = std::io::BufWriter::new(std::fs::File::create(output)?);
    let encoder = image::png::PngEncoder::new_with_quality(
        w,
        image::png::CompressionType::Best,
        image::png::FilterType::NoFilter,
    );

    // Clean up if file is created but there's a problem writing to it
    match encoder.encode(&imgbuf, width, height, image::ColorType::Rgb8) {
        Ok(_) => {}
        Err(err) => {
            eprintln!("{}", err);
            std::fs::remove_file(output)?;
        }
    }

    Ok(())
}

/// Save palette image file.
pub fn save_palette(
    res: &[[u8; 3]],
    height: u32,
    width: Option<u32>,
    title: &std::path::PathBuf,
) -> Result<(), Box<dyn std::error::Error>> {
    let len = res.len() as u32;
    let w = match width {
        Some(x) => {
            // Width must be at least the length of res
            if x < len {
                len
            } else {
                x
            }
        }
        None => height * len,
    };

    let mut imgbuf: image::RgbImage = image::ImageBuffer::new(w, height);

    for (x, _, pixel) in imgbuf.enumerate_pixels_mut() {
        let color = *res
            .get(
                (((x as f32 / w as f32) * len as f32 - 0.5)
                    .max(0.0)
                    .min(len as f32))
                .round() as usize,
            )
            .unwrap();
        *pixel = image::Rgb(color);
    }

    Ok(save_image(title, &imgbuf.to_vec(), w, height)?)
}

/// Save palette image file.
pub fn save_palette_lab(
    res: &[palette::Srgb<u8>],
    height: u32,
    width: Option<u32>,
    title: &std::path::PathBuf,
) -> Result<(), Box<dyn std::error::Error>> {
    let len = res.len() as u32;
    let w = match width {
        Some(x) => {
            // Width must be at least the length of res
            if x < len {
                len
            } else {
                x
            }
        }
        None => height * len,
    };

    let mut imgbuf: image::RgbImage = image::ImageBuffer::new(w, height);

    for (x, _, pixel) in imgbuf.enumerate_pixels_mut() {
        let color = *res
            .get(
                (((x as f32 / w as f32) * len as f32 - 0.5)
                    .max(0.0)
                    .min(len as f32))
                .round() as usize,
            )
            .unwrap();
        *pixel = image::Rgb([color.red, color.green, color.blue]);
    }

    Ok(save_image(title, &imgbuf.to_vec(), w, height)?)
}
