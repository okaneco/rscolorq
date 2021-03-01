use structopt::StructOpt;

#[derive(StructOpt, Debug)]
#[structopt(name = "rscolorq", about = "Spatial color quantizer")]
pub struct Opt {
    /// Input file.
    #[structopt(short, long, parse(from_os_str))]
    pub input: std::path::PathBuf,

    /// Output file.
    #[structopt(short, long, parse(from_os_str))]
    pub output: Option<std::path::PathBuf>,

    /// Number of colors to dither.
    #[structopt(short, long, default_value = "8", required = false)]
    pub n: u8,

    /// Dithering level, should be between 0.5 and 1.5. Higher numbers will have
    /// higher frequency noise.
    #[structopt(short, long, default_value = "0.8", required = false)]
    pub dither: f64,

    /// Calculate the dithering level automatically based on image size.
    #[structopt(long)]
    pub auto: bool,

    /// Disable saving a dithered image.
    #[structopt(long)]
    pub no_file: bool,

    /// Dither the image with supplied colors.
    #[structopt(
        short,
        long,
        min_values = 2,
        max_values = 255,
        value_delimiter = ",",
        required = false
    )]
    pub colors: Vec<String>,

    /// Use `Lab` color space for dithering.
    #[structopt(long)]
    pub lab: bool,

    /// Print the resulting dither colors.
    #[structopt(short, long)]
    pub print: bool,

    /// Random number generator seed.
    #[structopt(short, long)]
    pub seed: Option<u64>,

    /// Height of color palette image. If width is omitted, palette will be
    /// `height * n` pixels wide.
    #[structopt(long, default_value = "40")]
    pub height: u32,

    /// Width of color palette image. Will be at least `k` pixels wide.
    #[structopt(long)]
    pub width: Option<u32>,

    /// Output file.
    #[structopt(long = "op", parse(from_os_str))]
    pub palette_output: Option<std::path::PathBuf>,

    /// Filter size.
    #[structopt(long, default_value = "3", required = false, hidden = true)]
    pub filter: u8,

    /// Iterations per level.
    #[structopt(long, default_value = "3", required = false)]
    pub iters: usize,

    /// Repeated calculations per temperature.
    #[structopt(long, default_value = "1", required = false)]
    pub repeats: usize,

    /// Initial temperature for annealing.
    #[structopt(long, default_value = "1.0", required = false)]
    pub start: f64,

    /// Final temperature for annealing.
    #[structopt(long, default_value = "0.001", required = false)]
    pub end: f64,
}
