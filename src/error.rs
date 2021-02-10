/// Indicates that an error has occurred in the quantization calculation.
#[derive(Debug)]
pub enum QuantError {
    /// An error containing a `String` description of the invalid parameter.
    Parameter(String),
    /// An error containing a `String` description of the failed operation.
    Quantization(String),
}

impl std::fmt::Display for QuantError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            QuantError::Parameter(ref err) | QuantError::Quantization(ref err) => {
                write!(f, "{}", err)
            }
        }
    }
}

impl std::error::Error for QuantError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            QuantError::Parameter(_) | QuantError::Quantization(_) => None,
        }
    }
}

impl std::convert::From<&str> for QuantError {
    fn from(error: &str) -> Self {
        QuantError::Quantization(error.to_owned())
    }
}
