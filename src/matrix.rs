use crate::error::QuantError;

/// Two-dimensional matrix.
#[derive(Clone, Debug)]
pub struct Matrix2d<T> {
    data: Vec<T>,
    width: usize,
    height: usize,
}

impl<T: Clone + Default> Matrix2d<T> {
    /// Creates a new `Matrix2d` initialized with default values.
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            data: vec![T::default(); width * height],
            width,
            height,
        }
    }

    /// Creates a `Matrix2d` from a Vec.
    ///
    /// Panics if `width` by `height` does not equal the supplied Vec's length.
    ///
    /// ```should_panic
    /// let matrix = rscolorq::Matrix2d::from_vec(vec![0, 1], 1, 3);
    /// ```
    pub fn from_vec(data: Vec<T>, width: usize, height: usize) -> Self {
        assert!(width * height == data.len());
        Self {
            data,
            width,
            height,
        }
    }

    /// Returns the width of the `Matrix2d`.
    pub fn width(&self) -> usize {
        self.width
    }

    /// Returns the height of the `Matrix2d`.
    pub fn height(&self) -> usize {
        self.height
    }

    /// Returns a reference to the indexed element.
    pub fn get(&self, i: usize, j: usize) -> Option<&T> {
        self.data.get(j * self.width + i)
    }

    /// Returns a mutable reference to the indexed element.
    pub fn get_mut(&mut self, i: usize, j: usize) -> Option<&mut T> {
        self.data.get_mut(j * self.width + i)
    }

    /// Consume the `Matrix2d` and return the underlying Vec.
    pub fn into_raw_vec(self) -> Vec<T> {
        self.data
    }

    /// Returns an iterator to the underlying Vec.
    pub fn iter(&self) -> core::slice::Iter<'_, T> {
        self.data.iter()
    }

    /// Returns a mutable iterator to the underlying Vec.
    pub fn iter_mut(&mut self) -> core::slice::IterMut<'_, T> {
        self.data.iter_mut()
    }
}

impl<T: Clone + Default + Copy + core::ops::MulAssign<T>> Matrix2d<T> {
    /// Multiply a row in the matrix by a scalar value.
    pub fn multiply_row_by_scalar(&mut self, row: usize, factor: T) -> Result<(), QuantError> {
        for i in 0..self.width {
            *self
                .get_mut(i, row)
                .ok_or_else(|| "Could not multiply row by scalar")? *= factor;
        }

        Ok(())
    }
}

impl<T> Matrix2d<T>
where
    T: Clone + Default + Copy + core::ops::AddAssign<T> + core::ops::Mul<T, Output = T>,
{
    /// Add a multiple of one row to another.
    pub fn add_row_multiple(
        &mut self,
        from_row: usize,
        to_row: usize,
        factor: T,
    ) -> Result<(), QuantError> {
        for i in 0..self.width {
            let from_item = *self
                .get(i, from_row)
                .ok_or_else(|| "Index out of range in add_row_multiply")?
                * factor;
            *self
                .get_mut(i, to_row)
                .ok_or_else(|| "Index out of range in add_row_multiply")? += from_item;
        }

        Ok(())
    }
}

impl<T> Matrix2d<T>
where
    T: Clone
        + Default
        + Copy
        + core::ops::AddAssign<T>
        + core::ops::Mul<T, Output = T>
        + core::ops::MulAssign<T>
        + core::ops::Neg<Output = T>
        + crate::MatrixComponent,
{
    /// Returns the inverse of the matrix.
    pub fn matrix_inverse(&self) -> Result<Self, QuantError> {
        // Gaussian elimination, matrices are K x K where K is size of palette
        let mut a = self.clone();

        // Create identity matrix
        let mut result = Self::new(a.width, a.height);
        for i in 0..a.width {
            if let Some(identity) = result.get_mut(i, i) {
                *identity = T::identity();
            }
        }

        // Reduce to echelon form, mirroring in result
        for i in 0..a.width {
            result.multiply_row_by_scalar(
                i,
                T::inverse(a.get(i, i).ok_or_else(|| "Could not reduce matrix")?),
            )?;
            a.multiply_row_by_scalar(
                i,
                T::inverse(a.get(i, i).ok_or_else(|| "Could not reduce matrix")?),
            )?;
            for j in (i + 1)..a.height {
                result.add_row_multiple(
                    i,
                    j,
                    -*a.get(i, j)
                        .ok_or_else(|| "Could not reduce matrix, inner loop")?,
                )?;
                a.add_row_multiple(
                    i,
                    j,
                    -*a.get(i, j)
                        .ok_or_else(|| "Could not reduce matrix, inner loop")?,
                )?;
            }
        }

        // Back substitute
        for i in (0..a.width).rev() {
            let mut j = i as isize - 1;
            while j >= 0 {
                result.add_row_multiple(
                    i,
                    j as usize,
                    -*a.get(i, j as usize)
                        .ok_or_else(|| "Could not back substitute")?,
                )?;
                a.add_row_multiple(
                    i,
                    j as usize,
                    -*a.get(i, j as usize)
                        .ok_or_else(|| "Could not back substitute")?,
                )?;
                j -= 1;
            }
        }

        Ok(result)
    }
}

impl<T> core::ops::Mul<Vec<T>> for Matrix2d<T>
where
    T: Clone + Default + Copy + core::ops::AddAssign<T> + core::ops::Mul<T, Output = T>,
{
    type Output = Vec<T>;

    fn mul(self, other: Vec<T>) -> Self::Output {
        assert!(other.len() == self.width());

        let mut result = Vec::with_capacity(self.height());

        for row in 0..self.height() {
            let mut sum = T::default();
            for col in 0..self.width() {
                sum += *self.get(col, row).unwrap() * *other.get(col).unwrap();
            }
            result.push(sum);
        }

        result
    }
}

impl<T> core::ops::Mul<f64> for Matrix2d<T>
where
    T: Clone + Default + Copy + core::ops::Mul<f64, Output = T>,
{
    type Output = Self;

    fn mul(self, other: f64) -> Self {
        Self {
            data: self.iter().map(|&a| a * other).collect(),
            width: self.width,
            height: self.height,
        }
    }
}

/// Three-dimensional matrix.
#[derive(Clone, Debug)]
pub struct Matrix3d<T> {
    data: Vec<T>,
    width: usize,
    height: usize,
    depth: usize,
}

impl<T: Clone + Default> Matrix3d<T> {
    /// Creates a new `Matrix3d` initialized with default values.
    pub fn new(width: usize, height: usize, depth: usize) -> Self {
        Self {
            data: vec![T::default(); width * height * depth],
            width,
            height,
            depth,
        }
    }

    /// Creates a `Matrix3d` from a Vec.
    ///
    /// Panics if `width` by `height` by `depth` does not equal the supplied
    /// Vec's length.
    ///
    /// ```should_panic
    /// let matrix = rscolorq::Matrix3d::from_vec(vec![0, 1, 2], 1, 2, 1);
    /// ```
    pub fn from_vec(data: Vec<T>, width: usize, height: usize, depth: usize) -> Self {
        assert!(width * height * depth == data.len());
        Self {
            data,
            width,
            height,
            depth,
        }
    }

    /// Returns the width of the `Matrix3d`.
    pub fn width(&self) -> usize {
        self.width
    }

    /// Returns the height of the `Matrix3d`.
    pub fn height(&self) -> usize {
        self.height
    }

    /// Returns the depth of the `Matrix3d`.
    pub fn depth(&self) -> usize {
        self.depth
    }

    /// Returns a reference to the indexed element.
    pub fn get(&self, i: usize, j: usize, k: usize) -> Option<&T> {
        self.data
            .get(j * self.width * self.depth + i * self.depth + k)
    }

    /// Returns a mutable reference to the indexed element.
    pub fn get_mut(&mut self, i: usize, j: usize, k: usize) -> Option<&mut T> {
        self.data
            .get_mut(j * self.width * self.depth + i * self.depth + k)
    }

    /// Consume the `Matrix3d` and return the underlying Vec.
    pub fn into_raw_vec(self) -> Vec<T> {
        self.data
    }

    /// Returns an iterator to the underlying Vec.
    pub fn iter(&self) -> core::slice::Iter<'_, T> {
        self.data.iter()
    }

    /// Returns a mutable iterator to the underlying Vec.
    pub fn iter_mut(&mut self) -> core::slice::IterMut<'_, T> {
        self.data.iter_mut()
    }
}

#[cfg(test)]
mod tests {
    use crate::matrix::{Matrix2d, Matrix3d};
    #[test]
    fn array_into_raw() {
        let mat = Matrix2d::<i32>::new(2, 2);
        assert_eq!(mat.into_raw_vec(), vec![0; 4]);
    }

    #[test]
    fn array_indexing_2d() {
        let mat = Matrix2d::from_vec(vec![0, 1, 2, 3], 2, 2);
        assert_eq!(*mat.get(0, 0).unwrap(), 0);
        assert_eq!(*mat.get(1, 0).unwrap(), 1);
        assert_eq!(*mat.get(0, 1).unwrap(), 2);
        assert_eq!(*mat.get(1, 1).unwrap(), 3);
        let mat = Matrix2d::from_vec(vec![0, 1, 2, 3, 4, 5], 2, 3);
        assert_eq!(*mat.get(0, 0).unwrap(), 0);
        assert_eq!(*mat.get(1, 0).unwrap(), 1);
        assert_eq!(*mat.get(0, 1).unwrap(), 2);
        assert_eq!(*mat.get(1, 1).unwrap(), 3);
        assert_eq!(*mat.get(0, 2).unwrap(), 4);
        assert_eq!(*mat.get(1, 2).unwrap(), 5);
        let mat = Matrix2d::from_vec(vec![0, 1, 2, 3, 4, 5], 3, 2);
        assert_eq!(*mat.get(0, 0).unwrap(), 0);
        assert_eq!(*mat.get(1, 0).unwrap(), 1);
        assert_eq!(*mat.get(2, 0).unwrap(), 2);
        assert_eq!(*mat.get(0, 1).unwrap(), 3);
        assert_eq!(*mat.get(1, 1).unwrap(), 4);
        assert_eq!(*mat.get(2, 1).unwrap(), 5);
    }

    #[test]
    fn array_indexing_3d() {
        let mat = Matrix3d::from_vec(vec![0, 1, 2, 3, 4, 5, 6, 7], 2, 2, 2);
        assert_eq!(*mat.get(0, 0, 0).unwrap(), 0);
        assert_eq!(*mat.get(0, 0, 1).unwrap(), 1);
        assert_eq!(*mat.get(1, 0, 0).unwrap(), 2);
        assert_eq!(*mat.get(1, 0, 1).unwrap(), 3);
        assert_eq!(*mat.get(0, 1, 0).unwrap(), 4);
        assert_eq!(*mat.get(0, 1, 1).unwrap(), 5);
        assert_eq!(*mat.get(1, 1, 0).unwrap(), 6);
        assert_eq!(*mat.get(1, 1, 1).unwrap(), 7);
    }

    #[test]
    fn matrix_inverse() {
        let mat = Matrix2d::from_vec(vec![4.0, 7., 2., 6.], 2, 2);
        let inv_result = mat.matrix_inverse().unwrap();
        let inv_expected = Matrix2d::from_vec(vec![0.6, -0.7, -0.2, 0.4], 2, 2);
        for (result, expected) in inv_result.iter().zip(inv_expected.iter()) {
            assert!(f32::abs(result - expected) <= 0.005);
        }
    }
}
