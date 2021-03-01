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
        // j * self.width + i
        j.checked_mul(self.width)
            .map_or_else(|| None, |a| a.checked_add(i))
            .map_or_else(|| None, |index| self.data.get(index))
    }

    /// Returns a mutable reference to the indexed element.
    pub fn get_mut(&mut self, i: usize, j: usize) -> Option<&mut T> {
        j.checked_mul(self.width)
            .map_or_else(|| None, |a| a.checked_add(i))
            .map_or_else(|| None, move |index| self.data.get_mut(index))
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

    /// Returns an iterator over the rows of the `Matrix2d`.
    pub fn rows(&self) -> core::slice::ChunksExact<'_, T> {
        self.data.chunks_exact(self.width)
    }

    /// Returns a mutable iterator over the rows of the `Matrix2d`.
    pub fn rows_mut(&mut self) -> core::slice::ChunksExactMut<'_, T> {
        self.data.chunks_exact_mut(self.width)
    }
}

impl<T: Clone + Default + Copy + core::ops::MulAssign<T>> Matrix2d<T> {
    /// Multiply a row in the matrix by a scalar value.
    pub fn multiply_row_by_scalar(&mut self, row: usize, factor: T) -> Result<(), QuantError> {
        self.rows_mut()
            .nth(row)
            .ok_or("Could not multiply row by scalar")?
            .iter_mut()
            .for_each(|i| *i *= factor);

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
        /*
        To iterate over the Vec and mutate it, we have to use split_at_mut.
        Whichever row is smaller will be the last row in the left slice from
        split_at_mut.

        Sample matrix:
        2 2 2
        1 1 1
        0 0 0
        3 3 3
        4 4 4

        Less case (from_row < to_row):
            from 1 to 3 midpoint calculation:
            (from_row + 1) * self.width =>
                (1 + 1) * 3 = 6 // add 1 for start of the next row
            [2 2 2 1 1 1], [0 0 0 3 3 3 4 4 4]

            from_row_chunk_index = last chunk in split_left

            to_row - from_row =>
                3 - 1 = 2 as the new slice chunk index 1 beyond our target
            to_index_chunk_index for split_right =>
                2 - 1 = 1 as the nth chunk index we want for [3 3 3]
            therefore:
                nth(to_row - from_row - 1)
            [[0 0 0] [3 3 3] [4 4 4]]

        Greater case (from_row > to_row):
            from_row = nth(from_row - to_row - 1) chunk in split_right
            to_row = last chunk in split_left
        */

        match from_row.cmp(&to_row) {
            core::cmp::Ordering::Less => {
                let mid = from_row
                    .checked_add(1)
                    .ok_or("Index out of bounds in add add_row_multiple")?
                    .checked_mul(self.width)
                    .ok_or("Index out of bounds in mul add_row_multiple")?;

                let (split_left, split_right) = self.data.split_at_mut(mid);

                let from_row_chunk = split_left
                    .chunks_exact(self.width)
                    .last()
                    .ok_or("Could not split from slice in add_row_multiple")?;

                let to_row_chunk = split_right
                    .chunks_exact_mut(self.width)
                    .nth(to_row - from_row - 1)
                    .ok_or("Could not split to slice in add_row_multiple")?;

                for (to, from) in to_row_chunk.iter_mut().zip(from_row_chunk) {
                    *to += *from * factor;
                }
            }
            core::cmp::Ordering::Greater => {
                let mid = to_row
                    .checked_add(1)
                    .ok_or("Index out of bounds in add add_row_multiple")?
                    .checked_mul(self.width)
                    .ok_or("Index out of bounds in mul add_row_multiple")?;

                let (split_left, split_right) = self.data.split_at_mut(mid);

                let from_row_chunk = split_right
                    .chunks_exact(self.width)
                    .nth(from_row - to_row - 1)
                    .ok_or("Could not split from slice in add_row_multiple")?;

                let to_row_chunk = split_left
                    .chunks_exact_mut(self.width)
                    .last()
                    .ok_or("Could not split to slice in add_row_multiple")?;

                for (to, from) in to_row_chunk.iter_mut().zip(from_row_chunk) {
                    *to += *from * factor;
                }
            }
            core::cmp::Ordering::Equal => return Err("From and to rows cannot be equal".into()),
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
                T::inverse(a.get(i, i).ok_or("Could not reduce matrix")?),
            )?;
            a.multiply_row_by_scalar(i, T::inverse(a.get(i, i).ok_or("Could not reduce matrix")?))?;
            for j in (i + 1)..a.height {
                result.add_row_multiple(
                    i,
                    j,
                    -*a.get(i, j).ok_or("Could not reduce matrix, inner loop")?,
                )?;
                a.add_row_multiple(
                    i,
                    j,
                    -*a.get(i, j).ok_or("Could not reduce matrix, inner loop")?,
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
                    -*a.get(i, j as usize).ok_or("Could not back substitute")?,
                )?;
                a.add_row_multiple(
                    i,
                    j as usize,
                    -*a.get(i, j as usize).ok_or("Could not back substitute")?,
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

        self.rows()
            .map(|row| {
                row.iter()
                    .zip(&other)
                    .fold(T::default(), |mut sum, (&col, &other_entry)| {
                        sum += col * other_entry;
                        sum
                    })
            })
            .collect()
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
        // The base formula for access is:
        //     j * self.width * self.depth + i * self.depth + k
        //     ^---------(j_term)--------^   ^----(i_term)----^
        j.checked_mul(self.width)
            .map_or_else(|| None, |a| a.checked_mul(self.depth))
            .map_or_else(
                || None,
                |a| {
                    i.checked_mul(self.depth)
                        .map_or_else(|| None, |i_term| a.checked_add(i_term))
                },
            )
            .map_or_else(|| None, |a| a.checked_add(k))
            .map_or_else(|| None, |index| self.data.get(index))
    }

    /// Returns a mutable reference to the indexed element.
    pub fn get_mut(&mut self, i: usize, j: usize, k: usize) -> Option<&mut T> {
        j.checked_mul(self.width)
            .map_or_else(|| None, |a| a.checked_mul(self.depth))
            .map_or_else(
                || None,
                |a| {
                    i.checked_mul(self.depth)
                        .map_or_else(|| None, |i_term| a.checked_add(i_term))
                },
            )
            .map_or_else(|| None, |a| a.checked_add(k))
            .map_or_else(|| None, move |index| self.data.get_mut(index))
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
