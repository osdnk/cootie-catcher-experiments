use std::cmp::max;
use std::f32::RADIX;
use std::iter::{Sum, zip};
use rayon::prelude::*;
use std::ops::{Mul, Add, Sub, Rem};
use std::process::Command;
use std::sync::Mutex;
use num_traits::{One, Zero};
use rand::Rng;
use rayon::{spawn, ThreadPoolBuilder};
use rug::{Complete, Integer};
use crate::cpp::bindings;
use crate::cpp::bindings::{eltwise_add_mod, eltwise_reduce_mod};
use crate::prime_ring::ring::{DPrimeRingElement, poly_mul_mod, PrimeRing, PrimeRingElement};
use crate::prime_ring::r#static::{MAX_THREADS, MOD_Q, PHI};

#[derive(Clone, Debug, PartialEq)]

pub struct PowerSeries {
    pub expanded_layers: Vec<Vec<DPrimeRingElement>>,
    pub tensors: Vec<Vec<DPrimeRingElement>>
}


// Helper function to map a vector into PrimeRingElement::constant
pub fn map_vector_to_prime_ring(vector: Vec<u64>) -> Vec<DPrimeRingElement> {
    vector.into_iter().map(|v| PrimeRing::constant(v)).collect()
}

// Helper function to map a matrix into PrimeRingElement::constant
pub fn map_matrix_to_prime_ring(matrix: Vec<Vec<u64>>) -> Vec<Vec<DPrimeRingElement>> {
    matrix.into_iter().map(|row| map_vector_to_prime_ring(row)).collect()
}

/// Computes the dot product of a vector and a matrix in parallel.
///
/// # Arguments
///
/// * `matrix` - A reference to a matrix represented as a slice of vectors.
/// * `vector` - A reference to a vector.
///
/// # Returns
///
/// A new vector containing the result of the vector-matrix multiplication.
///
/// # Type Parameters
///
/// * `T` - The type of the elements in the matrix and the vector. It must implement the `Mul`, `Zero`, `Copy`, `Send`, `Sync`,
/// and `Add` traits.
///
/// # Panics
///
/// This function will panic if the number of rows in the matrix does not match the length of the vector.
///
/// # Examples
///
/// ```
/// # fn main() {
/// let matrix = vec![
///     vec![1, 2],
///     vec![3, 4],
///     vec![5, 6],
/// ];
/// let vector = vec![7, 8, 9];
/// let result = parallel_dot_vector_matrix(&matrix, &vector);
/// assert_eq!(result, vec![76, 100]);
/// # }
/// ```
pub fn parallel_dot_vector_matrix<T>(vector: &[T], matrix: &[Vec<T>]) -> Vec<T>
where
    T: Mul<Output = T> + Zero + Copy + Send + Sync + Add<Output = T>,
{
    assert!(
        matrix.len() == vector.len(),
        "Number of rows in the matrix must match the length of the vector"
    );

    (0..matrix[0].len())
        .into_par_iter()
        .map(|col| {
            matrix
                .iter()
                .zip(vector.iter())
                .map(|(row, &v)| row[col] * v)
                .fold(T::zero(), |acc, x| acc + x)
        })
        .collect()
}

/// Computes the dot product of a matrix and a vector in parallel.
///
/// # Arguments
///
/// * `matrix` - A reference to a matrix represented as a slice of vectors (each inner vector is a row).
/// * `vector` - A reference to a vector.
///
/// # Returns
///
/// A new vector containing the result of the matrix-vector multiplication.
///
/// # Type Parameters
///
/// * `T` - The type of the elements in the matrix and the vector. It must implement the `Mul`, `Zero`, `Copy`, `Send`, `Sync`,
/// and `Add` traits.
///
/// # Panics
///
/// This function will panic if the number of columns in the matrix does not match the length of the vector.
///
/// # Examples
///
/// ```
/// # fn main() {
/// let matrix = vec![
///     vec![1, 2, 3],
///     vec![4, 5, 6],
/// ];
/// let vector = vec![7, 8, 9];
/// let result = parallel_dot_matrix_vector(&matrix, &vector);
/// assert_eq!(result, vec![50, 122]);
/// # }
/// ```
pub fn parallel_dot_matrix_vector<T>(matrix: &[Vec<T>], vector: &[T]) -> Vec<T>
where
    T: Mul<Output = T> + Zero + Copy + Send + Sync + Add<Output = T>,
{
    assert!(
        matrix[0].len() == vector.len(),
        "Number of columns in the matrix must match the length of the vector"
    );

    matrix
        .par_iter()
        .map(|row| {
            row.iter()
                .zip(vector.iter())
                .map(|(&a, &b)| a * b)
                .fold(T::zero(), |acc, x| acc + x)
        })
        .collect()
}

#[test]
fn test_parallel_dot_vector_matrix_integers() {
    let matrix = vec![
        vec![1, 2],
        vec![3, 4],
        vec![5, 6],
    ];
    let vector = vec![7, 8, 9];
    let result = parallel_dot_vector_matrix(&vector, &matrix);
    assert_eq!(result, vec![76, 100]);
}

#[test]
fn test_parallel_dot_vector_matrix_rings() {
    let matrix = vec![
        vec![1, 2],
        vec![3, 4],
        vec![5, 6],
    ];
    let vector = vec![7, 8, 9];

    let mapped_matrix: Vec<Vec<DPrimeRingElement>> = map_matrix_to_prime_ring(matrix);
    let mapped_vector: Vec<DPrimeRingElement> = map_vector_to_prime_ring(vector);

    let result = parallel_dot_vector_matrix(&mapped_vector, &mapped_matrix);

    let expected_result = map_vector_to_prime_ring(vec![76, 100]); // Assuming these are the correct results

    assert_eq!(result, expected_result);
}


/// Multiplies each element in the given vector by a given scalar.
///
/// # Arguments
///
/// * `vector` - A reference to a vector of elements to be multiplied.
/// * `ell` - A reference to a scalar value by which each element of the vector will be multiplied.
///
/// # Returns
///
/// A new vector containing the result of element-wise multiplication.
///
/// # Type Parameters
///
/// * `T` - The type of elements in the vector and the scalar. It must implement the `Mul`, `Zero`, `Copy`, `Send`, `Sync`,
/// and `Add` traits.
///
/// # Examples
///
/// ```
/// # fn main() {
/// let vector = vec![1, 2, 3, 4];
/// let scalar = 2;
/// let result = vector_element_product(&vector, &scalar);
/// assert_eq!(result, vec![2, 4, 6, 8]);
/// # }
/// ```
pub fn vector_element_product<T>(vector: &Vec<T>, ell: &T) -> Vec<T>
where
    T: Mul<Output = T> + Zero + Copy + Send + Sync + Add<Output = T>,
{
    vector.par_iter().map(|x| *x * *ell).collect()
}

#[test]
fn test_vector_element_product_integers() {
    let vector = vec![1, 2, 3, 4];
    let scalar = 2;
    let result = vector_element_product(&vector, &scalar);
    assert_eq!(result, vec![2, 4, 6, 8]);
}

#[test]
fn test_vector_element_product_ring() {
    let vector = vec![1, 2, 3, 4];
    let scalar = 2;

    let mapped_vector = map_vector_to_prime_ring(vector);
    let mapped_scalar = PrimeRing::constant(scalar);

    let result = vector_element_product(&mapped_vector, &mapped_scalar);

    let expected_result = map_vector_to_prime_ring(vec![2, 4, 6, 8]);

    assert_eq!(result, expected_result);
}
#[test]
fn test_vector_element_product_floats() {
    let vector = vec![1.0, 2.0, 3.0, 4.0];
    let scalar = 0.5;
    let result = vector_element_product(&vector, &scalar);
    assert_eq!(result, vec![0.5, 1.0, 1.5, 2.0]);
}

#[test]
fn test_vector_element_product_zeros() {
    let vector = vec![0, 0, 0];
    let scalar = 999;
    let result = vector_element_product(&vector, &scalar);
    assert_eq!(result, vec![0, 0, 0]);
}

/// Computes the inner product of two vectors by splitting into M parallel chunks.
///
/// # Arguments
///
/// * `vec_a` - The first input vector.
/// * `vec_b` - The second input vector.
/// * `chunk_size` - The number of parallel chunks to split the computation into.
///
/// # Returns
///
/// * The inner product of the two vectors.
/// Computes the inner product of two vectors by splitting into a specified number of parallel chunks.
///
/// # Arguments
///
/// * `vec_a` - The first input vector.
/// * `vec_b` - The second input vector.
/// * `num_chunks` - The number of parallel chunks to split the computation into.
///
/// # Returns
///
/// * The inner product of the two vectors.
fn par_inner_product<T>(vec_a: &[T], vec_b: &[T], num_chunks: usize) -> T
where
    T: Mul<Output = T> + Zero + Copy + Send + Sync + Add<Output = T> + Sum,
{
    assert_eq!(vec_a.len(), vec_b.len(), "Vectors must be of the same length");

    // Compute the total length of the vectors
    let len = vec_a.len();

    // Calculate the approximate size of each chunk
    let chunk_size = if num_chunks == 0 { len } else { (len + num_chunks - 1) / num_chunks };

    // Compute the inner product in parallel
    (0..len)
        .into_par_iter()
        .with_min_len(chunk_size)
        .map(|i| vec_a[i] * vec_b[i])
        .sum()
}

fn inner_product<T>(vec_a: &[T], vec_b: &[T]) -> T
where
    T: Mul<Output = T> + Zero + Copy + Send + Sync + Add<Output = T> + Sum,
{
    assert_eq!(vec_a.len(), vec_b.len(), "Vectors must be of the same length");

    // Compute the total length of the vectors
    let len = vec_a.len();

    // Calculate the approximate size of each chunk
    (0..len)
        .into_iter()
        .map(|i| vec_a[i] * vec_b[i])
        .sum()
}

/// Computes the dot product of two matrices in parallel.
///
/// # Arguments
///
/// * `matrix_a` - A reference to a slice of vectors representing the first matrix.
/// * `matrix_b` - A reference to a slice of vectors representing the second matrix.
///
/// # Returns
///
/// A matrix (vector of vectors) containing the result of the dot product.
pub fn parallel_dot_matrix_matrix<T>(matrix_a: &[Vec<T>], matrix_b: &[Vec<T>]) -> Vec<Vec<T>>
where
    T: Mul<Output = T> + Zero + Copy + Send + Sync + Add<Output = T> + Sum,
{
    let nrows = matrix_a.len();
    let ncols = matrix_b[0].len();
    let inner_dim = matrix_b.len();

    // Ensure the matrices dimensions are compatible for multiplication
    assert!(matrix_a.first().map_or(true, |row| row.len() == inner_dim));

    // Transpose the second matrix for easier column access
    let matrix_b_t: Vec<Vec<T>> = transpose(&matrix_b.to_vec());

    let ip_threads = MAX_THREADS / (ncols * nrows);

    let res = (0..nrows).into_par_iter().map(|i| {
        (0..ncols).into_par_iter().map(|j| {
            par_inner_product(&matrix_a[i], &matrix_b_t[j], ip_threads)
        }).collect()
    }).collect();

    res
}

pub fn parallel_dot_series_matrix(matrix_a: &[PowerSeries], matrix_b: &[Vec<DPrimeRingElement>]) -> Vec<Vec<DPrimeRingElement>>
{
    let nrows = matrix_a.len();
    let ncols = matrix_b[0].len();
    let inner_dim = matrix_b.len();

    // Transpose the second matrix for easier column access
    let matrix_b_t= transpose(&matrix_b.to_vec());

    let ip_threads = MAX_THREADS / (ncols * nrows);

    let res = (0..nrows).into_par_iter().map(|i| {
        let row_a = matrix_a[i]
            .expanded_layers.iter()
            .filter(|l| l.len() == inner_dim)
            .take(1)
            .next().unwrap();

        (0..ncols).into_par_iter().map(|j| {
            par_inner_product(&row_a, &matrix_b_t[j], ip_threads)
        }).collect()
    }).collect();

    res
}

#[test]
fn test_parallel_dot_matrix_matrix() {
    let matrix_a = vec![
        vec![1, 2, 3],
        vec![4, 5, 6],
        vec![7, 8, 9],
    ];
    let matrix_b = vec![
        vec![1, 4, 7],
        vec![2, 5, 8],
        vec![3, 6, 9],
    ];
    let result = parallel_dot_matrix_matrix(&matrix_a, &matrix_b);
    let expected = vec![
        vec![14, 32, 50],
        vec![32, 77, 122],
        vec![50, 122, 194],
    ];
    assert_eq!(result, expected);
}

#[test]
fn test_parallel_dot_matrix_matrix_ring() {
    let matrix_a = vec![
        vec![1, 2, 3],
        vec![4, 5, 6],
        vec![7, 8, 9],
    ];
    let matrix_b = vec![
        vec![1, 4, 7],
        vec![2, 5, 8],
        vec![3, 6, 9],
    ];

    let mapped_matrix_a = map_matrix_to_prime_ring(matrix_a);
    let mapped_matrix_b = map_matrix_to_prime_ring(matrix_b);

    let result = parallel_dot_matrix_matrix(&mapped_matrix_a, &mapped_matrix_b);

    let expected = vec![
        vec![14, 32, 50],
        vec![32, 77, 122],
        vec![50, 122, 194],
    ];

    let mapped_expected = map_matrix_to_prime_ring(expected);

    assert_eq!(result, mapped_expected);
}

/// Transposes a matrix represented as a vector of vectors.
///
/// # Arguments
///
/// * `matrix` - A vector of vectors representing the matrix to be transposed.
///
/// # Returns
///
/// A new vector of vectors representing the transposed matrix.
///
/// # Example
///
/// ```rust
/// let matrix = vec![
///     vec![1, 2, 3],
///     vec![4, 5, 6],
///     vec![7, 8, 9],
/// ];
/// let result = transpose(matrix);
/// assert_eq!(result, vec![
///     vec![1, 4, 7],
///     vec![2, 5, 8],
///     vec![3, 6, 9],
/// ]);
/// ```
pub fn transpose<T>(matrix: &Vec<Vec<T>>) -> Vec<Vec<T>>
where
    T: Copy + Zero + Send + Sync,
{
    if matrix.is_empty() {
        return vec![];
    }

    let rows = matrix.len();
    let cols = matrix[0].len();
    let mut transposed = vec![vec![T::zero(); rows]; cols];

    // Transpose using parallel iteration
    transposed.par_iter_mut().enumerate().for_each(|(i, row)| {
        for (j, value) in row.iter_mut().enumerate() {
            *value = matrix[j][i];
        }
    });

    transposed
}

#[test]
fn test_transpose_square_matrix() {
    let matrix = vec![
        vec![1, 2, 3],
        vec![4, 5, 6],
        vec![7, 8, 9],
    ];
    let result = transpose(&matrix);
    let expected = vec![
        vec![1, 4, 7],
        vec![2, 5, 8],
        vec![3, 6, 9],
    ];
    assert_eq!(result, expected);
}

/// Extracts the first `n` columns from the input matrix.
///
/// # Arguments
///
/// * `matrix` - A 2D vector representing the input matrix.
/// * `n` - The number of columns to extract from the start.
///
/// # Returns
///
/// A new 2D vector (submatrix) containing the first `n` columns.
///
/// # Panics
///
/// Panics if the matrix is empty or if `n` is greater than the number of columns in the matrix.
///
/// # Examples
///
/// ```
/// let mat = vec![
///     vec![1, 2, 3, 4],
///     vec![5, 6, 7, 8],
///     vec![9, 10, 11, 12],
/// ];
///
/// let sub_mat = first_n_columns(mat, 2);
/// assert_eq!(sub_mat, vec![
///     vec![1, 2],
///     vec![5, 6],
///     vec![9, 10]
/// ]);
/// ```
pub fn first_n_columns<T: Copy>(matrix: &Vec<Vec<T>>, n: usize) -> Vec<Vec<T>> {
    if matrix.is_empty() {
        panic!("Matrix is empty");
    }
    let num_columns = matrix[0].len();
    if n > num_columns {
        panic!("Invalid range: `n` is greater than the number of columns in the matrix");
    }

    let mut submatrix = Vec::new();
    for row in matrix.iter() {
        let sub_row = row[0..n].to_vec();
        submatrix.push(sub_row);
    }

    submatrix
}

pub fn last_n_columns<T: Copy>(mat: &Vec<Vec<T>>, n: usize) -> Vec<Vec<T>> {
    if mat.is_empty() {
        panic!("Matrix is empty");
    }
    let num_columns = mat[0].len();
    if n > num_columns {
        panic!("Invalid range: `n` is greater than the number of columns in the matrix");
    }

    let from = num_columns - n;
    let to = num_columns - 1;

    let mut submatrix = Vec::new();
    for row in mat.iter() {
        let sub_row = row[from..=to].to_vec();
        submatrix.push(sub_row);
    }

    submatrix
}

pub fn random(len: usize, mod_q: u64) -> Vec<u64> {
    let mut rng = rand::thread_rng();
    let mut result = Vec::with_capacity(len);
    for i in 0..len {
        let number = rng.gen_range(0..mod_q);
        result.push(number);
    }
    result
}


/// Samples a random matrix of size n x m where each element is a random RingElement.
///
/// # Arguments
///
/// * `n` - The number of rows in the matrix.
/// * `m` - The number of columns in the matrix.
///
/// # Returns
///
/// A vector of vectors (matrix) where each element is a random RingElement.
pub fn sample_random_mat(n: usize, m: usize) -> Vec<Vec<DPrimeRingElement>> {
    // Create a matrix of size n x m
    (0..n).map(|_| {
        (0..m).map(|_| PrimeRing::random()).collect()
    }).collect()
}

pub fn sample_short_random_mat(n: usize, m: usize) -> Vec<Vec<DPrimeRingElement>> {
    // Create a matrix of size n x m
    (0..n).map(|_| {
        (0..m).map(|_| PrimeRing::random_short()).collect()
    }).collect()
}

pub fn sample_bin_random_mat(n: usize, m: usize) -> Vec<Vec<DPrimeRingElement>> {
    // Create a matrix of size n x m
    (0..n).map(|_| {
        (0..m).map(|_| PrimeRing::random_bin()).collect()
    }).collect()
}


/// Samples a random vector of the given size where each element is a random RingElement.
///
/// # Arguments
///
/// * `size` - The number of elements in the vector.
///
/// # Returns
///
/// A vector where each element is a random RingElement.
pub fn sample_random_vector(size: usize) -> Vec<DPrimeRingElement> {
    // Create a vector of the given size with random RingElement values
    (0..size).map(|_| PrimeRing::random()).collect()
}

pub fn sample_short_random_vector(size: usize) -> Vec<DPrimeRingElement> {
    // Create a vector of the given size with random RingElement values
    (0..size).map(|_| PrimeRing::random_short()).collect()
}

/// Computes the row-wise tensor product of two matrices `a` and `b`.
///
/// # Arguments
///
/// * `a` - A matrix represented as a vector of vectors of generic type `T`.
/// * `b` - A matrix represented as a vector of vectors of generic type `T`.
///
/// # Returns
///
/// A matrix represented as a vector of vectors of generic type `T`, where each row of the result
/// is the tensor product of the corresponding rows of `a` and `b`.
///
/// # Panics
///
/// Panics if the number of rows in `a` and `b` do not match.
///
/// # Examples
///
/// ```
/// let a = vec![
///     vec![1, 2],
///     vec![3, 4],
/// ];
///
/// let b = vec![
///     vec![5, 6],
///     vec![7, 8],
/// ];
///
/// let result = row_wise_tensor(a, b);
/// assert_eq!(result, vec![
///     vec![5, 6, 10, 12],
///     vec![21, 24, 28, 32],
/// ]);
/// ```
pub fn row_wise_tensor<T>(a: &Vec<Vec<T>>, b: &Vec<Vec<T>>) -> Vec<Vec<T>>
where
    T: Copy + Mul<Output = T> + Send + Sync,
{
    if a.len() != b.len() {
        panic!("Matrices `a` and `b` must have the same number of rows");
    }

    // Use the same length from vector `a` as they must be equal
    let num_rows = a.len();

    // Preallocate the result vector with the required size
    let mut result = Vec::with_capacity(num_rows);

    // Use par_iter to parallelize the outer loop
    result.par_extend(
        a.par_iter()
            .zip(b.par_iter())
            .map(|(row_a, row_b)| {
                let mut tensor_row = Vec::with_capacity(row_a.len() * row_b.len());
                for &elem_a in row_a {
                    for &elem_b in row_b {
                        tensor_row.push(elem_a * elem_b);
                    }
                }
                tensor_row
            })
    );

    result
}

pub fn kronecker_product<T>(a: &Vec<T>, b: &Vec<T>) -> Vec<T>
where
    T: Copy + Mul<Output = T> + Send + Sync,
{
    // Use par_iter to parallelize the outer loop
    a.par_iter()
        .flat_map(|&elem_a| {
            b.par_iter()
                .map(move |&elem_b| elem_a * elem_b)
                .collect::<Vec<T>>() // Collect the results into a temporary vector
        })
    .collect() // Collect all results into the final result vector
}
#[test]
fn test_kronecker_product_basic() {
    let a = vec![1, 2, 3];
    let b = vec![4, 5];
    let expected = vec![4, 5, 8, 10, 12, 15];
    let result = kronecker_product(&a, &b);
    assert_eq!(result, expected);
}
#[test]
fn test_row_wise_tensor_normal_case() {
    let a = vec![
        vec![1, 2],
        vec![3, 4],
    ];

    let b = vec![
        vec![5, 6],
        vec![7, 8],
    ];

    let result = row_wise_tensor(&a, &b);
    assert_eq!(result, vec![
        vec![5, 6, 10, 12],
        vec![21, 24, 28, 32],
    ]);
}

/// Adds two matrices element-wise.
///
/// # Arguments
///
/// * `matrix_a` - A reference to the first matrix.
/// * `matrix_b` - A reference to the second matrix.
///
/// # Returns
///
/// A new matrix which is the element-wise sum of `matrix_a` and `matrix_b`.
///
/// # Panics
///
/// This function will panic if the dimensions of the two matrices are not the same.
///
/// # Example
///
/// ```rust
/// let matrix_a = vec![
///     vec![1, 2, 3],
///     vec![4, 5, 6],
/// ];
/// let matrix_b = vec![
///     vec![7, 8, 9],
///     vec![10, 11, 12],
/// ];
/// let result = add_matrices(&matrix_a, &matrix_b);
/// assert_eq!(result, vec![
///     vec![8, 10, 12],
///     vec![14, 16, 18],
/// ]);
/// ```
pub fn add_matrices<T>(matrix_a: &[Vec<T>], matrix_b: &[Vec<T>]) -> Vec<Vec<T>>
where
    T: Add<Output = T> + Copy + Zero + Send + Sync,
{
    let nrows = matrix_a.len();
    let ncols = matrix_a[0].len();

    // Ensure the dimensions of the two matrices are the same
    assert_eq!(nrows, matrix_b.len(), "The number of rows in the matrices must be the same");
    assert_eq!(ncols, matrix_b[0].len(), "The number of columns in the matrices must be the same");

    // Allocate the resulting matrix with the required size
    let mut result = vec![vec![T::zero(); ncols]; nrows];

    result.par_iter_mut().enumerate().for_each(|(i, row)| {
        for j in 0..ncols {
            row[j] = matrix_a[i][j] + matrix_b[i][j];
        }
    });

    result
}

/// Adds corresponding elements of two vectors and returns the result as a new vector.
///
/// # Arguments
///
/// * `vector_a` - A reference to the first vector.
/// * `vector_b` - A reference to the second vector.
///
/// # Returns
///
/// A new vector containing the result of element-wise addition.
///
/// # Type Parameters
///
/// * `T` - The type of the elements in the vectors. It must implement the `Add`, `Zero`, `Copy`, `Send`, and `Sync` traits.
///
/// # Panics
///
/// This function will panic if the input vectors are not of the same length.
///
/// # Examples
///
/// ```
/// # fn main() {
/// let vector_a = vec![1, 2, 3, 4];
/// let vector_b = vec![5, 6, 7, 8];
/// let result = add_vectors(&vector_a, &vector_b);
/// assert_eq!(result, vec![6, 8, 10, 12]);
/// # }
/// ```
pub fn add_vectors<T>(vector_a: &Vec<T>, vector_b: &Vec<T>) -> Vec<T>
where
    T: Add<Output = T> + Copy + Zero + Send + Sync,
{
    assert_eq!(vector_a.len(), vector_b.len(), "Vectors must be of the same length");

    vector_a.iter().zip(vector_b).map(|(&a, &b)| a + b).collect()
}

#[test]
fn test_add_vectors_integers() {
    let vector_a = vec![1, 2, 3, 4];
    let vector_b = vec![5, 6, 7, 8];
    let result = add_vectors(&vector_a, &vector_b);
    assert_eq!(result, vec![6, 8, 10, 12]);
}

#[test]
fn test_add_vectors_floats() {
    let vector_a = vec![1.0, 2.0, 3.0, 4.0];
    let vector_b = vec![0.5, 1.5, 2.5, 3.5];
    let result = add_vectors(&vector_a, &vector_b);
    assert_eq!(result, vec![1.5, 3.5, 5.5, 7.5]);
}

#[test]
#[should_panic(expected = "Vectors must be of the same length")]
fn test_add_vectors_different_lengths() {
    let vector_a = vec![1, 2, 3];
    let vector_b = vec![1, 2, 3, 4];
    add_vectors(&vector_a, &vector_b); // This should panic
}

#[test]
fn test_add_vectors_zeros() {
    let vector_a = vec![0, 0, 0];
    let vector_b = vec![1, 2, 3];
    let result = add_vectors(&vector_a, &vector_b);
    assert_eq!(result, vec![1, 2, 3]);
}

#[test]
fn test_add_matrices_basic() {
    let matrix_a = vec![
        vec![1, 2, 3],
        vec![4, 5, 6],
    ];
    let matrix_b = vec![
        vec![7, 8, 9],
        vec![10, 11, 12],
    ];
    let result = add_matrices(&matrix_a, &matrix_b);
    let expected = vec![
        vec![8, 10, 12],
        vec![14, 16, 18],
    ];
    assert_eq!(result, expected);
}

pub fn extract_matrix_from_power_series(
    series: &Vec<PowerSeries>,
    size: usize
) -> Vec<Vec<DPrimeRingElement>> {
    series.iter().filter_map(|s| {
        s.expanded_layers.iter()
            .filter(|l| l.len() == size)
            .take(1)
            .cloned()
            .next()
    }).collect()
}

// Computes `a` raised to the power of `pow` using exponentiation by squaring.
///
/// This method works for `pow` being non-negative. If `pow` is 0, the result is 1.
///
/// # Arguments
///
/// * `a` - The base value of type `T`.
/// * `pow` - The exponent value of type `u32`.
///
/// # Returns
///
/// A value of type `T` representing `a` raised to the power of `pow`.
///
/// # Example
///
/// ```
/// let result = fast_power(2, 10); // result should be 1024
/// ```
pub fn fast_power<T>(a: T, pow: u32) -> T
where
    T: Mul<Output = T> + Copy + One,
{
    if pow == 0 {
        return T::one();
    }

    let mut base = a;
    let mut exponent = pow;
    let mut result = T::one();

    while exponent > 0 {
        if exponent % 2 != 0 {
            result = result * base;
        }
        base = base * base;
        exponent /= 2;
    }

    result
}

#[test]
fn test_fast_power() {
    assert_eq!(fast_power(2, 10), 1024);
    assert_eq!(fast_power(3, 0), 1); // a^0 should be 1
    assert_eq!(fast_power(5, 3), 125); // 5^3 = 5 * 5 * 5
    assert_eq!(fast_power(2, 1), 2); // 2^1 should be 2
}


pub fn sample_random_ss_mat(n: usize, m: usize) -> Vec<Vec<DPrimeRingElement>> {
    // Create a matrix of size n x m
    (0..n).map(|_| {
        (0..m).map(|_| PrimeRing::sample_subtractive()).collect()
    }).collect()
}


/// Joins two matrices horizontally.
///
/// This function takes two matrices and joins them horizontally, i.e., concatenates them column-wise.
///
/// # Arguments
///
/// * `mat1` - A reference to the first matrix.
/// * `mat2` - A reference to the second matrix.
///
/// # Returns
///
/// A single matrix that is the result of concatenating the input matrices column-wise.
///
/// # Panics
///
/// This function will panic if the matrices do not have the same number of rows.
///
/// # Example
///
/// ```rust
/// let mat1 = vec![
///     vec![1, 2],
///     vec![3, 4],
/// ];
/// let mat2 = vec![
///     vec![5, 6],
///     vec![7, 8],
/// ];
/// let result = join_matrices_horizontally(&mat1, &mat2);
/// assert_eq!(result, vec![
///     vec![1, 2, 5, 6],
///     vec![3, 4, 7, 8],
/// ]);
/// ```
pub fn join_matrices_horizontally<T>(mat1: &Vec<Vec<T>>, mat2: &Vec<Vec<T>>) -> Vec<Vec<T>>
where
    T: Clone,
{
    // Ensure both matrices have the same number of rows
    assert_eq!(mat1.len(), mat2.len(), "Both matrices must have the same number of rows");

    // Initialize the result matrix
    let mut result = Vec::with_capacity(mat1.len());

    // Extend each row of the result with rows from both matrices
    for (row1, row2) in mat1.iter().zip(mat2) {
        let mut new_row = row1.clone();
        new_row.extend(row2.clone());
        result.push(new_row);
    }

    result
}

// Unit testing
#[cfg(test)]
mod join_matrices_horizontally_tests {
    use super::*;

    #[test]
    fn test_join_matrices_horizontally_basic() {
        let mat1 = vec![
            vec![1, 2],
            vec![3, 4],
        ];
        let mat2 = vec![
            vec![5, 6],
            vec![7, 8],
        ];
        let result = join_matrices_horizontally(&mat1, &mat2);
        let expected = vec![
            vec![1, 2, 5, 6],
            vec![3, 4, 7, 8],
        ];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_join_matrices_horizontally_single_row() {
        let mat1 = vec![vec![1, 2, 3]];
        let mat2 = vec![vec![4, 5, 6]];
        let result = join_matrices_horizontally(&mat1, &mat2);
        let expected = vec![vec![1, 2, 3, 4, 5, 6]];
        assert_eq!(result, expected);
    }

    #[test]
    #[should_panic(expected = "Both matrices must have the same number of rows")]
    fn test_join_matrices_horizontally_diff_row_count() {
        let mat1 = vec![
            vec![1, 2],
            vec![3, 4],
        ];
        let mat2 = vec![
            vec![5, 6],
        ];
        join_matrices_horizontally(&mat1, &mat2);
    }

    #[test]
    fn test_join_matrices_horizontally_empty() {
        let mat1: Vec<Vec<i32>> = vec![];
        let mat2: Vec<Vec<i32>> = vec![];
        let result = join_matrices_horizontally(&mat1, &mat2);
        let expected: Vec<Vec<i32>> = vec![];
        assert_eq!(result, expected);
    }
}

/// Calls a Sage script to compute the inverse of a polynomial.
///
/// # Arguments
///
/// * `a` - A `RingElement` to find the inverse for.
///
/// # Returns
///
/// * `Option<RingElement>` - The inverse of the given polynomial if it exists, otherwise `None`.
pub fn call_sage_inverse_polynomial<
    const phi: usize,
    const two_phi_minus_one: usize,
    const mod_q: u64
>(a: &PrimeRingElement<phi, two_phi_minus_one, mod_q>) -> PrimeRingElement<phi, two_phi_minus_one, mod_q> {
        // Prepare the polynomial coefficients as a string.
        let coeffs_a_str = format!("{:?}", a.coeffs).replace(" ", "");

        // Call the Sage script with the required arguments.
        let output = Command::new("sage")
            .arg("inverse.sage")
            .arg(&coeffs_a_str)
            .arg(format!("{:?}", phi + 1).replace(" ", ""))
            .arg(mod_q.to_string())
            .output()
            .expect("Failed to execute Sage script");

        // Check if the Sage script execution was successful.
        if !output.status.success() {
            panic!("Error running Sage script: {:?}", output);
        }

        // Process the output from the Sage script.
        let stdout = String::from_utf8_lossy(&output.stdout);
        if stdout.trim() == "None" {
            panic!("Inverse does not exist!");
        }
        let output = stdout.trim();
        let result: Vec<u64> = output
            .trim_matches(&['[', ']'] as &[_])
            .split(',')
            .map(|s| s.trim().parse().expect("Invalid number"))
            .collect();

    PrimeRingElement {
        coeffs: <[u64; phi]>::try_from(result).unwrap()
    }
}

#[test]
pub fn test_inverse_sage() {
    let a = PrimeRing::random();
    let b = call_sage_inverse_polynomial(&a);
    assert_eq!(a * b, PrimeRing::constant(1));
}

/// Computes a power series for the given element, prefixing the series with a `1`.
///
/// # Arguments
///
/// * `element` - A `RingElement` which serves as the base for the power series.
/// * `len` - The length of the series.
///
/// # Returns
///
/// A vector containing a single vector of `RingElement`s, representing a power series prefixed with `1`.
pub fn compute_one_prefixed_power_series(element: &DPrimeRingElement, len: usize) -> PowerSeries {
    let mut series = Vec::with_capacity(len);
    series.push(PrimeRing::constant(1));
    series.push(element.clone());

    let mut power = element.clone();
    for _ in 2..len {
        power = power.clone() * element.clone();
        series.push(power.clone());
    }

    let mut ps = PowerSeries {
        expanded_layers: vec![],
        tensors: vec![],
    };
    let mut current_dim = len;
    while current_dim % 2 == 0 {
        ps.expanded_layers.push(series[0..current_dim].to_vec());
        current_dim /= 2;
        ps.tensors.push(vec![PrimeRingElement::one(), series[current_dim]]);
    }
    ps.expanded_layers.push(series[0..current_dim].to_vec());
    ps
}

pub fn compute_hp_power_series(elements: &Vec<DPrimeRingElement>) -> PowerSeries {
    let mut ps = PowerSeries {
        expanded_layers: vec![],
        tensors: vec![],
    };

    ps.expanded_layers.push(vec![PrimeRing::constant(1)]);

    for t in elements.iter().rev() {
        let l_factor = PrimeRing::constant(1) - t.clone();
        let r_factor = t.clone();
        ps.tensors.push(vec![l_factor.clone(), r_factor.clone()]);
        ps.expanded_layers.push(kronecker_product(ps.tensors.last().unwrap(), ps.expanded_layers.last().unwrap()));
    }
    PowerSeries {
        expanded_layers: ps.expanded_layers.iter().rev().cloned().collect(),
        tensors: ps.tensors.iter().rev().cloned().collect(),
    }
}


#[test]
fn test_compute_one_prefixed_power_series() {
    let element = PrimeRing::constant(2);
    let result = compute_one_prefixed_power_series(&element, 4);
    assert_eq!(
        result,
        PowerSeries {
            expanded_layers: map_matrix_to_prime_ring(vec![
                vec![1, 2, 4, 8],
                vec![1, 2],
                vec![1],
            ]),
            tensors: map_matrix_to_prime_ring(vec![
                vec![1, 4],
                vec![1, 2]
            ]),
        }
    );
}

pub fn ring_inner_product(a: &Vec<DPrimeRingElement<>>, b: &Vec<DPrimeRingElement>) -> DPrimeRingElement {
    assert_eq!(a.len(), b.len(), "Input vectors must have the same length");

    a.par_iter()
        .zip(b.par_iter())
        .map(|(a_i, b_i)| *a_i * *b_i)
        .reduce(DPrimeRingElement::zero, |acc, prod| acc + prod)
}






/// Computes the conjugate of each element in a `RingElement` vector.
///
/// # Arguments
///
/// * `series` - A vector of `RingElement`s.
///
/// # Returns
///
/// A new vector where each `RingElement` in the input has been replaced by its conjugate.
pub fn conjugate_vector(row: &Vec<DPrimeRingElement>) -> Vec<DPrimeRingElement> {
    row.iter().map(DPrimeRingElement::conjugate).collect()
}

// RNS: assume small elements.
fn determine_bit_width(matrix: &Vec<Vec<DPrimeRingElement>>) -> usize {
    let max_value = matrix.par_iter()
        .flat_map_iter(|row| row.iter())
        .flat_map_iter(|cell| cell.coeffs.iter())
        .cloned()
        .max()
        .unwrap_or(0);

    64 - max_value.leading_zeros() as usize
}

pub fn neg_matrix(
    matrix: &Vec<Vec<DPrimeRingElement>>,
) -> Vec<Vec<DPrimeRingElement>> {
    matrix.iter().map(|row| {
        row.iter().map(|el| {
            PrimeRingElement::zero() - *el
        }).collect()
    }).collect()
}
pub fn split_into_positive_and_negative(
    matrix: &Vec<Vec<DPrimeRingElement>>,
) -> (Vec<Vec<DPrimeRingElement>>, Vec<Vec<DPrimeRingElement>>){
    // Initialize two matrices for the results
    let mut positive_matrix: Vec<Vec<DPrimeRingElement>> = Vec::new();
    let mut negative_matrix: Vec<Vec<DPrimeRingElement>> = Vec::new();
    for row in matrix.iter() {
        let mut positive_row: Vec<DPrimeRingElement> = Vec::new();
        let mut negative_row: Vec<DPrimeRingElement> = Vec::new();
        for el in row {
            let mut positive_element = PrimeRingElement::zero();
            let mut negative_element = PrimeRingElement::zero();
            for i in 0..el.coeffs.len() {
                if 2 * el.coeffs[i] > MOD_Q { // assume
                    negative_element.coeffs[i] = MOD_Q - el.coeffs[i]
                } else {
                    positive_element.coeffs[i] = el.coeffs[i]
                }
            }
            positive_row.push(positive_element);
            negative_row.push(negative_element);
        }
        positive_matrix.push(positive_row);
        negative_matrix.push(negative_row);

    }
    (positive_matrix, negative_matrix)
}

pub fn decompose(
    matrix: &Vec<Vec<DPrimeRingElement>>,
    radix: u64,
    num_chunks: usize
) -> Vec<Vec<DPrimeRingElement>> {
        matrix.par_iter()
            .map(|row| {
                let mut new_row = Vec::with_capacity(row.len() * num_chunks);
                for cell in row {
                    let mut chunks: Vec<DPrimeRingElement> = vec![DPrimeRingElement::zero(); num_chunks];
                    for k in 0..PHI {
                        let mut val = cell.coeffs[k];
                        for chunk_idx in 0..num_chunks {
                            chunks[chunk_idx].coeffs[k] = (val & (radix - 1)) as u64;
                            val >>= radix.trailing_zeros();
                        }
                    }
                    new_row.extend_from_slice(&chunks);
                }
                new_row
            })
        .collect::<Vec<Vec<DPrimeRingElement>>>()
}

pub fn decompose_matrix_by_radix(
    matrix: &Vec<Vec<DPrimeRingElement>>,
    radix: u64
) -> (Vec<Vec<DPrimeRingElement>>, usize) {
    assert!(radix.is_power_of_two(), "Radix must be a power of two");


    let (positive_matrix, negative_matrix) = split_into_positive_and_negative(&matrix);

    let bit_width = max(
        determine_bit_width(&positive_matrix),
        determine_bit_width(&negative_matrix)
    );
    let num_chunks = (bit_width + (radix.trailing_zeros() as usize - 1)) / radix.trailing_zeros() as usize;
    let num_chunks = num_chunks as usize;




    let decomposed_matrix = add_matrices(
        &decompose(&positive_matrix, radix, num_chunks),
        &neg_matrix(&decompose(&negative_matrix, radix, num_chunks)),
    );

    (decomposed_matrix, num_chunks)
}

pub fn decompose_by_radix(
    element: &DPrimeRingElement,
    radix: u64
)  -> Vec<DPrimeRingElement> {

    let (decomposed_matrix, _) = decompose_matrix_by_radix(
        &vec![vec![element.clone()]],
        radix
    );
    decomposed_matrix[0].clone()
}


pub fn decompose_by_chunks(
    element: &DPrimeRingElement,
    num_chunks: usize
)  -> Vec<DPrimeRingElement> {

    let (decomposed_matrix, radix) = decompose_matrix_by_chunks(
        &vec![vec![element.clone()]],
        num_chunks
    );
    decomposed_matrix[0].clone()
}

pub fn decompose_matrix_by_chunks(matrix: &Vec<Vec<DPrimeRingElement>>, num_chunks: usize) -> (Vec<Vec<DPrimeRingElement>>, u64) {
    let (positive_matrix, negative_matrix) = split_into_positive_and_negative(&matrix);

    let bit_width = max(
        determine_bit_width(&positive_matrix),
        determine_bit_width(&negative_matrix)
    );

    let radix_log = (bit_width as f64 / num_chunks as f64).ceil() as u32;
    let radix = (2u64).pow(radix_log);

    let decomposed_matrix = add_matrices(
        &decompose(&positive_matrix, radix, num_chunks),
        &neg_matrix(&decompose(&negative_matrix, radix, num_chunks)),
    );

    (decomposed_matrix, radix)
}

#[test]
fn test_determine_bit_width() {
    let matrix = vec![
        vec![PrimeRing::constant(12345), PrimeRing::constant(67890)]
    ];
    assert_eq!(determine_bit_width(&matrix), 17);

    let matrix = vec![
        vec![PrimeRing::constant(1), PrimeRing::constant(2)],
        vec![PrimeRing::constant(3), PrimeRing::constant(4)]
    ];
    assert_eq!(determine_bit_width(&matrix), 3);

    let matrix = vec![
        vec![PrimeRing::constant(u64::MAX), PrimeRing::constant(u64::MAX)]
    ];
    assert_eq!(determine_bit_width(&matrix), 64);
}

#[test]
fn test_decompose_by_radix() {
    let matrix = map_matrix_to_prime_ring(vec![
        vec![17, 19],
        vec![17, 19],
    ]);
    let radix = 16;

    let (decomposed_matrix, num_chunks) = decompose_matrix_by_radix(&matrix, radix);
    assert_eq!(num_chunks, 2);
    assert_eq!(decomposed_matrix, map_matrix_to_prime_ring(vec![
        vec![1, 1, 3, 1],
        vec![1, 1, 3, 1]
    ]));
}

#[test]
fn test_decompose_by_radix_2() {
    let matrix = map_matrix_to_prime_ring(vec![
        vec![15, 15],
        vec![15, 15],
    ]);
    let radix = 16;

    let (decomposed_matrix, num_chunks) = decompose_matrix_by_radix(&matrix, radix);
    assert_eq!(num_chunks, 1);
    assert_eq!(decomposed_matrix, map_matrix_to_prime_ring(vec![
        vec![15, 15],
        vec![15, 15],
    ]));
}

#[test]
fn test_decompose_by_radix_3() {
    let matrix = map_matrix_to_prime_ring(vec![
        vec![15, 15],
        vec![15, 16],
    ]);
    let radix = 16;

    let (decomposed_matrix, num_chunks) = decompose_matrix_by_radix(&matrix, radix);
    assert_eq!(num_chunks, 2);
    assert_eq!(decomposed_matrix, map_matrix_to_prime_ring(vec![
        vec![15, 0, 15, 0],
        vec![15, 0, 0, 1],
    ]));
}

#[test]
fn test_decompose_by_chunks() {
    let matrix = map_matrix_to_prime_ring(vec![
        vec![15, 240],
        vec![255, 16]
    ]);
    let num_chunks = 4;

    let (decomposed_matrix, num_chunks) = decompose_matrix_by_chunks(&matrix, num_chunks);

    assert_eq!(num_chunks, 4);
    assert_eq!(decomposed_matrix, map_matrix_to_prime_ring(vec![
        vec![3, 3, 0, 0, 0, 0, 3, 3],
        vec![3, 3, 3, 3, 0, 0, 1, 0],
    ]));
}

pub fn compose_with_radix(
    decomposed_matrix: &Vec<Vec<DPrimeRingElement>>,
    radix: u64,
    num_chunks: usize // into how many parts split
) -> Vec<Vec<DPrimeRingElement>> {
    let radix_shift = radix.trailing_zeros() as usize;

    let (positive_matrix, negative_matrix) = split_into_positive_and_negative(&decomposed_matrix);

    let compose = |decomposed_matrix: &Vec<Vec<DPrimeRingElement>>| {
        decomposed_matrix
            .iter()
            .map(|chunked_row| {
                let num_elements_in_chunk = chunked_row.len() / num_chunks;
                (0..num_elements_in_chunk).map(|i| {
                    let mut coeffs = [0u64; PHI];
                    for k in 0..PHI {
                        coeffs[k] = (0..num_chunks).fold(0u64, |acc, el_chunk| {
                            acc | ((chunked_row[i * num_chunks + el_chunk].coeffs[k] as u64) << (radix_shift * el_chunk))
                        });
                    }
                    DPrimeRingElement { coeffs }
                }).collect::<Vec<_>>()
            }).collect::<Vec<_>>()
    };

    add_matrices(
        &compose(&positive_matrix),
        &neg_matrix(&compose(&negative_matrix)),
    )
}

pub fn compose_with_radix_mod(
    decomposed_matrix: &Vec<Vec<DPrimeRingElement>>,
    radix: u64,
    num_chunks: usize, // into how many parts split,
) -> Vec<Vec<DPrimeRingElement>> {
    let radix_shift = radix.trailing_zeros() as usize;
    decomposed_matrix.iter()
        .map(|chunked_row| {

            let num_elements_in_chunk = chunked_row.len() / num_chunks;
            (0..num_elements_in_chunk).map(|i| {
                let mut coeffs = [0u64; PHI];
                for k in 0..PHI {
                    coeffs[k] = (0..num_chunks).fold(0u64, |acc, el_chunk| {
                        #[cfg(any(target_arch = "aarch64", not(feature = "use-hardware")))]
                        {
                            return ((acc as u128 + ((chunked_row[i * num_chunks + el_chunk].coeffs[k] as u128) * 2u128.pow((radix_shift * el_chunk) as u32))) % (MOD_Q as u128)) as u64
                        }

                        #[cfg(all(target_arch = "x86_64", feature = "use-hardware"))]
                        unsafe {
                            return bindings::add_mod(
                                acc,
                                bindings::multiply_mod(
                                    chunked_row[i * num_chunks + el_chunk].coeffs[k],
                                    2u64.pow((radix_shift * el_chunk) as u32),
                                    MOD_Q
                                ),
                                MOD_Q
                            )
                        }
                    });
                }
                DPrimeRingElement { coeffs }
            }).collect::<Vec<_>>()
        }).collect::<Vec<_>>()
}

#[test]
fn test_compose_with_radix() {
    let matrix = sample_random_mat(2, 2);
    let radix = 4;
    let (decomposed_matrix, num_chunks) = decompose_matrix_by_radix(&matrix, radix);
    let composed_matrix = compose_with_radix_mod(&decomposed_matrix, radix, num_chunks);

    assert_eq!(composed_matrix, matrix);
}


pub fn reduce_mod_vec(a: &mut [u64], mod_q: u64) {
    #[cfg(all(target_arch = "x86_64", feature = "use-hardware"))]
    unsafe {
        eltwise_reduce_mod(a.as_mut_ptr(), a.as_mut_ptr(), a.len() as u64, mod_q);
    }
    #[cfg(any(target_arch = "aarch64", not(feature = "use-hardware")))]
    {
        for i in 0..a.len() {
            a[i] = a[i] % mod_q;
        }
    }
}

pub fn reduce_mod(a: &mut u64, mod_q: u64) {
    reduce_mod_vec(&mut [*a], mod_q)
}

pub fn add_mod(a: u64, b: u64, mod_q: u64) -> u64 {
    #[cfg(any(target_arch = "aarch64", not(feature = "use-hardware")))]
    {
        return ((a as u128 + b as u128) % (mod_q as u128)) as u64;
    }
    unsafe { bindings::add_mod(a, b, mod_q) }
}

pub fn sub_mod(a: u64, b: u64, mod_q: u64) -> u64 {
    #[cfg(any(target_arch = "aarch64", not(feature = "use-hardware")))]
    {
        return ((mod_q as u128 + a as u128 - b as u128) % (mod_q as u128)) as u64;
    }
    unsafe { bindings::sub_mod(a, b, mod_q) }
}

pub fn multiply_mod(a: u64, b: u64, mod_q: u64) -> u64 {
    #[cfg(any(target_arch = "aarch64", not(feature = "use-hardware")))]
    {
        return ((a as u128 * b as u128) % (mod_q as u128)) as u64;
    };
    unsafe { bindings::multiply_mod(a, b, mod_q) }
}

fn inverse_rns_slow(remainders: &[u64], primes: &[u64]) -> Integer {
    assert_eq!(remainders.len(), primes.len(), "Mismatched remainders and primes length");

    let product = primes.iter()
        .fold(Integer::from(1), |acc, &p| acc * Integer::from(p));

    let mut result = Integer::new();
    for (&r, &p) in remainders.iter().zip(primes.iter()) {
        let p_int = Integer::from(p);
        let n = (&product / &p_int).complete(); // Convert division result to Integer
        let m = n.clone().invert(&p_int).expect("Factors must be pairwise coprime");
        result += Integer::from(r) * &n * &m;
    }

    result % product
}


