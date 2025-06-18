extern crate criterion;
extern crate num_modular;

use criterion::{Criterion, criterion_group, criterion_main, BatchSize};
use rand::{thread_rng, Rng};

#[link(name = "hexl_wrapper")]
extern "C" {
    pub fn multiply_poly(result: *mut u64, operand1: *const u64, operand2: *const u64, n: u64, modulus: u64);
}

fn criterion_benchmark(c: &mut Criterion) {
    // let number: u128 = 1_234_567_890_123_456_789_012_345_678;
    // let modulus: u128 = 1_000_000_007;
    let n = 256;
    let mut rng = thread_rng();

    c.bench_function("HEXL with montgomery AVX-512", |b| {
        let modulus = 4546383823830515713; // Example modulus
        b.iter_batched(
            || {
                // Setup: not timed
                let operand1: Vec<u64> = (0..n).map(|_| rng.gen_range(0..modulus)).collect();
                let operand2: Vec<u64> = (0..n).map(|_| rng.gen_range(0..modulus)).collect();
                let result = vec![0u64; n];
                (operand1, operand2, result)
            },
            |(operand1, operand2, mut result)| {
                // Measurement: timed
                unsafe {
                    multiply_poly(result.as_mut_ptr(), operand1.as_ptr(), operand2.as_ptr(), n as u64, modulus);
                }
            },
            BatchSize::SmallInput,
        )
    });


    c.bench_function("HEXL with montgomery AVX-512-IFMA", |b| {
        let modulus = 4503599626682369; // Example modulus
        b.iter_batched(
            || {
                // Setup: not timed
                let operand1: Vec<u64> = (0..n).map(|_| rng.gen_range(0..modulus)).collect();
                let operand2: Vec<u64> = (0..n).map(|_| rng.gen_range(0..modulus)).collect();
                let result = vec![0u64; n];
                (operand1, operand2, result)
            },
            |(operand1, operand2, mut result)| {
                // Measurement: timed
                unsafe {
                    multiply_poly(result.as_mut_ptr(), operand1.as_ptr(), operand2.as_ptr(), n as u64, modulus);
                }
            },
            BatchSize::SmallInput,
        )
    });

    c.bench_function("HEXL with montgomery AVX-512-IFMA small", |b| {
        let modulus = 562949952700417; // Example modulus
        b.iter_batched(
            || {
                // Setup: not timed
                let operand1: Vec<u64> = (0..n).map(|_| rng.gen_range(0..modulus)).collect();
                let operand2: Vec<u64> = (0..n).map(|_| rng.gen_range(0..modulus)).collect();
                let result = vec![0u64; n];
                (operand1, operand2, result)
            },
            |(operand1, operand2, mut result)| {
                // Measurement: timed
                unsafe {
                    multiply_poly(result.as_mut_ptr(), operand1.as_ptr(), operand2.as_ptr(), n as u64, modulus);
                }
            },
            BatchSize::SmallInput,
        )
    });


    //     c.bench_function("Regular Modulus", |b| {
//         b.iter(|| regular_modulus(black_box(number), black_box(modulus)))
//     });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);


