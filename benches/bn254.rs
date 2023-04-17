// use ark_bn254::Fr;
// use ark_poly::univariate::DensePolynomial;
// use ark_poly::{EvaluationDomain, Polynomial, Radix2EvaluationDomain};
use ark_std::{rand::Rng, test_rng, time::Duration};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use ecfft_group::bn254::{Bn254EcFftParameters, F};
use ecfft_group::ecfft::EcFftParameters;
use ecfft_group::group_polynomial::DenseGroupPolynomial;

fn evaluations(c: &mut Criterion) {
    type P = Bn254EcFftParameters;
    let precomputation = P::precompute();
    let mut rng = test_rng();

    let mut group = c.benchmark_group("bn254");
    group.measurement_time(Duration::from_secs(30));

    for log_n in 1..=P::LOG_N {
        group.bench_with_input(BenchmarkId::new("ECFFT", log_n), &log_n, |b, _| {
            let coeffs: Vec<F> = (0..1 << log_n).map(|_| rng.gen()).collect();
            let poly = DenseGroupPolynomial { coeffs };
            b.iter(|| precomputation.evaluate_over_domain(&poly));
        });
        group.bench_with_input(BenchmarkId::new("Naive", log_n), &log_n, |b, _| {
            let coeffs: Vec<F> = (0..1 << log_n).map(|_| rng.gen()).collect();
            let poly = DenseGroupPolynomial { coeffs };
            let coset = P::sub_coset(P::LOG_N - log_n);
            b.iter(|| coset.iter().map(|x| poly.evaluate(x)).collect::<Vec<_>>());
        });
        // group.bench_with_input(BenchmarkId::new("Classic", log_n), &log_n, |b, _| {
        //     let coeffs: Vec<Fr> = (0..1 << log_n).map(|_| rng.gen()).collect();
        //     let poly = DenseGroupPolynomial { coeffs };
        //     let domain = Radix2EvaluationDomain::new(1 << log_n).unwrap();
        //     b.iter(|| poly.clone().evaluate_over_domain(domain));
        // });
    }
}

criterion_group!(benches, evaluations);
criterion_main!(benches);
