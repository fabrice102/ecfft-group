// use ark_bn254::Fr;
// use ark_poly::univariate::DensePolynomial;
// use ark_poly::{EvaluationDomain, Polynomial, Radix2EvaluationDomain};
use ark_std::{rand::Rng, test_rng, time::Duration};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use ecfft_bn254::ecfft::EcFftParameters;
use ecfft_bn254::group_polynomial::DenseGroupPolynomial;
use ecfft_bn254::ed25519sc::{Ed25519EcFftParameters};

fn evaluations(c: &mut Criterion) {
    type P = Ed25519EcFftParameters;
    type F = ark_ed25519::Fr;
    type G = ark_ed25519::EdwardsProjective;

    let precomputation = P::precompute();
    let mut rng = test_rng();

    let mut group = c.benchmark_group("ed25519sc");
    group.measurement_time(Duration::from_secs(5));

    for log_n in 1..=P::LOG_N {
        group.bench_with_input(BenchmarkId::new("ECFFT-sc", log_n), &log_n, |b, _| {
            let coeffs: Vec<F> = (0..1 << log_n).map(|_| rng.gen()).collect();
            let poly = DenseGroupPolynomial { coeffs };
            b.iter(|| precomputation.evaluate_over_domain(&poly));
        });
        group.bench_with_input(BenchmarkId::new("Naive-sc", log_n), &log_n, |b, _| {
            let coeffs: Vec<F> = (0..1 << log_n).map(|_| rng.gen()).collect();
            let poly = DenseGroupPolynomial { coeffs };
            let coset = P::sub_coset(P::LOG_N - log_n);
            b.iter(|| coset.iter().map(|x| poly.evaluate(x)).collect::<Vec<_>>());
        });
        group.bench_with_input(BenchmarkId::new("ECFFT-pt", log_n), &log_n, |b, _| {
            let coeffs: Vec<G> = (0..1 << log_n).map(|_| rng.gen()).collect();
            let poly = DenseGroupPolynomial { coeffs };
            b.iter(|| precomputation.evaluate_over_domain(&poly));
        });
        group.bench_with_input(BenchmarkId::new("Naive-pt", log_n), &log_n, |b, _| {
            let coeffs: Vec<G> = (0..1 << log_n).map(|_| rng.gen()).collect();
            let poly = DenseGroupPolynomial { coeffs };
            let coset = P::sub_coset(P::LOG_N - log_n);
            b.iter(|| coset.iter().map(|x| poly.evaluate(x)).collect::<Vec<_>>());
        });
    }
}

criterion_group!(benches, evaluations);
criterion_main!(benches);
