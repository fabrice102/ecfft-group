use ark_ec::{VariableBaseMSM, ScalarMul};
// use ark_bn254::Fr;
// use ark_poly::univariate::DensePolynomial;
// use ark_poly::{EvaluationDomain, Polynomial, Radix2EvaluationDomain};
use ark_std::{rand::Rng, test_rng, time::Duration};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use ecfft_bn254::ecfft::EcFftParameters;
use ecfft_bn254::group_polynomial::{DenseGroupPolynomial, smallest_range};
use ecfft_bn254::ed25519sc::{Ed25519EcFftParameters};

type P = Ed25519EcFftParameters;
type F = ark_ed25519::Fr;
type G = ark_ed25519::EdwardsProjective;

/// Benchmark of polynomial evaluations over ed25519sc
/// and extend
fn poly(c: &mut Criterion) {
    let precomputation = P::precompute();
    let mut rng = test_rng();

    let mut group = c.benchmark_group("ed25519-poly");
    //group.measurement_time(Duration::from_secs(5));
    group.measurement_time(Duration::from_secs(30));

    for log_n in 1..=P::LOG_N {
        // Test name = <coeff_type>-<domain>-<algo>
        // coeff_type =
        //    sc = scalar field = Fr
        //    pt = group element of the elliptic curve ED25519
        // domain =
        //    ecfftDm = ecfft domain points
        //    smallDm = small point from small_range centered on 0
        // algo =
        //    ecfft = ECFFT algorithm
        //    extend = extend part of the ECFFT (only algorithn that is not a poly evaluation)
        //    naive = Horner with no optimization
        //    smallHorner = Hornet optimized for small values (using negation)
        //    multiScalar = use multi-scalar multiplication with precomputation of powers of domain
        // Note that another algorithm consists in using multiscalar multiplication
        // This can be benchmarked using the benchmark below
        // Assuming the power of the evaluation points are pre-computed.
        // TODO: The ECFFT algorithm does not use double-scalar multiplication which could make operations much faster
        group.bench_with_input(BenchmarkId::new("sc-ecfftDm-ecfft", log_n), &log_n, |b, _| {
            let coeffs: Vec<F> = (0..1 << log_n).map(|_| rng.gen()).collect();
            let poly = DenseGroupPolynomial { coeffs };
            b.iter(|| precomputation.evaluate_over_domain(&poly));
        });
        group.bench_with_input(BenchmarkId::new("sc-ecfftDm-naive", log_n), &log_n, |b, _| {
            let coeffs: Vec<F> = (0..1 << log_n).map(|_| rng.gen()).collect();
            let poly = DenseGroupPolynomial { coeffs };
            let coset = P::sub_coset(P::LOG_N - log_n);
            b.iter(|| coset.iter().map(|x| poly.evaluate(x)).collect::<Vec<_>>());
        });
        if log_n < P::LOG_N {
            group.bench_with_input(BenchmarkId::new("pt-ecfftDm-extend", log_n), &log_n, |b, _| {
                let evals: Vec<G> = (0..1 << log_n).map(|_| rng.gen()).collect();
                let precomp = &precomputation.coset_precomputations[P::LOG_N - log_n - 1];
                b.iter(|| precomp.extend(&evals[..]));
            });
        }
        group.bench_with_input(BenchmarkId::new("pt-ecfftDm-ecfft", log_n), &log_n, |b, _| {
            let coeffs: Vec<G> = (0..1 << log_n).map(|_| rng.gen()).collect();
            let poly = DenseGroupPolynomial { coeffs };
            b.iter(|| precomputation.evaluate_over_domain(&poly));
        });
        group.bench_with_input(BenchmarkId::new("pt-ecfftDm-naive", log_n), &log_n, |b, _| {
            let coeffs: Vec<G> = (0..1 << log_n).map(|_| rng.gen()).collect();
            let poly = DenseGroupPolynomial { coeffs };
            let coset = P::sub_coset(P::LOG_N - log_n);
            b.iter(|| coset.iter().map(|x| poly.evaluate(x)).collect::<Vec<_>>());
        });
        group.bench_with_input(BenchmarkId::new("pt-smallDm-hornerSmall", log_n), &log_n, |b, _| {
            let coeffs: Vec<G> = (0..1 << log_n).map(|_| rng.gen()).collect();
            let poly = DenseGroupPolynomial { coeffs };
            let coset = smallest_range(1 << log_n);
            b.iter(|| coset.iter().map(|x| poly.evaluate_small(*x)).collect::<Vec<_>>());
        });
    }
}

/// Benchmark of base operations over ed25519
fn base(c: &mut Criterion) {
    let mut rng = test_rng();

    let mut group = c.benchmark_group("ed25519-base");
    group.measurement_time(Duration::from_secs(5));

    // Adding two random points
    group.bench_function("add-rng", |b| {
        // Multiply a random point by a random scalar
        let pt1: G = rng.gen();
        let pt2: G = rng.gen();
        b.iter(|| pt1 + pt2);
    });

    // Negating a random point
    group.bench_function("neg-rng", |b| {
        // Multiply a random point by a random scalar
        let pt1: G = rng.gen();
        b.iter(|| -pt1);
    });

    for log_n in 0..10 {
        let n = 1 << log_n;
        group.bench_with_input(BenchmarkId::new("convert-to-mul-base", log_n), &log_n, |b, _| {
            let pts: Vec<G> = (0..n).map(|_| rng.gen()).collect();
            b.iter(|| G::batch_convert_to_mul_base(&pts));
        });
        group.bench_with_input(BenchmarkId::new("multi-scalar-mult-rng", log_n), &log_n, |b, _| {
            let pts: Vec<G> = (0..n).map(|_| rng.gen()).collect();
            let pts = G::batch_convert_to_mul_base(&pts);
            let scs: Vec<F> = (0..n).map(|_| rng.gen()).collect();
            b.iter(|| G::msm(&pts, &scs).unwrap());
        });
    }

    group.bench_function("scalar-mult-rng", |b| {
        // Multiply a random point by a random scalar
        let pt: G = rng.gen();
        let sc: F = rng.gen();
        b.iter(|| pt * sc);
    });

    // Multiply a random point by a small scalar
    for sc_val in -1..=2 {
        group.bench_with_input(BenchmarkId::new("scalar-mult", sc_val), &sc_val, |b, _| {
            // Multiply a random point by 1
            let pt: G = rng.gen();
            let sc: F = F::from(sc_val);
            b.iter(|| pt * sc);
        });
    }
}

criterion_group!(benches, base, poly);
criterion_main!(benches);
