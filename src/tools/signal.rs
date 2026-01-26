//! Signal processing utilities.
//!
//! Provides functions for signal generation and analysis.

use num_complex::Complex64;
use rustfft::{FftPlanner, num_complex::Complex};
use std::f64::consts::PI;

/// Generate a Gaussian pulse.
///
/// # Arguments
/// * `n` - Number of samples
/// * `t0` - Center time index
/// * `sigma` - Standard deviation (in samples)
pub fn gaussian_pulse(n: usize, t0: f64, sigma: f64) -> Vec<f64> {
    (0..n)
        .map(|i| {
            let t = i as f64;
            let arg = (t - t0) / sigma;
            (-0.5 * arg * arg).exp()
        })
        .collect()
}

/// Generate a sinusoidal signal.
///
/// # Arguments
/// * `n` - Number of samples
/// * `freq` - Normalized frequency (cycles per sample)
/// * `phase` - Initial phase (radians)
pub fn sinusoid(n: usize, freq: f64, phase: f64) -> Vec<f64> {
    (0..n)
        .map(|i| (2.0 * PI * freq * i as f64 + phase).sin())
        .collect()
}

/// Compute the FFT of a real signal.
pub fn fft(signal: &[f64]) -> Vec<Complex64> {
    let n = signal.len();
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(n);

    let mut buffer: Vec<Complex<f64>> = signal
        .iter()
        .map(|&x| Complex::new(x, 0.0))
        .collect();

    fft.process(&mut buffer);

    buffer
        .into_iter()
        .map(|c| Complex64::new(c.re, c.im))
        .collect()
}

/// Compute the inverse FFT.
pub fn ifft(spectrum: &[Complex64]) -> Vec<f64> {
    let n = spectrum.len();
    let mut planner = FftPlanner::new();
    let ifft = planner.plan_fft_inverse(n);

    let mut buffer: Vec<Complex<f64>> = spectrum
        .iter()
        .map(|c| Complex::new(c.re, c.im))
        .collect();

    ifft.process(&mut buffer);

    // Normalize
    let scale = 1.0 / n as f64;
    buffer.iter().map(|c| c.re * scale).collect()
}

/// Compute the power spectrum (magnitude squared).
pub fn power_spectrum(signal: &[f64]) -> Vec<f64> {
    let spectrum = fft(signal);
    spectrum.iter().map(|c| c.norm_sqr()).collect()
}

/// Find the peak frequency in a signal.
///
/// Returns the normalized frequency (0 to 0.5) of the dominant frequency component.
pub fn peak_frequency(signal: &[f64]) -> f64 {
    let ps = power_spectrum(signal);
    let n = ps.len();
    let half_n = n / 2;

    // Find peak in positive frequencies
    let (peak_idx, _) = ps[1..half_n]
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .unwrap();

    (peak_idx + 1) as f64 / n as f64
}

/// Compute the analytic signal using Hilbert transform.
pub fn hilbert(signal: &[f64]) -> Vec<Complex64> {
    let n = signal.len();
    let mut spectrum = fft(signal);

    // Apply Hilbert transform in frequency domain
    // H(f) = -j*sign(f)
    spectrum[0] = Complex64::new(spectrum[0].re, 0.0); // DC component

    let half_n = n / 2;
    for elem in spectrum.iter_mut().take(half_n).skip(1) {
        *elem *= Complex64::new(2.0, 0.0);
    }
    if n.is_multiple_of(2) {
        spectrum[half_n] = Complex64::new(spectrum[half_n].re, 0.0);
    }
    for elem in spectrum.iter_mut().take(n).skip(half_n + 1) {
        *elem = Complex64::new(0.0, 0.0);
    }

    // Inverse FFT
    let mut planner = FftPlanner::new();
    let ifft = planner.plan_fft_inverse(n);

    let mut buffer: Vec<Complex<f64>> = spectrum
        .iter()
        .map(|c| Complex::new(c.re, c.im))
        .collect();

    ifft.process(&mut buffer);

    let scale = 1.0 / n as f64;
    buffer
        .into_iter()
        .map(|c| Complex64::new(c.re * scale, c.im * scale))
        .collect()
}

/// Compute the envelope of a signal.
pub fn envelope(signal: &[f64]) -> Vec<f64> {
    hilbert(signal).iter().map(|c| c.norm()).collect()
}

/// Interpolate a signal to a new sample rate.
pub fn interpolate(signal: &[f64], new_length: usize) -> Vec<f64> {
    let old_length = signal.len();
    if new_length == old_length {
        return signal.to_vec();
    }

    let ratio = (old_length - 1) as f64 / (new_length - 1) as f64;

    (0..new_length)
        .map(|i| {
            let pos = i as f64 * ratio;
            let idx = pos.floor() as usize;
            let frac = pos - idx as f64;

            if idx + 1 >= old_length {
                signal[old_length - 1]
            } else {
                signal[idx] * (1.0 - frac) + signal[idx + 1] * frac
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gaussian_pulse() {
        let pulse = gaussian_pulse(100, 50.0, 10.0);
        assert!((pulse[50] - 1.0).abs() < 1e-10); // Peak at center
        assert!(pulse[0] < 0.01); // Near zero at edges
    }

    #[test]
    fn test_fft_ifft_roundtrip() {
        let signal: Vec<f64> = (0..64).map(|i| (i as f64 * 0.1).sin()).collect();
        let spectrum = fft(&signal);
        let recovered = ifft(&spectrum);

        for (a, b) in signal.iter().zip(recovered.iter()) {
            assert!((a - b).abs() < 1e-10);
        }
    }

    #[test]
    fn test_peak_frequency() {
        let freq = 0.1; // 10% of Nyquist
        let signal = sinusoid(1000, freq, 0.0);
        let detected = peak_frequency(&signal);
        assert!((detected - freq).abs() < 0.01);
    }
}
