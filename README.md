# swift-fourier

Use this package when you need the frequencies a stretch of samples is made of, or the samples back from their frequencies: spectra, convolution, filtering in frequency. Start from `Fourier.Transform<Scalar>`, planned once for a length with `init(count:)` and then used on any number of signals of that length. `forward(_:)` and `inverse(_:)` take and give arrays of swift-complex's `Complex.Number<Scalar>`; `forward(real:)` and `inverse(real:)` take a real signal to its ⌊n/2⌋ + 1 bins from zero up to half the rate, and back. For the power at one frequency, on a bin or between bins, use `Fourier.Goertzel` with the frequency in cycles per sample. The scalar supplies `sin` and `cos` through `Trigonometry.Circular`; this package conforms no type. Values are dimensionless: a bin k of a transform of n samples taken at rate r is at frequency k·r/n, which is the caller's to compute.

## Laws

- The forward transform is X[k] = Σ x[t] e^(−2πi k t / n) over t from 0 to n − 1, unscaled. The inverse is x[t] = (1 / n) Σ X[k] e^(2πi k t / n), so the inverse of the forward transform is the signal, to rounding.
- Parseval: Σ |X[k]|² = n Σ |x[t]|².
- The transform is linear: the transform of a·x + b·y is a·X + b·Y.
- A pure tone cos(2π k₀ t / n) has n/2 at bins k₀ and n − k₀ and nothing elsewhere; at k₀ = 0 or n/2 the two are one bin of n.
- Any length from one up works. A length whose prime factors are all at most 13 runs by mixed radix (radix 4, then 2, then the odd primes); any other runs by Bluestein's chirp, as a circular convolution over the least power of two at least 2n − 1. Both agree with the defining sum.
- The roots of unity are brought into the first eighth of a turn in integers before any sine or cosine is taken, so every root at a multiple of an eighth of a turn is exact (0, ±1, ±√½), and a root and its conjugate, or two roots a quarter turn apart, differ only in sign and order. Bluestein's chirp index k² is reduced modulo 2n in integers.
- At a mixed-radix length, a transform that multiplies only by exact roots and whose additions are exact is exact: lengths one, two and four on whole numbers, an impulse at the start (one in every bin), and the zero bin of a constant whose n-fold sum is representable (n times it). Bluestein's lengths are within rounding, not exact.
- `forward(real:)` is exactly the first ⌊n / 2⌋ + 1 bins of the complex transform of the same signal. `inverse(real:)` takes those bins, mirrors the rest as conjugates, and gives the real part of the inverse; the imaginary part of bin zero, and for an even count of the last bin, is not read. Its length is the transform's count, since an odd and an even count can share a number of bins.
- `Fourier.Goertzel(frequency:).power(of:)` is |Σ x[t] e^(−2πi f t)|² over a real signal, f in cycles per sample, by Goertzel's recurrence: s₁² + s₂² − 2cos(2πf) s₁s₂. At a bin's frequency k / n it is |X[k]|² of the transform of those n samples; no samples give zero.
- A transform of no samples fails with `.count`, as does a signal or spectrum whose count is not the one the call takes.
