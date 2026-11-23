public import Trigonometry

extension Fourier {
    /// The transform of a real signal at one frequency, by Goertzel's recurrence: where `Transform` gives every bin
    /// of a length at once, this gives the power at any frequency, on a bin or between bins.
    public struct Goertzel<Scalar: BinaryFloatingPoint & Trigonometry.Circular> {
        /// The frequency in cycles per sample.
        public let frequency: Scalar

        /// 2 cos(2π × frequency).
        @usableFromInline
        internal let coefficient: Scalar

        @inlinable
        public init(frequency: Scalar) {
            self.frequency = frequency
            self.coefficient = 2 * Scalar.cos(2 * Scalar.pi * frequency)
        }
    }
}

extension Fourier.Goertzel: Sendable where Scalar: Sendable {}

extension Fourier.Goertzel {
    /// |Σ x[t] e^(−2πi f t)|² over the signal, f the frequency: s₁² + s₂² − 2cos(2πf) s₁s₂ after the recurrence
    /// s = x[t] + 2cos(2πf) s₁ − s₂ has run over every sample. At a bin's frequency k / n it is |X[k]|² of the
    /// transform of those n samples; a signal of no samples has none.
    @inlinable
    public func power(of signal: some Swift.Sequence<Scalar>) -> Scalar {
        var first: Scalar = 0
        var second: Scalar = 0
        signal.forEach { sample in
            let next = sample + coefficient * first - second
            second = first
            first = next
        }
        return first * first + second * second - coefficient * first * second
    }
}
