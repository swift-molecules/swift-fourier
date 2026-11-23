public import Trigonometry

extension Fourier {
    /// The discrete Fourier transform of one length, planned once and used on any number of signals of that length.
    ///
    /// The forward transform of x is X[k] = Σ x[t] e^(−2πi k t / n) over t from 0 to n − 1, with no scaling; the
    /// inverse is x[t] = (1 / n) Σ X[k] e^(2πi k t / n). Any length works: lengths whose prime factors are all at
    /// most 13 run by mixed radix, any other length by Bluestein's chirp over a power of two.
    public struct Transform<Scalar: BinaryFloatingPoint & Trigonometry.Circular & Sendable>: Sendable {
        @usableFromInline
        internal let storage: Storage

        /// Fails with `.count` for a transform of no samples.
        @inlinable
        public init(count: Cardinal) throws(Fourier.Error) {
            let length = Int(bitPattern: count)
            guard length > 0 else { throw .count }
            self.storage = Storage(count: length)
        }
    }
}

extension Fourier.Transform {
    /// The number of samples the transform takes, and of bins it gives for a complex signal.
    @inlinable
    public var count: Cardinal { Cardinal(UInt(storage.count)) }

    /// The number of bins of a real signal's spectrum: ⌊n / 2⌋ + 1, from zero up to half the rate.
    @inlinable
    public var bins: Cardinal { Cardinal(UInt(storage.count / 2 + 1)) }
}

extension Fourier.Transform {
    /// X[k] = Σ x[t] e^(−2πi k t / n), unscaled. Fails with `.count` unless the signal has the transform's count.
    @inlinable
    public func forward(_ signal: [Complex.Number<Scalar>]) throws(Fourier.Error) -> [Complex.Number<Scalar>] {
        guard signal.count == storage.count else { throw .count }
        return storage.forward(signal)
    }

    /// x[t] = (1 / n) Σ X[k] e^(2πi k t / n), so that the inverse of the forward transform is the signal. Fails with
    /// `.count` unless the spectrum has the transform's count.
    @inlinable
    public func inverse(_ spectrum: [Complex.Number<Scalar>]) throws(Fourier.Error) -> [Complex.Number<Scalar>] {
        guard spectrum.count == storage.count else { throw .count }
        return storage.inverse(spectrum)
    }

    /// The forward transform of a real signal: its first ⌊n / 2⌋ + 1 bins, from zero up to half the rate, which are
    /// the complex transform's bins of the same signal; the rest mirror them as conjugates. Fails with `.count`
    /// unless the signal has the transform's count.
    @inlinable
    public func forward(real signal: [Scalar]) throws(Fourier.Error) -> [Complex.Number<Scalar>] {
        guard signal.count == storage.count else { throw .count }
        let full = storage.forward(signal.map { Complex.Number($0) })
        return Array(full.prefix(storage.count / 2 + 1))
    }

    /// The real signal of a spectrum given as its first ⌊n / 2⌋ + 1 bins, the rest taken as their mirrored
    /// conjugates. The imaginary part of bin zero, and for an even count of the last bin, has no real signal and is
    /// not read. Fails with `.count` unless the spectrum has ⌊n / 2⌋ + 1 bins.
    @inlinable
    public func inverse(real spectrum: [Complex.Number<Scalar>]) throws(Fourier.Error) -> [Scalar] {
        let length = storage.count
        guard spectrum.count == length / 2 + 1 else { throw .count }
        var full = spectrum
        full.reserveCapacity(length)
        for bin in (length / 2 + 1)..<Swift.max(length, length / 2 + 1) {
            full.append(spectrum[length - bin].conjugate)
        }
        return storage.inverse(full).map { $0.realValue }
    }
}
