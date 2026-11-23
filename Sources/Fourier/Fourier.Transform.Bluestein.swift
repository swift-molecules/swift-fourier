extension Fourier.Transform {
    /// Bluestein's chirp: a transform of any length n as a circular convolution of length m, the least power of two
    /// at least 2n − 1, since k t = (k² + t² − (k − t)²) / 2.
    @usableFromInline
    internal struct Bluestein: Sendable {
        /// e^(−πi k² / n) for k from 0 to n − 1, its index k² taken modulo 2n in integers.
        @usableFromInline
        internal let chirp: [Complex.Number<Scalar>]

        /// The forward transform of the conjugate chirp laid circularly over m.
        @usableFromInline
        internal let response: [Complex.Number<Scalar>]

        @usableFromInline
        internal let inner: Fourier.Transform<Scalar>

        @usableFromInline
        internal init(count: Int) {
            var length = 1
            while length < 2 * count - 1 { length *= 2 }
            let chirp = (0..<count).map { k in Storage.root((k * k) % (2 * count), of: 2 * count) }
            var kernel = [Complex.Number<Scalar>](repeating: .zero, count: length)
            kernel[0] = chirp[0].conjugate
            for k in 1..<Swift.max(count, 1) {
                kernel[k] = chirp[k].conjugate
                kernel[length - k] = chirp[k].conjugate
            }
            let inner = Fourier.Transform<Scalar>(storage: Storage(count: length))
            self.chirp = chirp
            self.response = inner.storage.forward(kernel)
            self.inner = inner
        }

        @inlinable
        internal func forward(_ signal: [Complex.Number<Scalar>]) -> [Complex.Number<Scalar>] {
            let length = response.count
            var a = [Complex.Number<Scalar>](repeating: .zero, count: length)
            for k in signal.indices {
                a[k] = signal[k] * chirp[k]
            }
            var spectrum = inner.storage.forward(a)
            for index in spectrum.indices {
                spectrum[index] = spectrum[index] * response[index]
            }
            let convolution = inner.storage.inverse(spectrum)
            return signal.indices.map { convolution[$0] * chirp[$0] }
        }
    }
}

extension Fourier.Transform {
    @usableFromInline
    internal init(storage: Storage) {
        self.storage = storage
    }
}
