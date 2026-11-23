extension Fourier.Transform {
    /// The plan of one length: its radices and roots of unity, or Bluestein's chirp when a prime factor is too large
    /// for a radix.
    @usableFromInline
    internal final class Storage: Sendable {
        @usableFromInline
        internal let count: Int

        /// The radices, 4 first, then 2, then odd primes up to 13; empty for Bluestein.
        @usableFromInline
        internal let factors: [Int]

        /// e^(−2πi r / count) for every r from 0 to count − 1; empty for Bluestein.
        @usableFromInline
        internal let roots: [Complex.Number<Scalar>]

        @usableFromInline
        internal let bluestein: Bluestein?

        @usableFromInline
        internal init(count: Int) {
            self.count = count
            if let factors = Self.factors(of: count) {
                self.factors = factors
                self.roots = (0..<count).map { Self.root($0, of: count) }
                self.bluestein = nil
            } else {
                self.factors = []
                self.roots = []
                self.bluestein = Bluestein(count: count)
            }
        }
    }
}

extension Fourier.Transform.Storage {
    /// The largest prime a radix handles; a length with a larger prime factor goes through Bluestein.
    @usableFromInline
    internal static var largest: Int { 13 }

    @usableFromInline
    internal static func factors(of count: Int) -> [Int]? {
        var rest = count
        var result: [Int] = []
        while rest % 4 == 0 {
            result.append(4)
            rest /= 4
        }
        if rest % 2 == 0 {
            result.append(2)
            rest /= 2
        }
        var prime = 3
        while rest > 1, prime <= largest {
            while rest % prime == 0 {
                result.append(prime)
                rest /= prime
            }
            prime += 2
        }
        return rest == 1 ? result : nil
    }

    /// e^(−2πi r / d), exact at every eighth of a turn and symmetric under conjugation and quarter turns: the angle
    /// is brought into the first eighth of a turn in integers before any sine or cosine is taken.
    @usableFromInline
    internal static func root(_ r: Int, of d: Int) -> Complex.Number<Scalar> {
        let turn = ((r % d) + d) % d
        let quadrant = (4 * turn) / d
        let rest = (4 * turn) % d
        let cosine: Scalar
        let sine: Scalar
        if rest == 0 {
            cosine = 1
            sine = 0
        } else if 2 * rest == d {
            let half: Scalar = 0.5
            cosine = half.squareRoot()
            sine = cosine
        } else if 2 * rest < d {
            let angle = Scalar.pi / 2 * Scalar(rest) / Scalar(d)
            cosine = Scalar.cos(angle)
            sine = Scalar.sin(angle)
        } else {
            let angle = Scalar.pi / 2 * Scalar(d - rest) / Scalar(d)
            cosine = Scalar.sin(angle)
            sine = Scalar.cos(angle)
        }
        switch quadrant {
        case 0: return Complex.Number(cosine, -sine)
        case 1: return Complex.Number(-sine, -cosine)
        case 2: return Complex.Number(-cosine, sine)
        default: return Complex.Number(sine, cosine)
        }
    }
}

extension Fourier.Transform.Storage {
    @inlinable
    internal func forward(_ signal: [Complex.Number<Scalar>]) -> [Complex.Number<Scalar>] {
        if let bluestein {
            return bluestein.forward(signal)
        }
        var output = signal
        var scratch = [Complex.Number<Scalar>](repeating: .zero, count: factors.max() ?? 1)
        radix(signal, from: 0, stride: 1, into: &output, at: 0, count: count, factor: 0, scratch: &scratch)
        return output
    }

    @inlinable
    internal func inverse(_ spectrum: [Complex.Number<Scalar>]) -> [Complex.Number<Scalar>] {
        let length = Scalar(count)
        return forward(spectrum.map { $0.conjugate }).map { Complex.Number($0.realValue / length, -$0.imaginaryValue / length) }
    }

    /// One stage of the mixed-radix transform by decimation in time: the `count` samples of `input` from `start`,
    /// `stride` apart, transformed into `output` from `offset`.
    @inlinable
    internal func radix(
        _ input: [Complex.Number<Scalar>],
        from start: Int,
        stride: Int,
        into output: inout [Complex.Number<Scalar>],
        at offset: Int,
        count length: Int,
        factor: Int,
        scratch: inout [Complex.Number<Scalar>]
    ) {
        guard length > 1 else {
            output[offset] = input[start]
            return
        }
        let p = factors[factor]
        let m = length / p
        for j in 0..<p {
            radix(input, from: start + j * stride, stride: stride * p, into: &output, at: offset + j * m, count: m, factor: factor + 1, scratch: &scratch)
        }
        let step = count / length
        switch p {
        case 2:
            for k in 0..<m {
                let a = output[offset + k]
                let b = output[offset + m + k] * roots[k * step]
                output[offset + k] = a + b
                output[offset + m + k] = a - b
            }
        case 4:
            for k in 0..<m {
                let a0 = output[offset + k]
                let a1 = output[offset + m + k] * roots[k * step]
                let a2 = output[offset + 2 * m + k] * roots[2 * k * step]
                let a3 = output[offset + 3 * m + k] * roots[3 * k * step]
                let t0 = a0 + a2
                let t1 = a0 - a2
                let t2 = a1 + a3
                let d = a1 - a3
                let t3 = Complex.Number(d.imaginaryValue, -d.realValue)
                output[offset + k] = t0 + t2
                output[offset + m + k] = t1 + t3
                output[offset + 2 * m + k] = t0 - t2
                output[offset + 3 * m + k] = t1 - t3
            }
        default:
            let turn = count / p
            for k in 0..<m {
                for j in 0..<p {
                    scratch[j] = output[offset + j * m + k] * roots[j * k * step]
                }
                for q in 0..<p {
                    var sum = scratch[0]
                    for j in 1..<p {
                        sum += scratch[j] * roots[((j * q) % p) * turn]
                    }
                    output[offset + q * m + k] = sum
                }
            }
        }
    }
}
