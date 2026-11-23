import Fourier
import Numeric
import Testing

enum `Fourier.Transform Tests` {}

extension `Fourier.Transform Tests` {
    typealias C = Complex.Number<Double>

    /// The transform by its definition, a sum of n terms for each of n bins, with each angle reduced in integers.
    static func naive(_ x: [C]) -> [C] {
        let n = x.count
        return (0..<n).map { k in
            x.indices.reduce(C.zero) { sum, t in
                let angle = -2 * Double.pi * Double((k * t) % n) / Double(n)
                return sum + x[t] * C(Double.cos(angle), Double.sin(angle))
            }
        }
    }

    static func signal(_ n: Int, seed: Int = 1) -> [C] {
        (0..<n).map { t in
            let a = Double((t * 7_919 + seed * 104_729) % 1_013) / 1_013 - 0.5
            let b = Double((t * 6_037 + seed * 15_485) % 997) / 997 - 0.5
            return C(a, b)
        }
    }

    static func distance(_ x: [C], _ y: [C]) -> Double {
        zip(x, y).reduce(0) { Swift.max($0, ($1.0 - $1.1).magnitude().value) }
    }

    static func energy(_ x: [C]) -> Double {
        x.reduce(0) { $0 + $1.magnitude.squared.value }
    }

    static var lengths: [Int] { Array(1...100) + [127, 128, 204, 257, 900, 1_009] }
}

extension `Fourier.Transform Tests` {
    @Suite
    struct `The transform keeps its laws` {
        @Test
        func `The forward transform is the sum that defines it, at every length`() throws {
            for n in lengths {
                let x = signal(n)
                let transform = try Fourier.Transform<Double>(count: Cardinal(UInt(n)))
                let error = distance(try transform.forward(x), naive(x))
                #expect(error <= 1e-12 * Double(n), "length \(n): \(error)")
            }
        }

        @Test
        func `The inverse of the forward transform is the signal`() throws {
            for n in lengths {
                let x = signal(n, seed: 3)
                let transform = try Fourier.Transform<Double>(count: Cardinal(UInt(n)))
                let error = distance(try transform.inverse(try transform.forward(x)), x)
                #expect(error <= 1e-14 * Double(n), "length \(n): \(error)")
            }
        }

        @Test
        func `Parseval: the energy of the spectrum is n times the energy of the signal`() throws {
            for n in lengths {
                let x = signal(n, seed: 5)
                let transform = try Fourier.Transform<Double>(count: Cardinal(UInt(n)))
                let spectrum = energy(try transform.forward(x))
                let expected = Double(n) * energy(x)
                #expect(Swift.abs(spectrum - expected) <= 1e-12 * expected, "length \(n)")
            }
        }

        @Test
        func `The transform is linear`() throws {
            for n in [12, 64, 97, 210] {
                let x = signal(n, seed: 7)
                let y = signal(n, seed: 11)
                let a = C(0.75, -1.25)
                let b = C(-2, 0.5)
                let transform = try Fourier.Transform<Double>(count: Cardinal(UInt(n)))
                let mixed = try transform.forward(zip(x, y).map { a * $0 + b * $1 })
                let apart = zip(try transform.forward(x), try transform.forward(y)).map { a * $0 + b * $1 }
                #expect(distance(mixed, apart) <= 1e-12 * Double(n), "length \(n)")
            }
        }

        @Test
        func `A pure tone at bin k is n over two at k and at n − k, and nothing elsewhere`() throws {
            for (n, k) in [(64, 5), (100, 3), (97, 10), (257, 40)] {
                let x = (0..<n).map { t in C(Double.cos(2 * Double.pi * Double((k * t) % n) / Double(n))) }
                let transform = try Fourier.Transform<Double>(count: Cardinal(UInt(n)))
                let spectrum = try transform.forward(x)
                (0..<n).forEach { bin in
                    let expected = bin == k || bin == n - k ? Double(n) / 2 : 0
                    #expect((spectrum[bin] - C(expected)).magnitude().value <= 1e-11 * Double(n), "length \(n) bin \(bin)")
                }
            }
        }
    }
}

extension `Fourier.Transform Tests` {
    @Suite
    struct `The transform is exact where it can be` {
        @Test
        func `An impulse at the start gives one in every bin, exactly, at a radix length`() throws {
            for n in [1, 2, 3, 8, 12, 60, 128, 1_001] {
                var x = [C](repeating: .zero, count: n)
                x[0] = 1
                let transform = try Fourier.Transform<Double>(count: Cardinal(UInt(n)))
                #expect(try transform.forward(x) == [C](repeating: 1, count: n), "length \(n)")
            }
        }

        @Test
        func `A constant gives n times itself at bin zero, exactly`() throws {
            for n in [1, 2, 6, 16, 45, 128] {
                let transform = try Fourier.Transform<Double>(count: Cardinal(UInt(n)))
                let spectrum = try transform.forward([C](repeating: C(3, -2), count: n))
                #expect(spectrum[0] == C(3 * Double(n), -2 * Double(n)), "length \(n)")
            }
        }

        @Test
        func `Lengths one, two and four only add, so whole numbers come out exact`() throws {
            let four = try Fourier.Transform<Double>(count: 4)
            #expect(try four.forward([1, 2, 3, 4]) == [C(10, 0), C(-2, 2), C(-2, 0), C(-2, -2)])
            let two = try Fourier.Transform<Double>(count: 2)
            #expect(try two.forward([C(1, 2), C(3, 5)]) == [C(4, 7), C(-2, -3)])
            let one = try Fourier.Transform<Double>(count: 1)
            #expect(try one.forward([C(5, -1)]) == [C(5, -1)])
        }

        @Test
        func `Eighths of a turn are exact`() throws {
            let transform = try Fourier.Transform<Double>(count: 8)
            var x = [C](repeating: .zero, count: 8)
            x[1] = 1
            let spectrum = try transform.forward(x)
            let root = (0.5 as Double).squareRoot()
            #expect(spectrum[1] == C(root, -root))
            #expect(spectrum[2] == C(0, -1))
            #expect(spectrum[3] == C(-root, -root))
            #expect(spectrum[4] == C(-1, 0))
            #expect(spectrum[6] == C(0, 1))
            #expect(spectrum[7] == C(root, root))
        }
    }
}

extension `Fourier.Transform Tests` {
    @Suite
    struct `A real signal has half a spectrum` {
        @Test
        func `The real transform is the complex transform's first half`() throws {
            for n in [1, 2, 7, 16, 45, 97, 100] {
                let x = signal(n).map(\.realValue)
                let transform = try Fourier.Transform<Double>(count: Cardinal(UInt(n)))
                let half = try transform.forward(real: x)
                let bins = n / 2 + 1
                #expect(Cardinal(UInt(half.count)) == transform.bins)
                #expect(transform.bins == Cardinal(UInt(bins)))
                #expect(half == Array(try transform.forward(x.map { C($0) }).prefix(n / 2 + 1)))
            }
        }

        @Test
        func `The real inverse of the real forward transform is the signal`() throws {
            for n in [1, 2, 7, 16, 45, 97, 100] {
                let x = signal(n, seed: 9).map(\.realValue)
                let transform = try Fourier.Transform<Double>(count: Cardinal(UInt(n)))
                let back = try transform.inverse(real: try transform.forward(real: x))
                let error = zip(back, x).reduce(0.0) { Swift.max($0, Swift.abs($1.0 - $1.1)) }
                #expect(error <= 1e-14 * Double(n), "length \(n)")
            }
        }

        @Test
        func `The imaginary parts of bin zero and of the last bin of an even count are not read`() throws {
            let transform = try Fourier.Transform<Double>(count: 8)
            let spectrum = try transform.forward(real: [1, 0, -2, 3, 0.5, 0, 1, -1])
            var changed = spectrum
            changed[0] = C(spectrum[0].realValue, 9)
            changed[4] = C(spectrum[4].realValue, -9)
            #expect(try transform.inverse(real: changed) == (try transform.inverse(real: spectrum)))
        }
    }
}

extension `Fourier.Transform Tests` {
    @Suite
    struct `Counts are checked` {
        @Test
        func `A transform of no samples fails`() {
            #expect(throws: Fourier.Error.count) { try Fourier.Transform<Double>(count: 0) }
        }

        @Test
        func `A signal or spectrum of another count fails`() throws {
            let transform = try Fourier.Transform<Double>(count: 6)
            #expect(throws: Fourier.Error.count) { try transform.forward([C](repeating: 0, count: 5)) }
            #expect(throws: Fourier.Error.count) { try transform.inverse([C](repeating: 0, count: 7)) }
            #expect(throws: Fourier.Error.count) { try transform.forward(real: [0, 1]) }
            #expect(throws: Fourier.Error.count) { try transform.inverse(real: [C](repeating: 0, count: 6)) }
        }

        @Test
        func `Float works the same`() throws {
            let transform = try Fourier.Transform<Float>(count: 30)
            let x = (0..<30).map { Complex.Number<Float>(Float(($0 * 13) % 7) - 3, 0) }
            let back = try transform.inverse(try transform.forward(x))
            #expect(zip(back, x).allSatisfy { ($0 - $1).magnitude().value < 1e-5 })
        }
    }
}
