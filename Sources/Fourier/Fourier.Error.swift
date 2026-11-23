extension Fourier {
    public enum Error: Swift.Error, Hashable, Sendable {
        /// A transform of no samples, or a signal or spectrum whose count is not the one the transform takes.
        case count
    }
}
