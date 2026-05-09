using Distributions
using Statistics

Φ(y) = cdf(Normal(), y)

p = range(0, 1; length=Int(1e8))
y = quantile.(Normal(0, sqrt(1 / 3)), p)
_I = mean(Φ.(y) .* Φ.(√2y))

@show _I;