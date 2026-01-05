using AdiabaticInvariant
using Plots

A = AdiabaticInvariant.chebA(10)
At = transpose(A)

At*A