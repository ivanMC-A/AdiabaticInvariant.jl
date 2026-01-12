using AdiabaticInvariant
using Plots

A = AdiabaticInvariant.chebA(3)
At = A'

W = At*A
size(W)[1]