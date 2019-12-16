using RecursiveGreensFunction
using Test, Plots, FFTW, LinearAlgebra

L=2
#@testset DQMC.jl begin
m = HoppingBFModel(dims=2,L=2,J=1, h=0.2, α=0.,μ=1.)
mc=DQMC_bond(m,beta=1)
slices=  mc.p.slices
initialize_stack(mc) # redundant ?!
build_stack(mc)
kkk
