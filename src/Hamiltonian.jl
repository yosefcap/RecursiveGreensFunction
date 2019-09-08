@with_kw_noshow  struct HamiltonianParam
    #user mandatory
    slice_width::Int
    L::Int
    B::Float64
    #user optional

    #hopping parameters
    t_1::Float64 = 1.0
    t_2::Float64 = 0.0
    j_3::Float64 = 0.0
    mu::Float64  = 0.0

    #CDW order
    CDW_x::Float64 = 0.0
    CDW_y::Float64 = 0.0
    CDW_period_x::Float64 = 1.0
    CDW_period_y::Float64 = 1.0
    CDW_type::Int8 = 1 # 1 for s-wavw CDW and -1 for d-wave CDW

    #PDW order
    PDW_x::Float64 = 0.0
    PDW_y::Float64 = 0.0
    PDW_period_x::Float64 = 1.0
    PDW_period_y::Float64 = 1.0
    PDW_type::Int8 = 1 # 1 for s-wavw CDW and -1 for d-wave CDW

end
function Hamiltonian(L::int,slice_width::Int,B::Float64;  kwargs...)

    p = HamiltonianParam(; kwargs...)
    l = SquareLattice(L,slice_width)
end
