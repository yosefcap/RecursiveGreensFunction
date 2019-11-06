#=@with_kw_noshow  struct HamiltonianParam
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
=#
@with_kw_noshow  struct HamiltonianParam
    #user mandatory
    L::Int
    B::Float64
    #user optional

    #hopping parameters
    t::Float64 = 1.0
    μ::Float64  = 0.0
    η::Float64 = 0.01 # imagenary energy

    #CDW order
    CDW::Float64 = 0.0
    CDW_period::Float64 = 1.0
    CDW_type::Int8 = 1 # 1 for s-wavw CDW and -1 for d-wave CDW

    #PDW order
    Δ::Float64 = 0.0
    PDW_period::Float64 = 1.0
    PDW_type::Int8 = -1 # 1 for s-wavw PDW and -1 for d-wave PDW

end

struct Hamiltonian
    xi::Array{Float64,2}
    Δ_l::Array{Float64,1} #local sc terms
    Δ_h::Array{Float64,1} #hopping sc terms
end

struct Greens
    Left_Green::Array{Float64,3}
    Right_Green::Array{Float64,3}
    Full_Green::Array{Float64,3}
end

function init_Hamiltonian(H::Hamiltonian,L::Int,B::Float64, p_y::Float64;  kwargs...)

    p = HamiltonianParam(; kwargs...)
    Q_cdw=2*pi/p.CDW_period
    Q_pdw=2*pi/p.PDW_period
    for c in 1:L
        H.xi[c,1] = -2*p.t*cos(p_y+2*pi*B*c) + 2*p.CDW*cos(Q_cdw*x+ϕ[c]) - p.μ -
        H.xi[c,2] = -2*p.t*cos(-p_y+2*pi*B*c) + 2*p.CDW*cos(Q_cdw*x+ϕ[c]) - p.μ -
        H.Δ_l[c] =  p.PDW_type*2*p.Δ*cos(Q_pdw*(c-0.5))*cos(p_y)
        H.Δ_h[c] =  p.Δ*cos(Q_pdw*c)
    end
end

function Recursive_Greens_function(H::Hamiltonian,G::Greens,L::Int)
    h_1 = get_hamiltonian_block(1,1)
    h_L = get_hamiltonian_block(L,L)

    G.Left_Green[:,:,1] = inv(h_1)
    G.Right_Green[:,:,L]=inv(h_L)
    for c in 2:L
        h_c = get_hamiltonian_block(c,c)
        h_pm = get_hamiltonian_block(c,c-1)
        h_mp = get_hamiltonian_block(c-1,c)
        G.Left_Green[:,:,c] =  inv(h_c - h_pm*G.Left_Green[:,:,c-1]*h_mp)

        h_c = get_hamiltonian_block(L+1-c,L+1-c)
        v_pm = get_hamiltonian_block(L+2-c,L+1-c)
        v_mp = get_hamiltonian_block(L+1-c,L+2-c)
        G.Right_Green[:,:,L+1-c] = inv(h_c - v_mp*G.Right_Green[:,:,L+2-c]*v_pm)
    end

    h_1 = get_hamiltonian_block(1,1)
    v_pm = get_hamiltonian_block(2,1)
    v_mp = get_hamiltonian_block(1,2)
    G.Full_Green[:,:,c] = inv(h_c - v_mp*G.Right_Green[:,:,2]*v_pm)

    h_L = get_hamiltonian_block(L,L)
    h_pm = get_hamiltonian_block(L,L-1)
    h_mp = get_hamiltonian_block(L-1,L)

    G.Full_Green[:,:,c] = inv(h_c - h_pm*G.Left_Green[:,:,L-1]*h_mp)

    for c in 2:(L-1))
        h_c = get_hamiltonian_block(c,c)
        h_pm = get_hamiltonian_block(c,c-1)
        h_mp = get_hamiltonian_block(c-1,c)
        v_mp = get_hamiltonian_block(c,c+1)
        v_pm = get_hamiltonian_block(c+1,c)

        G.Full_Green[:,:,c] = inv(h_c - h_pm*G.Left_Green[:,:,c-1]*h_mp - v_mp*G.Right_Green[:,:,c+1]*v_pm)
    end

end
