module EccentricBurstTiming

# using Unitful, UnitfulAstro

export BurstTimingModel, evolve!, iterate, get_arrays, line_up_burst_times, get_e_array, get_p_array, get_t_array, get_w_array, get_V3_array

include("constants.jl")
include("utilities.jl")
include("postprocessing.jl")
include("equations.jl")

using .Constants

mutable struct BurstTimingModel{T, vecT}
    Ω₃::T 
    ι₃::T 
    m₁₂::T 
    m₃::T
    M::T
    η::T
    η₃::T # change to μ₃
    ω₃::T
    p₃::T
    e₃::T
    sqrt_Mp₃⁻³::T
    e::vecT
    p::vecT
    t::vecT
    ω::vecT
    Ω::vecT
    ι::vecT #
    V₃::vecT
    eᵢ₋₁::T
    pᵢ₋₁::T
    tᵢ₋₁::T
    ωᵢ₋₁::T
    Ωᵢ₋₁::T
    ιᵢ₋₁::T
    V₃ᵢ₋₁::T
    e₀::T
    p₀::T
    t₀::T
    ω₀::T
    Ω₀::T
    ι₀::T
    V₃₀::T

    function BurstTimingModel(;e0 = 0.99, p0 = 0.0038, t0 = 0, m12 = 60, eta = 0.20, e3=0.0,
                                w0 = π/2, m3 = 6e8, R3 = 1400, V3 = π/3, w3=0.0, i0=0.0, W0=0.0,
                                W3 = 0, iota3 = 0)
        T1 = typeof(e0)
        T2 = typeof(T1[])

        # m12 = SI_to_geometric_mass(m12)
        # m3 =  SI_to_geometric_mass(m3)
        # p0 = SI_to_geometric_length(p0)
        # t0 = SI_to_geometric_time(t0)

        p0 = p0*Rsun_to_m/(m12*Constants.Mconvert)
        R3 = R3*Rsun_to_m/(m12*Constants.Mconvert)
        m3 = m3/m12
        t0 = t0/(m12 * Constants.Msolsec)

        p3 = R3*(1 - e3)
        MT = m12 + m3
        sqrt_Mp₃⁻³ = sqrt(MT/p3^3)

        # eta ≡ μ/m₁₂
        return new{T1, T2}(W3, iota3, m12, m3, m12 + m3, eta, m3 / m12, w3, p3, e3, sqrt_Mp₃⁻³,
                            T1[e0], T1[p0], T1[t0], T1[w0], T1[W0], T1[i0], T1[V3],
                            e0, p0, t0, w0, W0, i0, V3, 
                            e0, p0, t0, w0, W0, i0, V3)
    end
end

SI_to_geometric_length(M, l) = l/(M*Constants.Mconvert)
SI_to_geometric_time(M, t) = t/(M*Constants.Msolsec)

geometric_to_SI_length(M, l) = l*M*Constants.Mconvert
geometric_to_SI_time(M, t) = t*(M*Constants.Msolsec)


function iterate!(model::BurstTimingModel)
    
    e_next = get_e_next(model)
    push!(model.e, e_next)
    
    p_next = get_p_next(model)
    push!(model.p, p_next)
    
    t_next = get_t_next(model)
    push!(model.t, t_next)
    
    ω_next = get_ω_next(model)
    push!(model.ω, ω_next)

    ι_next = get_ι_next(model)
    push!(model.ι, ι_next)

    Ω_next = get_Ω_next(model)
    push!(model.Ω, Ω_next)
    
    V₃_next = get_V3_next(model)
    push!(model.V₃, V₃_next)

    model.eᵢ₋₁  = e_next
    model.pᵢ₋₁ = p_next
    model.tᵢ₋₁ = t_next
    model.ωᵢ₋₁ = ω_next
    model.V₃ᵢ₋₁ = V₃_next
    model.ιᵢ₋₁ = ι_next
    nothing
end

function observed_burst_time_offsets_due_to_com_motion(model)
    # com_motion_induced_offset = @. model.R3/model.M*sin(model.iota)*cos(model.V3 + model.lambda)
    # com_motion_induced_offset = self.R3 / self.M * np.sin(self.iota) * np.cos(np.array(self.V3) + self.Lambda)

    return @. model.m₃/model.M*model.p₃*sin(model.ι₃)/(1 + model.e*cos(model.V₃))*cos(model.V₃ + model.ω₃)
end

function evolve!(model::BurstTimingModel, n_bursts)
    n = 0
    while (n <= n_bursts) && (model.eᵢ₋₁  > 0.7) && (model.eᵢ₋₁  < 1.0) && (model.pᵢ₋₁ > 0)
        if model.pᵢ₋₁ > 2*(3 + model.eᵢ₋₁ )
            #if model.m₃ == 0.0 || 0.1 > (np.power(
            #        (1 + model.eᵢ₋₁ ) ** 9 * model.η₃ /
            #        (1 - model.eᵢ₋₁ ) ** 2 /
            #        ((1 + model.eᵢ₋₁ ) * np.power(model.M / model.eᵢ₋₁ , 0.5)) ** 11,
            #         1/3)) / model.R3:
            R3 = get_R(model)
            if iszero(model.m₃) || 0.1 > (cbrt((1 + model.eᵢ₋₁ )^(13/2)*model.η₃/(1 - model.eᵢ₋₁ )^2/model.η/((1 + model.eᵢ₋₁ )*√(1 / model.eᵢ₋₁ ))^11)) / R3
                iterate!(model)
            end
            n += 1
        else
            n = n_bursts + 1
        end
    end
    # Cut off last point to ensure within region of validity
    model.e = model.e[1:end-1]
    model.p = model.p[1:end-1]
    model.t = model.t[1:end-1]
    model.ω = model.ω[1:end-1]
    model.ι = model.ι[1:end-1]
    model.Ω = model.Ω[1:end-1]
    model.V₃ = model.V₃[1:end-1]

    # Add observed time offset due to system inclination
    if model.m₃ > 0
        additional_offset = observed_burst_time_offsets_due_to_com_motion(model)
        model.t .-= additional_offset
    end
    # Normalise so that first burst happens at time 0
    #model.t = [x - model.t[0] for x in model.t]
end

function Base.getproperty(m::BurstTimingModel, s::Symbol)
    s_str = string(s)
    s_str = replace(s_str, "w" => "ω")
    s_str = replace(s_str, "W" => "Ω")
    s_str = replace(s_str, "0" => "₀")
    s_str = replace(s_str, "3" => "₃")
    s_str = replace(s_str, "1" => "₁")
    s_str = replace(s_str, "2" => "₂")
    return getfield(m, Symbol(s_str))
end

function Base.show(io::IO, m::BurstTimingModel)
    props_to_show = (:m₁₂, :η, :e₀, :p₀, :ω₀, :Ω₀, :ι₀, :m₃, :ω₃, :p₃, :e₃, :V₃₀, :Ω₃, :ι₃)
    for prop in props_to_show
        println(io, prop, " ", getproperty(m, prop))
    end
    nothing
end


end # module EccentricBurstTiming
