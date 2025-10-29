module EccentricBurstTiming

using Unitful, UnitfulAstro

export BurstTimingModel, evolve!, iterate, get_arrays, line_up_burst_times, get_e_array, get_p_array, get_t_array, get_w_array, get_V3_array

include("constants.jl")
include("utilities.jl")
include("postprocessing.jl")
include("equations.jl")
mutable struct BurstTimingModel{T, angT, massT, lengthT, timeT, floatT, vecangT, veclengthT, vectimeT, vecfloatT}
    Ω₃::angT # change to Ω₃
    ι₃::angT # change to ι₃
    m₁₂::massT 
    m₃::massT
    M::massT
    η::floatT
    η₃::floatT # change to μ₃
    ω₃::angT
    p₃::lengthT
    e₃::floatT
    sqrt_Mp₃⁻³::T
    e::vecfloatT
    p::veclengthT
    t::vectimeT
    ω::vecangT
    Ω::vecangT
    ι::vecangT #
    V₃::vecangT
    eᵢ₋₁::floatT
    pᵢ₋₁::lengthT
    tᵢ₋₁::timeT
    ωᵢ₋₁::angT
    Ωᵢ₋₁::angT
    ιᵢ₋₁::angT
    V₃ᵢ₋₁::angT
    e₀::floatT
    p₀::lengthT
    t₀::timeT
    ω₀::angT
    Ω₀::angT
    ι₀::angT
    V₃₀::angT

    # function BurstTimingModel(;e0 = 0.99, 
    #                            p0 = 30u"Rsun", 
    #                            t0 = 0u"s", 
    #                            m12 = 1u"Msun", 
    #                            eta = 0.20,
    #                            e3 = 0.0,
    #                            w0 = (π/2)u"rad", 
    #                            m3 = 1e7u"Msun", 
    #                            R3 = 1.1e7u"Rsun", 
    #                            V3 = (π/3)u"rad", 
    #                            w3 = 0.0u"rad", 
    #                            i0 = 0.0u"rad", W0=0.0u"rad",
    #                            W3 = 0u"rad", iota3 = 0u"rad")
    #     T1 = typeof(e0)
    #     T2 = typeof(T1[])

    #     p3 = R3*(1 - e3)
    #     MT = m12 + m3
    #     sqrt_Mp₃⁻³ = sqrt(MT/p3^3)

    #     # T, angT, massT, lengthT, timeT, floatT, vecangT, 
    #     # vecangT, veclengthT, vectimeT, vecfloatT
    #     vecs = ([w0], [p0], [t0], [e0], [w0], [W0], [i0], [V3])
    #     Ts = typeof.((sqrt_Mp₃⁻³, w0, m3, p0, t0, e0))
    #     vecTs = typeof.(vecs)

    #     # eta ≡ μ/m₁₂
    #     return new{Ts..., vecTs...}(W3, iota3, m12, m3, m12 + m3, eta, m3 / m12, w3, p3, e3, sqrt_Mp₃⁻³,
    #                         [e0], [p0], [t0], [w0], [W0], [i0], [V3],
    #                         e0, p0, t0, w0, W0, i0, V3, 
    #                         e0, p0, t0, w0, W0, i0, V3)
    # end
end

    function BurstTimingModel(;e0 = 0.99, 
                               p0 = 30.0u"Rsun", 
                               t0 = 0.0u"s", 
                               m12 = 1.0u"Msun", 
                               eta = 0.20,
                               e3 = 0.0,
                               w0 = (π/2)u"rad", 
                               m3 = 1e7u"Msun", 
                               R3 = 1.1e7u"Rsun", 
                               V3 = (π/3)u"rad", 
                               w3 = 0.0u"rad", 
                               i0 = 0.0u"rad", W0=0.0u"rad",
                               W3 = 0.0u"rad", iota3 = 0.0u"rad")
        T1 = typeof(e0)
        T2 = typeof(T1[])

        p3 = R3*(1 - e3)
        MT = m12 + m3
        sqrt_Mp₃⁻³ = sqrt(MT/p3^3)

        # T, angT, massT, lengthT, timeT, floatT, vecangT, 
        # vecangT, veclengthT, vectimeT, vecfloatT
        # vecs = ([w0], [p0], [t0], [e0], [w0], [W0], [i0], [V3])
        # Ts = typeof.((sqrt_Mp₃⁻³, w0, m3, p0, t0, e0))
        # vecTs = typeof.(vecs)

        # eta ≡ μ/m₁₂
        return BurstTimingModel(W3, iota3, m12, m3, m12 + m3, 
                                eta, m3 / m12, w3, p3, e3, sqrt_Mp₃⁻³,
                                [e0], [p0], [t0], [w0], [W0], [i0], [V3],
                                e0, p0, t0, w0, W0, i0, V3, 
                                e0, p0, t0, w0, W0, i0, V3)
    end


# function get_e_next(model)
#     two_body_update = (608*π/15*model.η*sqrt((model.M/model.pᵢ₋₁)^5))*(1 + 121/304*model.eᵢ₋₁ ^2)
#     γ₃³ = γ(e3, V3)^3
#     three_body_perturbation = if model.m₃ > zero(model.m₃)
#         15π/2*model.η₃*model.Cₚ³*γ₃³/(model.M / model.pᵢ₋₁)^3*sin(2*(model.V₃ᵢ₋₁ - model.ωᵢ₋₁))/sqrt((1 - model.eᵢ₋₁^2))^5
#     else
#         0.0
#     end

#     return model.eᵢ₋₁*(1 - two_body_update - three_body_perturbation)
# end

# function get_p_next(model)
#     two_body_update = (128π/5*model.η*sqrt((model.M/model.pᵢ₋₁)^5))*(1 + 7/8*model.eᵢ₋₁^2)
#     γ₃³ = γ(e3, V3)^3    
#     three_body_perturbation = if model.m₃ > zero(model.m₃)
#         15π*model.η₃*model.Cₚ³*γ₃³/(model.M/model.pᵢ₋₁)^3*model.eᵢ₋₁^2*sin(2*(model.V₃ᵢ₋₁ - model.ωᵢ₋₁))/sqrt((1 - model.eᵢ₋₁^2))^7
#     else
#         0.0
#     end

#     return model.pᵢ₋₁*(1 - two_body_update + three_body_perturbation)
# end

# function get_t_next(model)
#     Acoeff = Utils.get_Acoeff(model.η, model.pᵢ₋₁, model.M, model.eᵢ₋₁ )
#     Bcoeff = Utils.get_Bcoeff(model.eᵢ₋₁ )
#     # println(Acoeff, " ", Bcoeff)
#     two_body_update = 2π/√model.M*sqrt((model.pᵢ₋₁ + model.η*sqrt(model.M^5/model.pᵢ₋₁^3)*Acoeff) / 
#                                         (1 - model.eᵢ₋₁^2 + model.η*sqrt((model.M/model.pᵢ₋₁)^5)*Bcoeff))^3
#     γ₃³ = γ(e3, V3)^3    
#     three_body_perturbation = if model.m₃ > zero(model.m₃)
#         1 + model.η₃*model.Cₚ³*γ₃³/(model.M/model.pᵢ₋₁)^3*(5*(4 + (3*model.eᵢ₋₁^2)) + 
#         (96 + (51 * model.eᵢ₋₁^2))*cos(2*(model.V₃ᵢ₋₁ - model.ωᵢ₋₁)))/(16*(1 - model.eᵢ₋₁ ^2)^3)
#     else
#         1.0
#     end

#     return model.tᵢ₋₁ + (two_body_update * three_body_perturbation)
# end

# function get_w_next(model)
#     γ₃³ = γ(e3, V3)^3
#     three_body_perturbation = if model.m₃ > zero(model.m₃)
#         3π/2*model.η₃*model.Cₚ³*γ₃³/(model.M/model.pᵢ₋₁)^3*(1 + (5cos(2*(model.V₃ᵢ₋₁ - model.ωᵢ₋₁))))/sqrt((1 - model.eᵢ₋₁^2))^5
#     else
#         0
#     end

#     return model.ωᵢ₋₁ + three_body_perturbation
# end

# function get_V3_next(model)
#     return model.V₃ᵢ₋₁ + model.sqrt_Mp₃⁻³*(γ(model.e₃, model.V₃ᵢ₋₁)^2*(model.t[end] - model.tᵢ₋₁))
# end

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
            if iszero(model.m₃) || 0.1 > (cbrt((1 + model.eᵢ₋₁ )^(13/2)*model.η₃/(1 - model.eᵢ₋₁ )^2/model.η/((1 + model.eᵢ₋₁ )*√(model.m₁₂ / model.eᵢ₋₁ ))^11)) / R3
                iterate!(model)
            end
            n += 1
        else
            n = Nbursts + 1
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
