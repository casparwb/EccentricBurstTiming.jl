module EccentricBurstTiming

# using Unitful, UnitfulAstro

export BurstTimingModel, evolve!, iterate, get_arrays, line_up_burst_times, get_e_array, get_p_array, get_t_array, get_w_array, get_V3_array, get_a_array
export get_e, get_p, get_t, get_w, get_V3, get_a, get_masses

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

    function BurstTimingModel(;
            e0 = 0.99, a0 = 0.19, t0 = 0, m12 = 60, eta = 0.2, e3 = 0.0,
            w0 = 0.0, m3 = 6.0e8, R3 = 1400, V3 = π / 3, w3 = 0.0, i0 = 0.0, W0 = 0.0,
            W3 = 0, iota3 = 0
        )
        T1 = typeof(e0)
        T2 = typeof(T1[])

        p0 = (1 - e0^2) * a0 * Rsun_to_m / (m12 * Constants.Mconvert)
        R3 = R3 * Rsun_to_m / (m12 * Constants.Mconvert)
        m3 = m3 / m12
        t0 = t0 / (m12 * Constants.Msolsec)

        p3 = R3 * (1 - e3)
        MT = m12 + m3
        sqrt_Mp₃⁻³ = sqrt(MT / p3^3)

        # eta ≡ μ/m₁₂
        return new{T1, T2}(
            W3, iota3, m12, m3, m12 + m3, eta, m3 / m12, w3, p3, e3, sqrt_Mp₃⁻³,
            T1[e0], T1[p0], T1[t0], T1[w0], T1[W0], T1[i0], T1[V3],
            e0, p0, t0, w0, W0, i0, V3,
            e0, p0, t0, w0, W0, i0, V3
        )
    end
end

function get_masses(model::BurstTimingModel)
    M = model.m₁₂
    μ = model.η*M

    fac = √(M^2 - 4μ*M)/2

    return M/2 + fac, M/2 - fac
end

SI_to_geometric_length(M, l) = l / (M * Constants.Mconvert)
SI_to_geometric_time(M, t) = t / (M * Constants.Msolsec)

geometric_to_SI_length(M, l) = l * M * Constants.Mconvert
geometric_to_SI_time(M, t) = t * (M * Constants.Msolsec)

# const SYMS = (:η, :ιᵢ₋₁, :V₃ᵢ₋₁, :Ωᵢ₋₁, :ωᵢ₋₁, :ω₃, :eᵢ₋₁, :pᵢ₋₁, :m₃, :tᵢ₋₁, :p₃, :e₃)
# get_params(m) = map(x -> getfield(m, x), SYMS)

function iterate!(model::BurstTimingModel; do_outer_orbit = true, save_params=true)

    (; η, ιᵢ₋₁, V₃ᵢ₋₁, Ωᵢ₋₁, ωᵢ₋₁, ω₃, eᵢ₋₁, pᵢ₋₁, m₃, tᵢ₋₁, p₃, e₃) = model
    # η, ιᵢ₋₁, V₃ᵢ₋₁, Ωᵢ₋₁, ωᵢ₋₁, ω₃, eᵢ₋₁, pᵢ₋₁, m₃, tᵢ₋₁, p₃, e₃ = get_params(model)

    R = get_R(p₃, e₃, V₃ᵢ₋₁)

    e_next = get_e_next(η, ιᵢ₋₁, V₃ᵢ₋₁, Ωᵢ₋₁, ωᵢ₋₁, ω₃, eᵢ₋₁, pᵢ₋₁, m₃, R)
    p_next = get_p_next(η, ιᵢ₋₁, V₃ᵢ₋₁, Ωᵢ₋₁, ωᵢ₋₁, ω₃, eᵢ₋₁, pᵢ₋₁, m₃, R)
    t_next = get_t_next(η, ιᵢ₋₁, V₃ᵢ₋₁, Ωᵢ₋₁, ωᵢ₋₁, ω₃, eᵢ₋₁, pᵢ₋₁, m₃, R, tᵢ₋₁)
    ω_next = get_ω_next(η, ιᵢ₋₁, V₃ᵢ₋₁, Ωᵢ₋₁, ωᵢ₋₁, ω₃, eᵢ₋₁, pᵢ₋₁, m₃, R)
    ι_next = get_ι_next(η, ιᵢ₋₁, V₃ᵢ₋₁, Ωᵢ₋₁, ωᵢ₋₁, ω₃, eᵢ₋₁, pᵢ₋₁, m₃, R)
    Ω_next = get_Ω_next(η, ιᵢ₋₁, V₃ᵢ₋₁, Ωᵢ₋₁, ωᵢ₋₁, ω₃, eᵢ₋₁, pᵢ₋₁, m₃, R)

    V₃_next = if do_outer_orbit
        get_V3_next(model.sqrt_Mp₃⁻³, e₃, V₃ᵢ₋₁, model.t[end], tᵢ₋₁)
    else
        model.V₃ᵢ₋₁
    end
    
    if save_params
        push!(model.e, e_next)
        push!(model.p, p_next)
        push!(model.t, t_next)
        push!(model.ω, ω_next)
        push!(model.ι, ι_next)
        push!(model.Ω, Ω_next)
        push!(model.V₃, V₃_next)
    end

    model.eᵢ₋₁ = e_next
    model.pᵢ₋₁ = p_next
    model.tᵢ₋₁ = t_next
    model.ωᵢ₋₁ = ω_next
    model.Ωᵢ₋₁ = Ω_next
    model.V₃ᵢ₋₁ = V₃_next
    model.ιᵢ₋₁ = ι_next
    return nothing
end

function observed_burst_time_offsets_due_to_com_motion(model)
    # com_motion_induced_offset = @. model.R3/model.M*sin(model.iota)*cos(model.V3 + model.lambda)
    # com_motion_induced_offset = self.R3 / self.M * np.sin(self.iota) * np.cos(np.array(self.V3) + self.Lambda)

    return @. model.m₃ / model.M * model.p₃ * sin(model.ι₃) / (1 + model.e * cos(model.V₃)) * cos(model.V₃ + model.ω₃)
end

function evolve!(model::BurstTimingModel, n_bursts; t_max=Inf, e_min=0.7, f_GW_max=Inf, do_outer_orbit = true, verbose=false, save_params=true)
    n = 0

    # to_s = (model.m12 * Constants.Msolsec)
    t_convert = (model.m₁₂ * Constants.Msolsec)
    s_to_code_units = 1/t_convert
    t_max = t_max*s_to_code_units

    π_inv = 1/π

    f_GW_max_code_units = t_convert*f_GW_max
    while (n <= n_bursts) 

        if model.tᵢ₋₁ >= t_max 
            verbose && @info "Stopping condition: t >= $t_max"
            break
        end

        if (model.eᵢ₋₁ < e_min) 
            verbose && @info "Stopping condition: Eccentricity < $e_min"
            break
        end

        if (model.eᵢ₋₁ > 1.0) 
            verbose && @info "Stopping condition: Eccentricity > 1"
            break
        end

        if (model.pᵢ₋₁ < 0)
            verbose && @info "Stopping condition: SLR > 0"
            break
        end

        if (π_inv*(1 + model.eᵢ₋₁)^(1.195)/sqrt(model.pᵢ₋₁^3)) >= f_GW_max_code_units
            verbose && @info "Stopping condition: fGW > fGW_max"
            break
        end

        if model.pᵢ₋₁ > 2 * (3 + model.eᵢ₋₁)
            R3 = get_R(model)
            if iszero(model.m₃) || 0.1 > (cbrt((1 + model.eᵢ₋₁)^(13 / 2) * model.η₃ / (1 - model.eᵢ₋₁)^2 / model.η / ((1 + model.eᵢ₋₁) * √(1 / model.eᵢ₋₁))^11)) / R3
                iterate!(model, do_outer_orbit = do_outer_orbit, save_params=save_params)
            end
            n += 1
        else
            n = n_bursts + 1
        end
    end

    # Cut off last point to ensure within region of validity

    region_of_validity = min(findlast(x -> x < one(x), model.e), findlast(x -> x > zero(x), model.p))

    model.e = model.e[1:region_of_validity]
    model.p = model.p[1:region_of_validity]
    model.t = model.t[1:region_of_validity]
    model.ω = model.ω[1:region_of_validity]
    model.ι = model.ι[1:region_of_validity]
    model.Ω = model.Ω[1:region_of_validity]
    model.V₃ = model.V₃[1:region_of_validity]

    # Add observed time offset due to system inclination
    if model.m₃ > 0
        additional_offset = observed_burst_time_offsets_due_to_com_motion(model)
        model.t .-= additional_offset
    end
    
    nothing
end

function evolve_to_fGW!(model::BurstTimingModel, fGW_max; max_bursts=1_000, do_outer_orbit = true)
    M_SI = model.m₁₂ * Constants.Msol

    burst_number = 1
    while burst_number < max_bursts
        iterate!(model, do_outer_orbit = do_outer_orbit, save_params=false)

        model.e[1] = model.eᵢ₋₁
        model.p[1] = model.pᵢ₋₁
        model.t[1] = model.tᵢ₋₁
        model.ω[1] = model.ωᵢ₋₁
        model.Ω[1] = model.Ωᵢ₋₁
        model.V₃[1] = model.V₃ᵢ₋₁
        model.ι[1] = model.ιᵢ₋₁

        fGW = peak_f_GW(model)
        if fGW >= fGW_max
            break
        end

        burst_number += 1
    end

    model.e₀ = model.eᵢ₋₁
    model.p₀ = model.pᵢ₋₁
    model.t₀ = model.tᵢ₋₁
    model.ω₀ = model.ωᵢ₋₁
    model.Ω₀ = model.Ωᵢ₋₁
    model.V₃₀ = model.V₃ᵢ₋₁
    model.ι₀ = model.ιᵢ₋₁

end

# function Base.getproperty(m::BurstTimingModel, s::Symbol)
#     s_str = string(s)
#     s_str = replace(s_str, "w" => "ω")
#     s_str = replace(s_str, "W" => "Ω")
#     s_str = replace(s_str, "0" => "₀")
#     s_str = replace(s_str, "3" => "₃")
#     s_str = replace(s_str, "1" => "₁")
#     s_str = replace(s_str, "2" => "₂")
#     return getfield(m, Symbol(s_str))
# end

function Base.show(io::IO, m::BurstTimingModel)
    props_to_show = (:m₁₂, :η, :e₀, :p₀, :ω₀, :Ω₀, :ι₀, :m₃, :ω₃, :p₃, :e₃, :V₃₀, :Ω₃, :ι₃)
    for prop in props_to_show
        println(io, prop, " ", getproperty(m, prop))
    end

    return nothing
end

function peak_f_GW(M, e, p)
    return sqrt(Constants.G*M)/π*(1 + e)^(1.195)/sqrt(p^3)
end

peak_f_GW(m::BurstTimingModel, idx=1) = peak_f_GW(m.m₁₂ * Constants.Msol, get_e(m, idx), get_p(m, idx))

end # module EccentricBurstTiming
