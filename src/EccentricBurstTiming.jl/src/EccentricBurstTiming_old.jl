module EccentricBurstTiming

export BurstTimingModel, evolve!, iterate, get_arrays, line_up_burst_times, get_e_array, get_p_array, get_t_array, get_w_array, get_V3_array

include("constants.jl")
include("utilities.jl")
include("postprocessing.jl")
mutable struct BurstTimingModel{T, vecT}
    lambda::T
    iota::T
    M::T # M
    eta::T
    eta3::T
    m3::T
    R3::T
    e::vecT
    p::vecT
    t::vecT
    w::vecT
    V3::vecT
    eprev::T
    pprev::T
    tprev::T
    wprev::T
    V3prev::T
    Omega3::T
    e0::T
    p0::T
    t0::T
    w0::T
    V30::T

    function BurstTimingModel(;
            e0 = 0.99, p0 = 30, t0 = 0, M = 60, eta = 0.2,
            w0 = 0, m3 = 1.0e7, R3 = 1.1e7, V3 = π / 3,
            lambda = 0, iota = 0
        )
        T1 = typeof(e0)
        T2 = typeof(T1[])

        return new{T1, T2}(
            lambda, iota, M, eta, m3 / M, m3, R3,
            T1[e0], T1[p0], T1[t0], T1[w0], T1[V3],
            e0, p0, t0, w0, V3,
            √((M + m3) / R3^3),
            e0, p0, t0, w0, V3
        )
    end
end

function get_e_next(model::BurstTimingModel)
    two_body_update = (608 * π / 15 * model.eta * (model.M / model.pprev)^(5 / 2)) * (1 + 121 / 304 * model.eprev^2)

    three_body_perturbation = if model.m3 > zero(model.m3)
        15π / 2 * model.eta3 * (model.M / model.R3)^3 / (model.M / model.pprev)^3 * sin(2 * (model.V3prev - model.wprev)) / (1 - model.eprev^2)^(5 / 2)
    else
        0.0
    end

    return model.eprev * (1 - two_body_update - three_body_perturbation)
end

function get_p_next(model)
    two_body_update = (128π / 5 * model.eta * (model.M / model.pprev)^(5 / 2)) * (1 + 7 / 8 * model.eprev^2)

    three_body_perturbation = if model.m3 > zero(model.m3)
        15π * model.eta3 * (model.M / model.R3)^3 / (model.M / model.pprev)^3 * model.eprev^2 * sin(2 * (model.V3prev - model.wprev)) / (1 - model.eprev^2)^(7 / 2)
    else
        0.0
    end

    return model.pprev * (1 - two_body_update + three_body_perturbation)
end

function get_t_next(model)
    Acoeff = Utils.get_Acoeff(model.eta, model.pprev, model.M, model.eprev)
    Bcoeff = Utils.get_Bcoeff(model.eprev)
    # println(Acoeff, " ", Bcoeff)
    two_body_update = 2π / √model.M * (
        (model.pprev + model.eta * (model.M^(5 / 2) / model.pprev^(3 / 2)) * Acoeff) /
            (1 - model.eprev^2 + model.eta * (model.M / model.pprev)^(5 / 2) * Bcoeff)
    )^(3 / 2)

    three_body_perturbation = if model.m3 > zero(model.m3)
        1 + model.eta3 * (model.M / model.R3)^3 / (model.M / model.pprev)^3 * (5 * (4 + (3 * model.eprev^2)) + (96 + (51 * model.eprev^2)) * cos(2 * (model.V3prev - model.wprev))) / (16 * (1 - model.eprev^2)^3)
    else
        1.0
    end

    return model.tprev + (two_body_update * three_body_perturbation)
end

function get_w_next(model)
    three_body_perturbation = if model.m3 > zero(model.m3)
        3π / 2 * model.eta3 * (model.M / model.R3)^3 / (model.M / model.pprev)^3 * (1 + (5cos(2 * (model.V3prev - model.wprev)))) / (1 - model.eprev^2)^(5 / 2)
    else
        0
    end

    return model.wprev + three_body_perturbation
end

function get_V3_next(model)
    return model.V3prev + (model.Omega3 * (model.t[end] - model.tprev))
end

function iterate!(model::BurstTimingModel)

    e_next = get_e_next(model)
    push!(model.e, e_next)

    p_next = get_p_next(model)
    push!(model.p, p_next)

    t_next = get_t_next(model)
    push!(model.t, t_next)

    w_next = get_w_next(model)
    push!(model.w, w_next)

    V3_next = get_V3_next(model)
    push!(model.V3, V3_next)

    model.eprev = e_next
    model.pprev = p_next
    model.tprev = t_next
    model.wprev = w_next
    model.V3prev = V3_next
    return nothing
end

function observed_burst_time_offsets_due_to_com_motion(model)
    com_motion_induced_offset = @. model.R3 / model.M * sin(model.iota) * cos(model.V3 + model.lambda)
    return com_motion_induced_offset
end

function evolve!(model::BurstTimingModel, n_bursts)
    n = 0
    while (n <= n_bursts) && (model.eprev > 0.7) && (model.eprev < 1.0) && (model.pprev > 0)
        if model.pprev > 2 * (3 + model.eprev)
            #if model.m3 == 0.0 || 0.1 > (np.power(
            #        (1 + model.eprev) ** 9 * model.eta3 /
            #        (1 - model.eprev) ** 2 /
            #        ((1 + model.eprev) * np.power(model.M / model.eprev, 0.5)) ** 11,
            #         1/3)) / model.R3:
            if iszero(model.m3) || 0.1 > (cbrt((1 + model.eprev)^(13 / 2) * model.eta3 / (1 - model.eprev)^2 / model.eta / ((1 + model.eprev) * √(model.M / model.eprev))^11)) / model.R3
                iterate!(model)
            end
            n += 1
        else
            n = Nbursts + 1
        end
    end
    # Cut off last point to ensure within region of validity
    model.e = model.e[1:(end - 1)]
    model.p = model.p[1:(end - 1)]
    model.t = model.t[1:(end - 1)]
    model.w = model.w[1:(end - 1)]
    model.V3 = model.V3[1:(end - 1)]
    # Add observed time offset due to system inclination
    return if model.m3 > 0
        additional_offset = observed_burst_time_offsets_due_to_com_motion(model)
        model.t .-= additional_offset
    end
    # Normalise so that first burst happens at time 0
    #model.t = [x - model.t[0] for x in model.t]
end


function main(;
        e0 = 0.99, p0 = 30, t0 = 0, M = 60, eta = 0.2,
        w0 = 0, m3 = 1.0e7, R3 = 1.1e7, V3 = π / 3,
        N = 73, A = 1.0e-21, geocent_time = 0,
        ra = 0, dec = 0, psi = 0, lambda = 0, iota = 0
    )

    burst_model_perturbed = BurstTimingModel(
        e0 = e0, p0 = p0, t0 = 0, w0 = w0, M = 1, eta = eta,
        m3 = m3, R3 = R3, V3 = V30, lambda = lambda, iota = iota
    )
    evolve!(burst_model_perturbed, N_request)

    burst_model_perturbed_tilted = BurstTimingModel(
        e0 = e0, p0 = p0, t0 = 0, w0 = w0, M = 1, eta = eta,
        m3 = m3, R3 = R3, V3 = V30, lambda = np.pi / 6, iota = np.pi / 4
    )
    evolve!(burst_model_perturbed_tilted, N_request)

    burst_model_unperturbed = BurstTimingModel(
        e0 = e0, p0 = p0, t0 = 0, w0 = w0, M = 1, eta = eta,
        m3 = 0, R3 = R3, V3 = V30, lambda = lambda, iota = iota
    )
    evolve!(burst_model_unperturbed, N_request)

    e_array_perturbed = burst_model_perturbed.e[(end - N):end]
    p_array_perturbed = burst_model_perturbed.p[(end - N):end] * M * Constants.Mconvert / 1.0e-5 # [au / 10^-5]
    t_array_perturbed = burst_model_perturbed.t[(end - N):end] * M * Constants.Msolsec
    w_array_perturbed = burst_model_perturbed.w[(end - N):end]
    V3_array_perturbed = burst_model_perturbed.V3[(end - N):end]

    e_array_perturbed_tilted = burst_model_perturbed_tilted.e[(end - N):end]
    p_array_perturbed_tilted = burst_model_perturbed_tilted.p[(end - N):end] * M * Constants.Mconvert / 1.0e-5 # [au / 10^-5]
    t_array_perturbed_tilted = burst_model_perturbed_tilted.t[(end - N):end] * M * Constants.Msolsec
    w_array_perturbed_tilted = burst_model_perturbed_tilted.w[(end - N):end]
    V3_array_perturbed_tilted = burst_model_perturbed_tilted.V3[(end - N):end]

    e_array_unperturbed = burst_model_unperturbed.e[(end - N):end]
    p_array_unperturbed = burst_model_unperturbed.p[(end - N):end] * M * Constants.Mconvert / 1.0e-5 # [au / 10^-5]
    t_array_unperturbed = burst_model_unperturbed.t[(end - N):end] * M * Constants.Msolsec
    w_array_unperturbed = burst_model_unperturbed.w[(end - N):end]

    println(t_array_unperturbed[end] - t_array_unperturbed[end - 1])

    println(t_array_unperturbed)
    println(t_array_perturbed)

    # Make final burst times line up
    time_shift = t_array_perturbed[end] - t_array_unperturbed[end]
    t_array_unperturbed += time_shift
    time_shift_2 = t_array_perturbed[end] - t_array_perturbed_tilted[end]
    t_array_perturbed_tilted += time_shift_2
    println(length(t_array_unperturbed))
    println(length(t_array_perturbed))
    println(length(t_array_perturbed_tilted))

    sampling_frequency = 512
    # t = np.arange(-np.max(t_array_perturbed)-0.05, np.max(t_array_perturbed)+0.05, 2 * (np.max(t_array_perturbed) + 0.05) / sampling_frequency)
    t = range(-maximum(t_array_perturbed) - 0.05, stop = maximum(t_array_perturbed) + 0.05, step = 2 * (maximum(t_array_perturbed) + 0.05) / sampling_frequency)
    A = 1

    burst_signal_perturbed = Utils.burst_signal(t, t_array_perturbed, A, M)
    return burst_signal_unperturbed = Utils.burst_signal(t, t_array_unperturbed, A, M)
end

function test() end

end # module EccentricBurstTiming
