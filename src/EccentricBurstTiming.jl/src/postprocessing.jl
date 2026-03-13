get_e_array(model, N) = model.e[(end - min(N, length(model.e) - 1)):end]
get_p_array(model, N) = model.p[(end - min(N, length(model.p) - 1)):end] .* model.m₁₂ * Constants.Mconvert #.|> x -> geometric_to_SI_length(model.m₁₂, x)#.* Constants.Mconvert
get_t_array(model, N) = model.t[(end - min(N, length(model.t) - 1)):end] .* model.m₁₂ * Constants.Msolsec #.|> x -> geometric_to_SI_time(model.m₁₂, x)#.* Constants.Msolsec
get_w_array(model, N) = model.ω[(end - min(N, length(model.ω) - 1)):end]
get_W_array(model, N) = model.Ω[(end - min(N, length(model.Ω) - 1)):end]
get_i_array(model, N) = model.ι[(end - min(N, length(model.ι) - 1)):end]
get_V3_array(model, N) = model.V₃[(end - min(N, length(model.V₃) - 1)):end]
get_a_array(model, N) = get_p_array(model, N) ./ (1 .- get_e_array(model, N) .^ 2)

get_e(model, idx) = model.e[idx]
get_p(model, idx) = model.p[idx] * model.m₁₂ * Constants.Mconvert
get_t(model, idx) = model.t[idx] * model.m₁₂ * Constants.Msolsec
get_w(model, idx) = model.w[idx]
get_W(model, idx) = model.Ω[idx]
get_i(model, idx) = model.i[idx]
get_V3(model, idx) = model.V3[idx]
get_a(model, idx) = get_p(model, idx) / (1 .- get_e(model, idx) ^ 2)

to_meters(model, x::Number) = x * model.m₁₂ * Constants.Mconvert
to_seconds(model, x::Number) = x * model.m₁₂ * Constants.Msolsec

to_meters(model, x::AbstractVector) = x .* (model.m₁₂ * Constants.Mconvert)
to_seconds(model, x::AbstractVector) = x .* (model.m₁₂ * Constants.Msolsec)

function get_arrays(model, N)
    e_array = get_e_array(model, N)
    p_array = get_p_array(model, N)
    t_array = get_t_array(model, N)
    w_array = get_w_array(model, N)
    V3_array = get_V3_array(model, N)

    W_array = get_W_array(model, N)
    i_array = get_i_array(model, N)

    return e_array, p_array, t_array, w_array, V3_array, W_array, i_array
end

function line_up_burst_times(t_array_perturbed, t_array_unperturbed)
    time_shift = t_array_perturbed[end] - t_array_unperturbed[end]
    # t_array_unperturbed .+= time_shift

    # return t_array_perturbed, t_array_unperturbed .+ time_shift
    return time_shift
end

function get_burst_signal(model, sampling_frequency; A = 1)
    return t = range(-maximum(t_array_perturbed) - 0.05, maximum(t_array_perturbed) + 0.05, step = 2 * (maximum(t_array_perturbed) + 0.05) / sampling_frequency)
end
