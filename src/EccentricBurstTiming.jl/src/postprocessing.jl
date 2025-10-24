

get_e_array(model, N) = model.e[end-min(N, length(model.e)-1):end]
get_p_array(model, N) = model.p[end-min(N, length(model.p)-1):end] .* model.m₁₂*Constants.Mconvert
get_t_array(model, N) = model.t[end-min(N, length(model.t)-1):end] .* model.m₁₂*Constants.Msolsec
get_w_array(model, N) = model.w[end-min(N, length(model.w)-1):end]
get_W_array(model, N) = model.Ω[end-min(N, length(model.Ω)-1):end]
get_i_array(model, N) = model.i[end-min(N, length(model.i)-1):end]
get_V3_array(model, N) = model.V3[end-min(N, length(model.V3)-1):end]

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

    return t_array_perturbed, t_array_unperturbed .+ time_shift
end 

function get_burst_signal(model, sampling_frequency; A=1)
    t = range(-maximum(t_array_perturbed)-0.05, maximum(t_array_perturbed)+0.05, step=2 * (maximum(t_array_perturbed) + 0.05) / sampling_frequency)


end