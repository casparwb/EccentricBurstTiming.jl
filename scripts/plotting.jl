using CairoMakie

function example_run(;N=70, N_plot=20, args...)

    m3 = get(args, :m3, 1e4)

    m = BurstTimingModel(m3=0, w0=0)
    m_p = BurstTimingModel(m3=m3, w0=0, i0=0);

    evolve!(m, N)
    evolve!(m_p, N)
    plot_e_p_evolution(m_p, m, N_plot)
end

function plot_e_p_evolution(model_perturbed, model_unperturbed, N)

    e_perturbed = get_e_array(model_perturbed, N)
    p_perturbed = get_p_array(model_perturbed, N) .* EccentricBurstTiming.Constants.m_to_Rsun
    t_perturbed = get_t_array(model_perturbed, N)

    e_unperturbed = get_e_array(model_unperturbed, N)
    p_unperturbed = get_p_array(model_unperturbed, N) .* EccentricBurstTiming.Constants.m_to_Rsun
    t_unperturbed = get_t_array(model_unperturbed, N)

    # t_perturbed ./= sqrt(model_perturbed.m3*model_perturbed.p3)
    # t_perturbed ./= sqrt(EccentricBurstTiming.Constants.c/EccentricBurstTiming.Constants.G)

    # t_perturbed, t_unperturbed = line_up_burst_times(t_perturbed, t_unperturbed)
    time_shift = line_up_burst_times(t_perturbed, t_unperturbed)
    # return t_perturbed, t_unperturbed
    # t_perturbed = let t = t_perturbed
    #     t = t .- t[1]
    #     t ./ t[end]
    # end

    # t_unperturbed = let t = t_unperturbed
    #     t = t .- t[1]
    #     t ./ t[end]
    # end



    # aout = model_perturbed.R3*model_perturbed.M*Constants.Mconvert
    fig = Figure(size=(600, 800))
    ax_e = Axis(fig[1, 1], ylabel="Eccentricity", yticks=LinearTicks(7))
    ax_p = Axis(fig[2, 1], xlabel="Time [s]", ylabel=L"Semi-latus rectum [R$_\odot$]")

    hidexdecorations!(ax_e)
    linkxaxes!(ax_e, ax_p)

    cs = Makie.wong_colors()
    sc = scatter!(ax_e, t_unperturbed .+ time_shift, e_unperturbed, color=cs[1], label="unperturbed", markersize=15)
    scatter!(ax_e, t_perturbed, e_perturbed, color=cs[2], label="perturbed", markersize=10)

    sc2 = scatter!(ax_p, t_unperturbed .+ time_shift, p_unperturbed, color=cs[1], markersize=15)
    scatter!(ax_p, t_perturbed, p_perturbed, color=cs[2], markersize=10)

    translate!(sc, 0, 0, -10)
    translate!(sc2, 0, 0, -10)

    axislegend(ax_e, position=:lb, labelsize=20)
    return fig
end

function get_total_shift(model_unperturbed, model_perturbed, N)
    t_perturbed = get_t_array(model_perturbed, N)
    t_unperturbed = get_t_array(model_unperturbed, N)
    t_perturbed, t_unperturbed = line_up_burst_times(t_perturbed, t_unperturbed)
    
    N_perturbed = length(t_perturbed)
    N_unperturbed = length(t_unperturbed)
    if N_perturbed > N_unperturbed
        dt = @. abs(t_perturbed[1:N_unperturbed] - t_unperturbed)
    elseif N_perturbed < N_unperturbed
        dt = @. abs(t_perturbed - t_unperturbed[1:N_perturbed])
    else
        dt = @. abs(t_perturbed - t_unperturbed)
    end
    
    return dt
end

function plot_total_shift(burst_model_unperturbed, burst_model_perturbed, N)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Burst number", ylabel="Δt [s]")

    dt = get_total_shift(burst_model_unperturbed, burst_model_perturbed, N)

    scatter!(ax, dt)
    return fig
end