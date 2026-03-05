module Utils

include("constants.jl")

using Distributions: Normal, pdf

const a₀₁ = 14.17740665941005
const a₀₂ = -0.2369034013387491
const a₀₃ = 0.9624394843324173
const a₀₄ = -2.415911957053192
const a₁ = -16.587168717092286

const b₀ = 170π / 3
const b₁ = -139.376624104913
const b₂ = -1.088644959382641


function get_A₀(eta, pprev, m = 1)
    return a₀₁ + a₀₂ * (4eta)^a₀₃ * (0.1pprev / m)^a₀₄
end

function get_Acoeff(eta, pprev, eprev, m = 1)
    return get_A₀(eta, pprev, m) + a₁ * (1 - eprev^2)
end

function get_Bcoeff(eprev)
    return b₀ + b₁ * (1 - eprev^2) + b₂ * (1 - eprev^2)^2
end

function convolve!(w, u, v)
    if length(u) < length(v)
        return convolve!(w, v, u)
    end

    n = length(u)
    m = length(v)
    for i in 1:(n + m - 1)
        s = zero(eltype(u))
        @simd for j in max(i - m, 0):(min(n, i) - 1)
            s += u[j + 1] * v[i - j]
        end
        w[i] = s
    end
    return
end

function sine_gaussian_chirplet_frequency_domain(f, μ, σ)
    # A frequency domain signal comprised of sine-Gaussian bursts
    # is a sum over the convolution of the Fourier transform of a
    # Sine wave and the Fourier transform of a Gaussian
    fourier_domain_gaussian = 1 / √(2π) * exp(-(2π * σ * f)^2 / (2 * (2 + 1im * μ * π * f)))
    index = findfirst(x -> (x > 2π * (1 / μ)), f)
    fourier_domain_sine = 1im * √(π / 2) * (unit_impulse(length(f), -index) - unit_impulse(length(f), index))
    conv = Vector{Float64}(undef, max(length(fourier_domain_gaussian), length(fourier_domain_sine)))
    convolve!(conv, fourier_domain_gaussian, fourier_domain_sine)
    return conv
end

function sine_gaussian_chirplet(x, mu, sigma)
    # Gaussian
    # gaussian = norm.pdf(x, mu, sigma)
    gaussian = pdf(Normal(mu, sigma), x)
    # Sine wave should have one full cycle visible within Gaussian
    # The period of a sin wave is 2pi
    if iszero(mu)
        sine_wave = sin.(x)
    else
        sine_wave = @. sin(2π / mu * x)
    end

    chirplet = @. gaussian * sine_wave

    max_val = findmax(x -> abs(x), chirplet)
    if max_val > zero(max_val)
        return chirplet / max_val
    else
        return chirplet
    end
end


function negative_gaussian_chirplet(x, mu, sigma)
    # Gaussian
    gaussian = pdf(Normal(mu, sigma), x)
    chirplet = gaussian
    max_val = findmax(x -> abs(x), chirplet)

    if max_val > zero(max_val)
        chirplet_normalised = chirplet / max_val
    else
        chirplet_normalised = chirplet
    end

    return chirplet_normalised
end

function burst_signal(t, burst_times, amp, M)
    burst_list_cross = [sine_gaussian_chirplet(t, t_b, 0.00015 * M) for t_b in burst_times]
    burst_list_plus = [negative_gaussian_chirplet(t, t_b, 0.00015 * M) for t_b in burst_times]
    plus = @. amp * sum(burst_list_plus)
    cross = @. amp * sum(burst_list_cross)
    # return {'plus': plus, 'cross': cross}
    return Dict("plus" => plus, "cross" => cross)
end

end
