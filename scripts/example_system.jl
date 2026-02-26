
const au_to_Rsun = 214.4772339784946
function make_system_from_paper()

    R = 6.55*au_to_Rsun
    e_in = 0.99
    a_in = 1.8e-5*au_to_Rsun/(1 - e_in^2)
    m12 = 60
    η = 0.2
    m3 = 6e8
    w0 = 0.0
    V3 = π/3

    m = BurstTimingModel(m12=m12, a0=a_in, R3=R, eta=η, m3=m3, w0=w0, V3=V3)
    m_up = BurstTimingModel(m12=m12, a0=a_in, R3=R, eta=η, m3=0, w0=w0, V3=V3)
    return m, m_up
end