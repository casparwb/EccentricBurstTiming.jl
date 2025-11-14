
module Constants

export G, c, c⁻¹, Gc⁻², G⁻¹c², Msun_to_kg, kg_to_Msun, Rsun_to_m, m_to_Rsun

const a01 = 14.17740665941005
const a02 = -0.2369034013387491
const a03 = 0.9624394843324173
const a04 = -2.415911957053192
const a1 = -16.587168717092286

const b0 = 170π/3
const b1 = -139.376624104913
const b2 = -1.088644959382641

const G = 6.67408e-11#u"m^3*kg^-1*s^-1"
const c = 299792458.0#u"Msun"

# const G = 3.9413556368747037e-7 # in Msun, Rsun, s
# const c = 0.4309220324852666 # in Rsun/s
# 
const c⁻¹ = 1/c
const sqrt_G = sqrt(G)
const Msol = 1.998e30
const au = 149597870700.0
const pc_au = 206265.0
const pc = pc_au * au


const Gc⁻² = G/c^2
const G⁻¹c² = 1/Gc⁻²

const Msun_to_kg = 1.998e30
const kg_to_Msun = 1/1.998e30
const Rsun_to_m = 6.957e8
const m_to_Rsun = 1/6.957e8

const Msolsec = G*Msol/c^3
const Mconvert = Msolsec*c

end