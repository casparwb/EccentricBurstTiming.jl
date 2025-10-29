
module Constants

using Unitful, UnitfulAstro

const a01 = 14.17740665941005
const a02 = -0.2369034013387491
const a03 = 0.9624394843324173
const a04 = -2.415911957053192
const a1 = -16.587168717092286

const b0 = 170π/3
const b1 = -139.376624104913
const b2 = -1.088644959382641

const G = 6.67408e-11u"m^3*kg^-1*s^-1"
const c = 299792458.0u"Msun"
const sqrt_G = sqrt(G)
const Msol = 1.998e30
const au = 149597870700.0
const pc_au = 206265.0
const pc = pc_au * au


const Msolsec = G*Msol/c^3
const Mconvert = Msolsec*c/au

end