import FNTexture
import arraymancer, complex, math, random, sequtils
import FNLinalg
import tables
import std/strutils

const I = complex64(0,1)
const sqrt2 = 1.414214
const vH = 246  # GeV

const mass_factor* = vH / sqrt2

# C vector specification:
# C[0-8] = re_Cu
# C[9-17] = im_Cu
# C[18-26] = re_Cd
# C[27-35] = im_Cd
# C[36] = log_eps

type UDHierarchy* = enum
    normal, inverted, degen

# NOTE: mass hierarchy enum only affects comparisons to SM, rather
# than inverting the return values themselves.  (So, it does nothing in this
# module!)
type FNTheory* = object 
    tex*: FNTexture
    C*: array[37, float]
    hierarchy*: UDHierarchy


#[
    # Work in progress...
proc `$`(self: FNTheory): string = 
    let Cu, Cd = self.Cmat()
    result = "FNTheory("
    result &= """{
    {$#, $#, $#},
    {$#, $#, $#},
    {$#, $#, $#}}""" % [ Cu[0,0], Cu[0,1], Cu[0,2],
                         Cu[1,0], Cu[1,1], Cu[1,2],
                         Cu[2,0], Cu[2,1], Cu[2,2] ]

    result &= ")"
]#

type TheoryPriorType* = enum
    uniform, lognormal, wide_lognormal

proc makeRandomTheory*(tex: FNTexture, prior_type: TheoryPriorType = uniform,
    ud_hierarchy: UDHierarchy = normal): FNTheory =

    var C: array[37, float]

    case prior_type
    of uniform:
        # [-3, 3] uniform dist
        for i in 0 ..< 36:
            C[i] = 6.0*(rand(1.0) - 0.5)
    of lognormal:
        var phase, mag: array[27, float]

        for i in concat(toSeq(0 ..< 9), toSeq(18 ..< 27)):
            phase[i] = 2 * 3.1415926535 * rand(1.0)
            mag[i] = gauss(0, 0.691)

            C[i] = exp(mag[i]) * cos(phase[i])
            C[i+9] = exp(mag[i]) * sin(phase[i])
    of wide_lognormal:
        var phase, mag: array[27, float]

        for i in concat(toSeq(0 ..< 9), toSeq(18 ..< 27)):
            phase[i] = 2 * 3.1415926535 * rand(1.0)
            mag[i] = gauss(0, 1.382)

            C[i] = exp(mag[i]) * cos(phase[i])
            C[i+9] = exp(mag[i]) * sin(phase[i])

    
    C[36] = ln(0.165959)
    result = FNTheory(tex: tex, C: C, hierarchy: ud_hierarchy)


proc Cmat*(thy: FNTheory): (Tensor[Complex64], Tensor[Complex64]) =
    let Cu_re = thy.C[0..<9].toTensor()
    let Cu_im = thy.C[9..<18].toTensor()
    let Cd_re = thy.C[18..<27].toTensor()
    let Cd_im = thy.C[27..<36].toTensor()

    var Cu = Cu_re.asType(Complex[float64]) +  I * Cu_im.asType(Complex[float64])
    var Cd = Cd_re.asType(Complex[float64]) +  I * Cd_im.asType(Complex[float64])

    Cu = Cu.reshape(3,3)
    Cd = Cd.reshape(3,3)

    return (Cu, Cd)


proc Yukawa*(thy: FNTheory): (Tensor[Complex64], Tensor[Complex64]) =
    let (Cu, Cd) = thy.Cmat()

    var eps = exp(thy.C[36])

    var Yu : array[3, array[3, Complex64]]
    var Yd : array[3, array[3, Complex64]]

    var eps_u: Complex64
    var eps_d: Complex64

    for i in 0..2:
        for j in 0..2:
            eps_u = complex(pow(eps,abs(thy.tex.XQ[i] - thy.tex.Xu[j]).float), 0)
            eps_d = complex(pow(eps,abs(thy.tex.XQ[i] - thy.tex.Xd[j]).float), 0)
            Yu[i][j] = Cu[i,j] * eps_u
            Yd[i][j] = Cd[i,j] * eps_d

    return (Yu.toTensor(), Yd.toTensor())

proc conjugate(a: Tensor[Complex64]): Tensor[Complex64] =
    var flat = a.toFlatSeq()

    var flat_star: seq[Complex64] 

    for x in flat:
        flat_star.add(conjugate(x))

    result = flat_star.toTensor().reshape(a.shape)

proc JarlskogInv(VCKM: Tensor[Complex64]): float64 =
    var J_tmp = VCKM[0,1]
    J_tmp *= VCKM[1,2]
    J_tmp *= conjugate(VCKM[0,2])
    J_tmp *= conjugate(VCKM[1,1])

    result = abs(J_tmp.im)

proc getMassMatrix*(Y: Tensor[Complex64]): (array[3, float64], array[9, Complex64]) =
    var Y_star = conjugate(Y)
    let Y_sq = einsum(Y, Y_star):
        Y[i,j] * Y_star[k,j]
    
    var Ysq_array: array[9, Complex[float64]]

    # LAPACK interface seems to be transposing matrix
    # compared to the shape here...
    for i, x in Y_sq.transpose().toFlatSeq():
        Ysq_array[i] = x

    result = eigh3(Ysq_array)


proc getSMParams*(thy: FNTheory): array[25, float64] =
    var Yu, Yd: Tensor[Complex64]
    (Yu, Yd) = thy.Yukawa()

    var Ud, Uu: array[9, Complex64]
    var Md_sq, Mu_sq: array[3, float64]

    (Mu_sq, Uu) = getMassMatrix(Yu)
    (Md_sq, Ud) = getMassMatrix(Yd)

    # Calculate CKM matrix
    # Eigenvectors returned by LAPACK are stored in ROWS,
    # so transpose yet again!
    var Ud_tens = Ud.toTensor().reshape(3,3).transpose()
    var Uu_tens = Uu.toTensor().reshape(3,3).transpose()

    var Uu_dagger = conjugate(Uu_tens).transpose()

    let VCKM = einsum(Uu_dagger, Ud_tens):
        Uu_dagger[i,j] * Ud_tens[j,k]

    # Convert Yukawa couplings to quark masses
    var Mu, Md: array[3, float]

    for i, M2 in Mu_sq:
        Mu[i] = mass_factor * sqrt(M2)
    for i, M2 in Md_sq:
        Md[i] = mass_factor * sqrt(M2)

    # Package results for return
    let VCKM_flat = VCKM.toFlatSeq()

    var VCKM_re, VCKM_im : array[9, float64]

    for i, V in VCKM_flat:
        VCKM_re[i] = V.re
        VCKM_im[i] = V.im

    result[0 .. 2] = Mu
    result[3 .. 5] = Md
    result[6 .. 14] = VCKM_re
    result[15 .. 23] = VCKM_im
    result[24] = JarlskogInv(VCKM)


proc reportMixing*(thy: FNTheory): array[8, float64] =
    var Yu, Yd: Tensor[Complex64]
    (Yu, Yd) = thy.Yukawa()

    var Ud, Uu: array[9, Complex64]
    var Md_sq, Mu_sq: array[3, float64]

    (Mu_sq, Uu) = getMassMatrix(Yu)
    (Md_sq, Ud) = getMassMatrix(Yd)

    # Calculate CKM matrix
    # Eigenvectors returned by LAPACK are stored in ROWS,
    # so transpose yet again!
    var Ud_tens = Ud.toTensor().reshape(3,3).transpose()
    var Uu_tens = Uu.toTensor().reshape(3,3).transpose()

    let tr_Ud = abs(Ud_tens[0,0]) + abs(Ud_tens[1,1]) + abs(Ud_tens[2,2])
    let tr_Uu = abs(Uu_tens[0,0]) + abs(Uu_tens[1,1]) + abs(Uu_tens[2,2])

    result[0] = tr_Ud
    result[1] = tr_Uu
    result[2] = abs(Ud_tens[0,1]) / tr_Ud
    result[3] = abs(Ud_tens[0,2]) / tr_Ud
    result[4] = abs(Ud_tens[1,2]) / tr_Ud
    result[5] = abs(Uu_tens[0,1]) / tr_Uu
    result[6] = abs(Uu_tens[0,2]) / tr_Uu
    result[7] = abs(Uu_tens[1,2]) / tr_Uu

proc reportRotationMatrices*(thy: FNTheory): 
    (Tensor[Complex64], Tensor[Complex64], Tensor[Complex64], Tensor[Complex64]) =
    var Yu, Yd, Uu, Ud, Kd, Ku: Tensor[Complex64]
    (Yu, Yd) = thy.Yukawa()

    var Ud_raw, Uu_raw, Kd_raw, Ku_raw: array[9, Complex64]
    var Md_sq, Mu_sq, Mu_sq_2, Md_sq_2: array[3, float64]

    (Mu_sq, Uu_raw) = getMassMatrix(Yu)
    (Md_sq, Ud_raw) = getMassMatrix(Yd)

    Ud = Ud_raw.toTensor().reshape(3,3).transpose()
    Uu = Uu_raw.toTensor().reshape(3,3).transpose()

    # Find the Ku/Kd right-handed rotation matrices
    # From the definition Yu = Uu Mu Ku^\dagger, we have:
    # Yu^\dagger Yu = Ku Mu Uu^\dagger Uu Mu Ku^\dagger
    # Yu^\dagger = Ku Mu^2 Ku^\dagger
    # i.e. we can use the same routine we're using to get Uu
    # to calculate Ku, except we pass in Yu^\dagger.

    # (Mu_sq_2, Ku_raw) = getMassMatrix(conjugate(Yu).transpose())
    # (Md_sq_2, Kd_raw) = getMassMatrix(conjugate(Yd).transpose())

    # Kd = Kd_raw.toTensor().reshape(3,3).transpose()
    # Ku = Ku_raw.toTensor().reshape(3,3).transpose()

    # This should have worked, but mysteriously failed...  
    # Another approach: K and U are unitary, inverse is trivial;
    # M is diagonal, inverse is easy.  So instead, find K as follows:
    # Yu^\dagger = Ku Mu Uu^\dagger, so
    # Ku = Yu^\dagger Uu Mu^{-1}.

    let Mu_inv = [
        [1/sqrt(Mu_sq[0]), 0.0, 0.0],
        [0.0, 1/sqrt(Mu_sq[1]), 0.0],
        [0.0, 0.0, 1/sqrt(Mu_sq[2])]
    ].toTensor().asType(Complex[float64])
    
    let Md_inv = [
        [1/sqrt(Md_sq[0]), 0.0, 0.0],
        [0.0, 1/sqrt(Md_sq[1]), 0.0],
        [0.0, 0.0, 1/sqrt(Md_sq[2])]
    ].toTensor().asType(Complex[float64])

    Ku = Yu.conjugate().transpose() * Uu * Mu_inv
    Kd = Yd.conjugate().transpose() * Ud * Md_inv

    result = (Uu, Ud, Ku, Kd)




proc reportSMParams*(thy: FNTheory): TableRef[string, seq[float]] =
    let SM_params = getSMParams(thy)

    result = newTable[string, seq[float]]()

    result["Mu"] = SM_params[0..2]
    result["Md"] = SM_params[3..5]

    let J = SM_params[24]
    let V12 = sqrt(SM_params[7]^2 + SM_params[16]^2)
    let V13 = sqrt(SM_params[8]^2 + SM_params[17]^2)
    let V23 = sqrt(SM_params[11]^2 + SM_params[20]^2)

    result["VCKM_slice"] = @[V12, V23, V13, J]

proc getRawParamString*(thy: FNTheory): string =
    for i in 0 ..< 36:
        result &= $thy.C[i] & ","

    result &= $thy.C[^1]

proc getSMParamString*(thy: FNTheory): string =
    let SM_report = reportSMParams(thy)

    for M in SM_report["Mu"]:
        result &= $M & ","
    for M in SM_report["Md"]:
        result &= $M & ","
    for V in SM_report["VCKM_slice"]:
        result &= $V & ","
    
    result = result[0.. ^2]

proc getFullSMParamString*(thy: FNTheory): string = 
    let SM_report = getSMParams(thy)

    for x in SM_report:
        result &= $x & ","

    result = result[0.. ^2]

proc getRotationString*(thy: FNTheory): string =
    let (Uu, Ud, Ku, Kd) = reportRotationMatrices(thy)

    result &= "Uu,"
    for x in Uu.toFlatSeq():
        result &= $x.re & "," & $x.im & ","
    result = result[0.. ^2]
    result &= "\n"

    result &= "Ud,"
    for x in Ud.toFlatSeq():
        result &= $x.re & "," & $x.im & ","
    result = result[0.. ^2]
    result &= "\n"

    result &= "Ku,"
    for x in Ku.toFlatSeq():
        result &= $x.re & "," & $x.im & ","
    result = result[0.. ^2]
    result &= "\n"

    result &= "Kd,"
    for x in Kd.toFlatSeq():
        result &= $x.re & "," & $x.im & ","
    result = result[0.. ^2]
#    result &= "\n"

#    result = result[0.. ^3]




