import unittest
import tables
import arraymancer
import FNTexture
import FNTheory
import FNCost
import random
import math
import FNLinalg

test "simple eigenvalue test":
    let A = [
        [1.0, -2.0, 0.0],
        [-2.0, 1.0, 0.0],
        [0.0, 0.0, 2.0]
    ].toTensor().asType(Complex[float64])

    var evals: array[3, float64]
    var U: array[9, Complex[float64]]
    
    var A_array: array[9, Complex[float64]]

    for i, x in A.toFlatSeq():
        A_array[i] = x

    echo "eigh3:"

    (evals, U) = eigh3(A_array)

    echo evals
    echo U

    echo "eigh3r:"

    for i, x in A.toFlatSeq():
        A_array[i] = x

    (evals, U) = eigh3r(A_array)

    echo evals
    echo U

test "small eigenvalue test":
    let Lambda = [
        [1e-10, 0.0, 0.0],
        [0.0, 1e-3, 0.0],
        [0.0, 0.0, 1.0]
    ].toTensor().asType(Complex[float64])

    let R = [
        [complex(0.707107,0.0),  complex(0.0, 0.707107),   complex(-0.0, 0.0)],
        [complex(0.0,0.707107),  complex(0.707107, 0.0),  complex(0.0, -0.0)],
        [complex(0.0,0.0),   complex(-0.0, -0.0), complex(1.0, 0.0)]
    ].toTensor()

    let Rdag = conjugate(R.transpose())
    let A = Rdag * Lambda * R

    var evals: array[3, float64]
    var U: array[9, Complex[float64]]
    var A_array: array[9, Complex[float64]]

    for i, x in A.toFlatSeq():
        A_array[i] = x

    echo "eigh3:"

    (evals, U) = eigh3(A_array)

    echo evals
    echo U

    echo "eigh3r:"

    for i, x in A.toFlatSeq():
        A_array[i] = x

    (evals, U) = eigh3r(A_array)

    echo evals
    echo U

test "unstable Yukawa test":
    let Yu = [
        [ complex(4.03138179e-10, -3.44622002e-10), complex(-1.94658409e-06, 1.10900415e-06),
  complex(-3.76007814e-05, -1.09248637e-04)],
        [ complex(1.44240996e-08,1.01225597e-08), complex(-2.66018931e-06,1.00103939e-05),
  complex(-2.24398670e-03,1.86175377e-03)],
        [ complex(-3.72600328e-06,1.85063816e-06), complex(5.18083402e-04,-1.43250974e-03),
  complex(-4.88727283e-01,1.06113560e+00)]
    ].toTensor()

    let Yudag = conjugate(Yu.transpose())

    let Ysq = (Yu * Yudag).transpose()

    var evals: array[3, float64]
    var U: array[9, Complex[float64]]
    var A_array: array[9, Complex[float64]]

    for i, x in Ysq.toFlatSeq():
        A_array[i] = x

    echo "eigh3:"

    (evals, U) = eigh3(A_array)

    echo evals
    echo U

    for i, z in evals:
        echo sqrt(z) * mass_factor


    echo "eigh3r:"

    for i, x in Ysq.toFlatSeq():
        A_array[i] = x

    (evals, U) = eigh3r(A_array)

    echo evals
    echo U

    for i, z in evals:
        echo sqrt(z) * mass_factor






test "create valid FNTexture":
    let testTex = initFNTexture(
        XQ = [3,3,0],
        Xu = [-2,-1,0],
        Xd = [-2,-2,-1])
    echo testTex
    check testTex.isValid == true

test "FNTexture string conversion":
    let testTex = initFNTexture(
        XQ = [3,3,0],
        Xu = [-2,-1,0],
        Xd = [-2,-2,-1])
    
    let texStr = testTex.toCompactString()
    echo texStr

    let rebuildTex = texStr.fromCompactString()
    echo rebuildTex

    check testTex == rebuildTex

test "count textures":
    let testMaxQ = 4
    let allValid = getAllValidTextures(testMaxQ)

    echo "Number of textures with max charge = 4: ", len(allValid)
    check len(allValid) == 167125

test "create random FNTheory":
    let testTex = initFNTexture(
        XQ = [3,3,0],
        Xu = [-2,-1,0],
        Xd = [-2,-2,-1])
    let testThy = makeRandomTheory(testTex, prior_type=TheoryPriorType.uniform)
    echo testThy

    discard testThy.Yukawa()
    discard getSMParams(testThy)

test "rotation matrices":
    let testTex = initFNTexture(
        XQ = [3,3,0],
        Xu = [-2,-1,0],
        Xd = [-2,-2,-1]
    )
    let testThy = makeRandomTheory(testTex, prior_type=TheoryPriorType.lognormal)

    var Uu, Ud, Ku, Kd: Tensor[Complex64]

    (Uu, Ud, Ku, Kd) = testThy.reportRotationMatrices()

    echo "Uu: "
    echo Uu
    echo "Rotation string: "
    echo testThy.getRotationString()

    # Test unitarity
    let U_test = Uu * Uu.conjugate().transpose()
    let K_test = Ku * Ku.conjugate().transpose()

    check abs(sum(U_test) - 3) <= 1e-6
    check abs(sum(K_test) - 3) <= 1e-6    

    # Verify that product of rotations gives the Yukawa matrices back
    let (Yu, Yd) = testThy.Yukawa()
    let SMPars = testThy.getSMParams()

    let Mu = [
        [SMpars[0], 0.0, 0.0],
        [0.0, SMpars[1], 0.0],
        [0.0, 0.0, SMpars[2]]
    ].toTensor().asType(Complex[float64]) / FNTheory.mass_factor.complex64
    let Md = [
        [SMpars[3], 0.0, 0.0],
        [0.0, SMpars[4], 0.0],
        [0.0, 0.0, SMpars[5]]
    ].toTensor().asType(Complex[float64]) / FNTheory.mass_factor.complex64

    let Yu_rebuild = Uu * Mu * Ku.conjugate().transpose()
    let Yd_rebuild = Ud * Md * Kd.conjugate().transpose()

    let Yu_diff = Yu - Yu_rebuild
    let Yd_diff = Yd - Yd_rebuild

    check abs(sum(Yu_diff)) <= 1e-6
    check abs(sum(Yd_diff)) <= 1e-6

    echo "Yu diff: "
    echo Yu_diff
    echo "Yd diff: "
    echo Yd_diff

test "check Fedele fixed model":
    # arXiv:2009.05587, appendix A
    let fTex = initFNTexture(
        XQ = [1, 1, 0],
        Xu = [-2, -2, 0],
        Xd = [-2, -2, -2]
    )

    # CU_re, CU_im, CD_re, CD_im, ln(eps)
    let fC = [
        0.296996, -0.171824, -0.218137,
        -0.288732, 0.178763, 0.225587,
        -0.173489, 0.270735, -0.185740,
        -0.962979, 0.987427, 0.965404,
        0.962114, -0.987484, -0.995365,
        0.992865, -0.967206, 0.820763,
        -0.444605, -0.479893, 0.0586865,
        0.554765, 0.438668, -0.129055,
        -0.786523, -0.483583, 0.417669,
        -0.892952, -0.882641, 1.00937,
        0.830214, 0.892558, -0.983295,
        -0.605194, -0.866336, 0.895754,
        ln(0.0894772)
    ]

    let fThy = FNTheory(tex: fTex, C: fC)

    let fSMPars = reportSMParams(fThy)
    echo "Fedele SM params: ", fSMPars

    let fedeleXi =  score_SM_deviations(fThy, max_exp = true)

    echo "Fedele SM score, Xi = ", $fedeleXi

    check fedeleXi < 2.0


test "optimize random FNTheory":
#[    
    let testTex = initFNTexture(
        XQ = [3,2,0],
        Xu = [-2,-4,0],
        Xd = [-3,-3,-3])
]#    
    let testTex = initFNTexture(
        XQ = [-4,-4,0],
        Xu = [-4,-3,0],
        Xd = [-3,1,4])


    randomize()
    let testThy = makeRandomTheory(testTex, prior_type=TheoryPriorType.uniform)
    echo testThy

    echo "Random starting point: ", reportSMParams(testThy)
    echo "Initial cost: ", score_SM_deviations(testThy)
    echo "Starting eps: ", testThy.C[^1]

    var opt_eps = optimize_eps(testThy, brute_force = false, max_exp = true)
    echo "Opt eps: ", opt_eps
    check opt_eps != 0.2

    var opt_C = testThy.C
    opt_C[36] = ln(opt_eps)
    var optThy = FNTheory(tex: testTex, C: opt_C)

    echo "Optimized theory: ", reportSMParams(optThy)
    echo "Optimized cost: ", score_SM_deviations(optThy)

test "optimize to SM":
    let testTex = initFNTexture(
        XQ = [3,2,0],
        Xu = [-2,-4,0],
        Xd = [-3,-3,-3])

    randomize()
    let testThy = makeRandomTheory(testTex, prior_type=TheoryPriorType.uniform)
    echo testThy

    echo "Initial cost: ", score_SM_deviations(testThy)

    var optThy = optimize_SM(testThy)

    echo "Optimized theory: ", reportSMParams(optThy)
    echo "Optimized cost: ", score_SM_deviations(optThy)




