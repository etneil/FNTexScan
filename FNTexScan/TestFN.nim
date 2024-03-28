import unittest
import tables
import arraymancer
import FNTexture
import FNTheory
import FNCost
import random

test "create valid FNTexture":
    let testTex = initFNTexture(
        XQ = [3,3,0],
        Xu = [-2,-1,0],
        Xd = [-2,-2,-1])
    echo testTex
    check testTex.isValid == true

test "count textures":
    let testMaxQ = 2
    let allValid = getAllValidTextures(testMaxQ)

    echo len(allValid)

    echo allValid[0]
    echo allValid[^1]
#    echo allValid[160000]

test "create invalid FNTexture":
    let testTex2 = initFNTexture(
        XQ = [3,3,1],
        Xu = [-2,-1,0],
        Xd = [-2,-2,-1])
    echo testTex2
    check testTex2.isValid == false

test "create random FNTheory":
    let testTex = initFNTexture(
        XQ = [3,3,0],
        Xu = [-2,-1,0],
        Xd = [-2,-2,-1])
    let testThy = makeRandomTheory(testTex, prior_type=TheoryPriorType.uniform)
    echo testThy
    
    let testThy2 = makeRandomTheory(testTex, prior_type=TheoryPriorType.lognormal)
    echo testThy2

test "fixed FNTheory tests":
    let testTex = initFNTexture(
        XQ = [3,2,0],
        Xu = [-2,-4,0],
        Xd = [-3,-3,-4])
    var test_C: array[37, float]
    for i in 0..<36:
        test_C[i] = (i.float+1) / 36
    test_C[36] = 0.2

    let testThy = FNTheory(tex: testTex, C: test_C)

    echo "Yukawas: ", testThy.Yukawa()

    echo testThy.C
    echo "getSM: ", getSMParams(testThy)
    echo reportSMParams(testThy)

test "check Fedele fixed model":
    # arXiv:2009.05587, appendix A
    let fTex = initFNTexture(
        XQ = [1, 1, 0],
        Xu = [-2, -2, 0],
        Xd = [-2, -2, -2]
    )

    # CU_re, CU_im, CD_re, CD_im, eps
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
    echo "Fedele SM params: ", reportSMParams(fThy)

    echo "Fedele SM score, Xi = ", score_SM_deviations(fThy, max_exp = true)




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




