import FNTexture
import FNTheory
import FNCost
import FNTexScan
import random, tables, suru, arraymancer

#[
const thy_data_file = "thy_data.csv"
const xi_data_file = "xi_data.csv"
const eps_data_file = "eps_data.csv"
]#

randomize()

const max_Q = 4
const test_tex_i = 351439
#const test_tex_i = 351642
const scan_prior = "lognormal"
#const scan_prior = "wide_lognormal"
#const scan_prior = "uniform"

#var all_textures = getAllValidTextures(max_Q)
#var thisTex = all_textures[test_tex_i]

#[
let thisTex = initFNTexture(
    XQ = [1, 1, 0],
    Xu = [0, -2, 0],
    Xd = [-2, -2, -2]
)
]#

#[
let thisTex = initFNTexture(
    XQ = [3, 2, 0],
    Xu = [-2, -4, 0],
    Xd = [-3, -3, -3]
)
]#

let thisTex = initFNTexture(
    XQ = [3, 2, 0],
    Xu = [-2, -3, 0],
    Xd = [-2, -3, -3]
)

#let mirrorTex = mirror(thisTex)

#[
let tex_i = all_textures.find(thisTex)
echo tex_i

let mtex_i = all_textures.find(mirrorTex)
echo mirrorTex
echo mtex_i
]#

#[
let thisTex = initFNTexture(
    XQ = [3, -2, 0],
    Xu = [5, -2, 0],
    Xd = [-1, 0, 0]
)
]#

echo "This texture is: ", thisTex

const thy_per_tex = 10000
var all_costs, all_costs_inv, all_xi, all_xi_inv, all_nat: seq[float] = @[]
var all_costs_copyinv, all_xi_copyinv: seq[float] = @[]
var num_good, num_good_inv, num_good_copyinv: int
var randThy, randThyInv, optThy, optThyInv: FNTheory
#var optSMpars: TableRef
var opt_eps: float
var all_eps, all_eps_inv: seq[float] = @[]
var num_ud_flip = 0

const do_opt_SM = true
const do_priors = false

# Setup data file names
proc make_file_name(base_name: string): string =
    result = base_name
    if do_opt_SM:
        result = result & "_SMopt"
    if do_priors:
        result = result & "_prior"
    
    result = result & ".csv"

let raw_thy_data_file = make_file_name("raw_thy_data")
let thy_data_file = make_file_name("thy_data")
let thy_full_data_file = make_file_name("thy_full_data")
let xi_data_file = make_file_name("xi_data")
let eps_data_file = make_file_name("eps_data")
let cost_data_file = make_file_name("cost_data")
let mix_data_file = make_file_name("mix_data")
let tuning_data_file = make_file_name("tuning_data")

# Clobber data files
let f1 = open(thy_data_file, fmWrite)
f1.writeLine(" ")
f1.close()
let f2 = open(xi_data_file, fmWrite)
f2.writeLine(" ")
f2.close()
let f3 = open(eps_data_file, fmWrite)
f3.writeLine(" ")
f3.close()
let f4 = open(raw_thy_data_file, fmWrite)
f4.writeLine(" ")
f4.close()
let f5 = open(cost_data_file, fmWrite)
f5.writeLine(" ")
f5.close()
let f6 = open(mix_data_file, fmWrite)
f6.writeLine(" ")
f6.close()
let f7 = open(thy_full_data_file, fmWrite)
f7.writeLine(" ")
f7.close()
let f8 = open(tuning_data_file, fmWrite)
f8.writeLine(" ")
f8.close()

let raw_thy_f = open(raw_thy_data_file, fmAppend)
let thy_f = open(thy_data_file, fmAppend)
let xi_f = open(xi_data_file, fmAppend)
let eps_f = open(eps_data_file, fmAppend)
let cost_f = open(cost_data_file, fmAppend)
let mix_f = open(mix_data_file, fmAppend)
let thy_full_f = open(thy_full_data_file, fmAppend)
let tuning_f = open(tuning_data_file, fmAppend)

for k in suru(0 ..< thy_per_tex):
    randThy = makeRandomTheory(tex=thisTex, prior_type=scan_prior)
    randThyInv = makeRandomTheory(tex=thisTex, prior_type=scan_prior, ud_hierarchy=inverted)

    try:
        if do_priors:
            optThy = randThy
            optThyInv = randThyInv
            opt_eps = exp(randThy.C[^1])

        else:
            opt_eps = optimize_eps(randThy, brute_force = false, max_exp = true)
            optThy = FNTheory(tex: randThy.tex, C: randThy.C, hierarchy: randThy.hierarchy)
            optThy.C[^1] = ln(opt_eps)

            opt_eps = optimize_eps(randThyInv, brute_force = false, max_exp = true)
            optThyInv = FNTheory(tex: randThyInv.tex, C: randThyInv.C, hierarchy: randThyInv.hierarchy)
            optThyInv.C[^1] = ln(opt_eps)

        if do_opt_SM:
            optThy = optimize_SM(optThy)
            optThyInv = optimize_SM(optThyInv)

        all_eps.add opt_eps
        all_eps_inv.add opt_eps


    except Exception as e:
        raise
#        continue

    var copyThy = FNTheory(tex: optThy.tex, C: optThy.C, hierarchy: UDHierarchy.inverted)

    #[
    if k == 0:
        let (Cu, Cd) = optThy.Cmat()
        echo optThy
        echo "Cu = ", Cu
        echo "Cd = ", Cd
        echo "eps = ", exp(optThy.C[^1])
        echo "xi = ", score_SM_deviations(optThy, max_exp=true)
    ]#

    all_costs.add score_SM_deviations(optThy, max_exp=false)
    all_xi.add score_SM_deviations(optThy, max_exp=true)

    all_costs_inv.add score_SM_deviations(optThyInv, max_exp=false)
    all_xi_inv.add score_SM_deviations(optThyInv, max_exp=true)

    all_costs_copyinv.add score_SM_deviations(copyThy, max_exp=false)
    all_xi_copyinv.add score_SM_deviations(copyThy, max_exp=true)

    all_nat.add naturalness_cost(optThy, sigma_weight=0.691)

    let optSMpars = optThy.reportSMParams()
    if all_xi[^1] <= 2 and optSMpars["Mu"][0] > optSMpars["Md"][0]:
        num_ud_flip += 1

    raw_thy_f.writeLine(optThy.getRawParamString())
    thy_f.writeLine(optThy.getSMParamString())
    xi_f.writeLine($all_xi[^1])
    cost_f.writeLine($all_costs[^1])
#    if not do_opt_SM:
    eps_f.writeLine($opt_eps)

    mix_f.writeLine($optThy.reportMixing())
    thy_full_f.writeLine(optThy.getFullSMParamString())
    tuning_f.writeLine($all_nat[^1])

raw_thy_f.close()
thy_f.close()
xi_f.close()
eps_f.close()
cost_f.close()
mix_f.close()
thy_full_f.close()
tuning_f.close()

if not do_opt_SM:
    var mean_eps = 0.0
    for eps in all_eps:
        mean_eps += eps
    mean_eps /= float64(len(all_eps))

    echo "Eps = ", mean_eps

for target_xi in @[2.0, 5.0]:
    num_good = count_below(all_xi, target_score=target_xi)
    num_good_inv = count_below(all_xi_inv, target_score=target_xi)
    num_good_copyinv = count_below(all_xi_copyinv, target_score=target_xi)

    echo "Xi < ", target_xi, " (normal): ",  num_good
    echo "Xi < ", target_xi, " (inverted): ", num_good_inv
    echo "Xi < ", target_xi, " (normal --> inv): ", num_good_copyinv

echo "Num_ud_flip: ", num_ud_flip