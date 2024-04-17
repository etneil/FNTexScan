import arraymancer
import FNTexture, FNTheory, FNCost
import random
import tables
import std/algorithm
import docopt
import docopt/dispatch
import suru
import strutils
import os

const doc = """
FNTexScan - scan and study Frogatt-Nielsen textures.
Companion code for the paper arXiv:2306.08026.

Commands:
    scan_all -  Scan all textures; create a file 
                listing all textures which are "good",
                defined by <good_frac> being within a
                factor of 5 of the SM.  (Called "F_5" in the paper.)
                This scan is *unranked* (i.e. unordered.)
    rescan -    Re-run scan of all textures in the given
                texture file, producing a new file, optionally
                restricting output to just the top <top_N>.
                This is meant to be "high resolution",
                i.e. set thy_per_tex higher here.
                This scan is *ranked* by Xi score,
                i.e. output texture file is ordered,
                and contains more detail (F2, F5, eps_avg).
    study -     Detailed study of a given texture.
                Data files are output to <report_dir>; 
                created if it does not exist, files
                are overwritten if present.
                Texture should be specified by a sequence
                of 7 numbers, separated by commas:
                XQ1,XQ2,Xu1,Xu2,Xd1,Xd2,Xd3

                Note: in our conventions XQ3=Xu3=0.

Usage:
    FNTexScan scan_all <max_Q> <texture_file> [--thy_per_tex=<tpt>] [--good_frac=<gf>] [--prior_uniform | --prior_wide ]
    FNTexScan rescan <texture_file> <report_file> [--top_N=<top_N>] [--thy_per_tex=<tpt>] [--opt_SM] [--prior_uniform | --prior_wide]
    FNTexScan study <texture> <report_dir> [--thy_per_tex=<tpt>] [--opt_SM] [--do_priors] [--inverted] [--prior_uniform | --prior_wide]
    FNTexScan -h | -help

Options:
    --thy_per_tex=<tpt>  Random theories per texture [default: 1000].
    --good_frac=<gf>     Fraction of theories within 5x to consider a texture "good" [default: 0.05].
    --top_N=<top_N>      Report only this number of top-ranked textures, -1 to report all [default: -1].
    --opt_SM             Optimize all parameters (not just epsilon) to fit the Standard Model.
    --prior_uniform      Use uniform prior instead of default lognormal.
    --prior_wide         Use a wider-than-normal lognormal prior.
    -h --help            Show this info.
"""

proc count_below*(scores: seq[float], target_score: float = 5.0): int =
    result = 0
    for s in scores:
        if s <= target_score:
            result += 1

# Returns: F2, F5, eps_avg
proc scanTexture(tex: FNTexture, thy_per_tex: int, prior_type: TheoryPriorType, optSM: bool):
    (float, float, float) =

    var all_costs, all_eps: seq[float] = @[]
    var optThy: FNTheory
    var opt_eps, opt_score: float

    for k in 0 .. thy_per_tex:
        optThy = makeRandomTheory(tex=tex, prior_type=prior_type)

        try:
            opt_eps = optimize_eps(optThy, brute_force = false, max_exp = true)
        except Exception as e:
            opt_eps = optimize_eps(optThy, brute_force = true, max_exp = true)

        optThy.C[^1] = ln(opt_eps)

        if optSM:
            optThy = optimize_SM(optThy)

        all_eps.add opt_eps
        opt_score = optThy.scoreSMDeviations(max_exp = true)
        all_costs.add opt_score

    let num_good_2 = count_below(all_costs, target_score = 2.0)
    let num_good_5 = count_below(all_costs, target_score = 5.0)

    let F2 = (num_good_2 / thy_per_tex).float
    let F5 = (num_good_5 / thy_per_tex).float
    let eps_avg = sum(all_eps.toTensor()) / thy_per_tex.float

    result = (F2, F5, eps_avg)


proc run_scan_all(max_Q: int, thy_per_tex: int, texture_file: string, good_frac: float, 
                    prior_uniform: bool, prior_wide: bool) =
    var all_textures = getAllValidTextures(max_Q)

    # Set prior type
    var prior_type: TheoryPriorType
    if prior_wide:
        prior_type = TheoryPriorType.wide_lognormal
    elif prior_uniform:
        prior_type = TheoryPriorType.uniform
    else:
        prior_type = TheoryPriorType.lognormal


    let f = open(texture_file, fmWrite)
    defer: f.close()

    var good_tex_indices: seq[int] = @[]
    var F2, F5, eps_avg: float

    for i, tex in suru(all_textures):
        (F2, F5, eps_avg) = scanTexture(
            tex=tex, 
            thy_per_tex=thy_per_tex, 
            prior_type=prior_type, 
            optSM=false)

        if F5 > good_frac: 
            good_tex_indices.add i
            f.writeLine(tex.toCompactString())

    echo $len(good_tex_indices), " good textures found."
    echo "Written to file: ", texture_file


proc run_rescan(texture_file: string, report_file: string, top_N: int, thy_per_tex: int, 
                opt_SM: bool, prior_uniform: bool, prior_wide: bool) =

    # Set prior type
    var prior_type: TheoryPriorType
    if prior_wide:
        prior_type = TheoryPriorType.wide_lognormal
    elif prior_uniform:
        prior_type = TheoryPriorType.uniform
    else:
        prior_type = TheoryPriorType.lognormal


    let strAllTex = readFile(texture_file)
    var allTex: seq[FNTexture] = @[]
    for texStr in strAllTex.split('\n'):
        allTex.add texStr.fromCompactString()

    var F2, F5, eps_avg: float

    var tableF2 = initOrderedTable[int, float]()
    var tableF5 = initOrderedTable[int, float]()
    var tableEpsAvg = initOrderedTable[int, float]()

    # Run the scan
    for i, tex in suru(allTex):
        (F2, F5, eps_avg) = scanTexture(
            tex=tex, 
            thy_per_tex=thy_per_tex, 
            prior_type=prior_type, 
            optSM=optSM)

        tableF2[i] = F2
        tableF5[i] = F5
        tableEpsAvg[i] = eps_avg

    # Sort tables
    proc sortValues(x,y: (int, float)): int = cmp(x[1], y[1])

    tableF2.sort(sortValues, order=SortOrder.Descending)
    
    var tex_i: int
    var this_tex: FNTexture
    var texStr: string

    let f_report = open(report_file, fmWrite)

    for key, val in tableF2:
        if tex_i < top_N:
            this_tex = allTex[key]

            texStr = this_tex.toCompactString()
            f_report.writeLine(texStr, " ", val, " ", tableF5[key], " ", tableEpsAvg[key])

        tex_i += 1


proc run_study(texture: string, report_dir: string, thy_per_tex: int, 
                opt_SM: bool, do_priors: bool, inverted: bool, 
                prior_uniform: bool, prior_wide: bool) =

    echo "run_study"
    let this_tex = fromCompactString(texture)
    echo this_tex

    # Set prior type
    var prior_type: TheoryPriorType
    if prior_wide:
        prior_type = TheoryPriorType.wide_lognormal
    elif prior_uniform:
        prior_type = TheoryPriorType.uniform
    else:
        prior_type = TheoryPriorType.lognormal

    # Make directory if it doesn't exist
    if not dirExists(report_dir):
        createDir(report_dir)

    os.setCurrentDir(report_dir)

    # Open report file paths, clobber if they exist
    # Setup data file names
    proc make_file_name(base_name: string): string =
        result = base_name
        if opt_SM:
            result = result & "_SMopt"
        if do_priors:
            result = result & "_prior"
        
        result = result & ".csv"

    var files = initTable[string, File]()

    let data_keys = ["thy_data", "SM_data", "xi_data", "eps_data", "tuning_data", "rotation_data"]

    for key in data_keys:
        # Clobber files
        let f = open(make_file_name(key), fmWrite)
        var header: string = ""
        if key == "thy_data":
            header = "Cu_re0, Cu_re1, Cu_re2, Cu_re3, Cu_re4, Cu_re5, Cu_re6, Cu_re7, Cu_re8, "
            header &= "Cu_im0, Cu_im1, Cu_im2, Cu_im3, Cu_im4, Cu_im5, Cu_im6, Cu_im7, Cu_im8, "
            header &= "Cd_re0, Cd_re1, Cd_re2, Cd_re3, Cd_re4, Cd_re5, Cd_re6, Cd_re7, Cd_re8, "
            header &= "Cd_im0, Cd_im1, Cd_im2, Cd_im3, Cd_im4, Cd_im5, Cd_im6, Cd_im7, Cd_im8, ln_eps"
        elif key == "SM_data":
            header = "Mu1, Mu2, Mu3, Md1, Md2, Md3, V_12, V_23, V_13, J"
        
        f.writeLine(header)
        f.close()

        # Open files for writing data
        files[key] = open(make_file_name(key), fmAppend)

    # Run texture scan and gather results
    var randThy, optThy, copyThy: FNTheory
    var opt_eps: float

    var all_eps, all_xi, all_nat: seq[float] = @[]
    var num_ud_flip = 0

    var Uu, Ud, Ku, Kd: Tensor[Complex64]

    for k in suru(0 ..< thy_per_tex):
        randThy = makeRandomTheory(tex=this_tex, prior_type=prior_type)

        if do_priors:
            optThy = randThy
        else:
            opt_eps = optimize_eps(randThy, brute_force = false, max_exp = true)
            optThy = FNTheory(tex: randThy.tex, C: randThy.C, hierarchy: randThy.hierarchy)
            optThy.C[^1] = ln(opt_eps)

        if opt_SM:
            optThy = optimize_SM(optThy)

        all_eps.add opt_eps
        all_xi.add score_SM_deviations(optThy, max_exp=true)
        all_nat.add naturalness_cost(optThy, sigma_weight=0.691)

        (Uu, Ud, Ku, Kd) = optThy.reportRotationMatrices()

        let optSMpars = optThy.reportSMParams()
        
        # Check for up/down mass inversion
        if all_xi[^1] <= 2 and optSMpars["Mu"][0] > optSMpars["Md"][0]:
            num_ud_flip += 1
    
        # Write data out
        files["thy_data"].writeLine(optThy.getRawParamString())
        files["SM_data"].writeLine(optThy.getSMParamString())
        files["xi_data"].writeLine($all_xi[^1])
        files["eps_data"].writeLine($opt_eps)
        files["tuning_data"].writeLine($all_nat[^1])

        files["rotation_data"].writeLine(optThy.getRotationString())        
#        files["rotation_data"].writeLine($Uu.toFlatSeq())
#        files["rotation_data"].writeLine($Ud.toFlatSeq())
#        files["rotation_data"].writeLine($Ku.toFlatSeq())
#        files["rotation_data"].writeLine($Kd.toFlatSeq())

    for filename, f in files:
        f.close()

    for target_xi in @[2.0, 5.0]:
        let num_good = count_below(all_xi, target_score=target_xi)
        echo "Xi < ", target_xi, " (normal): ",  num_good

    echo "Num_ud_flip: ", num_ud_flip


# Apparently this patch is needed for docopt to work here...?
proc toFloat(v: Value): float =
    parseFloat($v)


when isMainModule:
    # docopt hooks to be added here

    let args = docopt(doc, version = "FNTexScan v1.0")

    var ranCommand: bool = false

    ranCommand = ranCommand or args.dispatchProc(run_scan_all, "scan_all")
    ranCommand = ranCommand or args.dispatchProc(run_rescan, "rescan")
    ranCommand = ranCommand or args.dispatchProc(run_study, "study")

    if not ranCommand:
        echo doc






