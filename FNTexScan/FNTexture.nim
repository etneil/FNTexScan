import std/algorithm
import tables

type FNTexture* = object
    XQ*: array[3, int]
    Xu*: array[3, int]
    Xd*: array[3, int]
    isValid*: bool

proc isMonotonicInc(X: array[3, int]): bool =
    result = X[0] <= X[1] 
    result = result and X[1] <= X[2]

proc isMonotonicDecX2(X: array[3, int]): bool =
    result = X[0] >= X[1]

proc isMonotonicDec(X: array[3, int]): bool =
    result = X[0] >= X[1] 
    result = result and X[1] >= X[2]

proc isMonotonicAbsDec(X: array[3, int]): bool =
    if abs(X[2]) == 0:
        result = abs(X[0]) >= abs(X[1])
    else:
        result = abs(X[0]) >= abs(X[1])
        result = result and abs(X[1]) <= abs(X[2])

# Sorting algorithms for charges within a texture
proc sortAbs(x, y: int): int = cmp(abs(x), abs(y))

proc sortAbsPos(x, y: int): int = 
    if abs(x) == abs(y):
        cmp(x, y)
    else:
        cmp(abs(x), abs(y))

type
    ChargeSortAlgorithm = proc (x, y: int): int

proc sortCharges(tex: FNTexture, sortAlgo: ChargeSortAlgorithm = sortAbsPos): FNTexture =
    let XQ_sort = tex.XQ[0 .. 1].sorted(sortAlgo, Descending)
    result.XQ = [XQ_sort[0], XQ_sort[1], 0]

    let Xu_sort = tex.Xu[0 .. 1].sorted(sortAlgo, Descending)
    result.Xu = [Xu_sort[0], Xu_sort[1], 0]

    let Xd_sort = tex.Xd.sorted(sortAlgo, Descending)
    result.Xd = [Xd_sort[0], Xd_sort[1], Xd_sort[2]]

proc `==`(tex1, tex2: FNTexture): bool =
    let tex1_sort = tex1.sortCharges()
    let tex2_sort = tex2.sortCharges()

    result = tex1_sort.XQ == tex2_sort.XQ
    result = result and (tex1_sort.Xu == tex2_sort.Xu)
    result = result and (tex1_sort.Xd == tex2_sort.Xd)

proc mirror*(tex: FNTexture): FNTexture =
    result.XQ = [-tex.XQ[0], -tex.XQ[1], -tex.XQ[2]]
    result.Xu = [-tex.Xu[0], -tex.Xu[1], -tex.Xu[2]]
    result.Xd = [-tex.Xd[0], -tex.Xd[1], -tex.Xd[2]]

    result = result.sortCharges()


proc validateTexture(tex: FNTexture): bool =
    # Valid textures should have XQ[2] = 0 and Xu[2] = 0
    result = tex.XQ[2] == 0
    result = result and tex.Xu[2] == 0

    # Valid textures should be properly sorted already
    let texSort = tex.sortCharges()
    result = result and (tex.XQ == texSort.XQ)
    result = result and (tex.Xu == texSort.Xu)
    result = result and (tex.Xd == texSort.Xd)
    

proc initFNTexture*(XQ, Xu, Xd: array[3, int]): FNTexture =
    # Verify that the texture is valid
    result = FNTexture(XQ: XQ, Xu: Xu, Xd: Xd, isValid: false)
    result = result.sortCharges()
    result.isValid = result.validateTexture()

proc getAllValidTextures*(max_Q: int = 4): seq[FNTexture] =
    var allXQ: seq[array[3,int]] = @[]
    var allXu: seq[array[3,int]] = @[]
    var allXd: seq[array[3,int]] = @[]

    for a in -max_Q .. max_Q:
        for b in -max_Q .. max_Q:
           allXQ.add([a, b, 0])
           allXu.add([a, b, 0])
           for c in -max_Q .. max_Q:
               allXd.add [a, b, c]  

    var trialTex: FNTexture
    var resultAbsXQ = newTable[array[7,int], seq[FNTexture]]()

    for XQ in allXQ:
        for Xu in allXu:
            for Xd in allXd:
                let absXQ = [abs(XQ[0]), abs(XQ[1]), abs(Xu[0]), abs(Xu[1]), abs(Xd[0]), abs(Xd[1]), abs(Xd[2])]
                if not resultAbsXQ.hasKey(absXQ):
                    resultAbsXQ[absXQ] = @[]

                trialTex = initFNTexture(XQ=XQ, Xu=Xu, Xd=Xd)

                # Keep only if it's already canonically sorted
                if trialTex.XQ == XQ and trialTex.Xu == Xu and trialTex.Xd == Xd and trialTex.isValid:
                    # Check for mirror
                    if not (mirror(trialTex) in resultAbsXQ[absXQ]):
                        resultAbsXQ[absXQ].add trialTex

    for absXQ in resultAbsXQ.keys:
        for tex in resultAbsXQ[absXQ]:
            result.add tex
