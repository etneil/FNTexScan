import complex
import nimlapack
import arraymancer/private/sequninit


converter toLCF(x: Complex[float32]): lapack_complex_float =
    result.re = x.re
    result.im = x.im

converter toComplex(x: lapack_complex_float): Complex[float32] =
    result.re = x.re
    result.im = x.im

converter toLCF(x: Complex[float64]): lapack_complex_double =
    result.re = x.re
    result.im = x.im

converter toComplex(x: lapack_complex_double): Complex[float64] =
    result.re = x.re
    result.im = x.im

template first[A](x: openArray[A]): ptr A = addr(x[0])

proc eigh3*(a: array[9, Complex[float64]]): (array[3, float64], array[9, Complex[float64]]) =
    var aWork: array[9, lapack_complex_double]

    # Compute the eigensystem of the complex Hermitian 3x3 matrix a, using
    # LAPACK's "zheev" subroutine.

    # IMPORTANT NOTE: As of version 0.2.1 of nimlapack, this code requires
    # a hack to the package source files: exposing the lapack_complex_ types as 
    # public using *.

    # Now fixed in current version of nimlapack!

    for i, x in a:
        aWork[i] = x.toLCF()

    var
        jobz = cstring"V"
        uplo = cstring"L"
        N = 3.cint
        LDA = 3.cint
        info: cint
        w: array[3,cdouble]
        work: array[18, lapack_complex_double]
        lwork = -1.cint
        rwork: array[3*3-2, cdouble]

    # Workspace size query first    
    zheev(
        jobz,
        uplo,
        addr N,
        aWork.first,
        addr LDA,
        w.first,
        work.first,
        addr lwork,
        rwork.first,
        addr info
    )

    doAssert info == 0
    lwork = work[0].re.cint

    var work_opt = newSeqUninit[lapack_complex_double](lwork)

    zheev(
        jobz,
        uplo,
        addr N,
        aWork.first,
        addr LDA,
        w.first,
        work_opt.first,
        addr lwork,
        rwork.first,
        addr info
    )


    doAssert info == 0

    var U: array[9, Complex[float64]]

    for i, x in aWork:
        U[i] = x.toComplex()

    result = (w, U)


proc eigh3d*(a: array[9, Complex[float64]]): (array[3, float64], array[9, Complex[float64]]) =
    var aWork: array[9, lapack_complex_double]

    # Compute the eigensystem of the complex Hermitian 3x3 matrix a, using
    # LAPACK's "zheevd" subroutine.

    # IMPORTANT NOTE: As of version 0.2.1 of nimlapack, this code requires
    # a hack to the package source files: exposing the lapack_complex_ types as 
    # public using *.

    # Now fixed in current version of nimlapack!

    for i, x in a:
        aWork[i] = x.toLCF()

    var
        jobz = cstring"V"
        uplo = cstring"L"
        N = 3.cint
        LDA = 3.cint
        info: cint
        w: array[3,cdouble]
        work: array[18, lapack_complex_double]
        lwork = 18.cint
        rwork: array[36, cdouble]
        lrwork = 36.cint
        iwork: array[18,cint]
        liwork = 18.cint
    
    zheevd(
        jobz,
        uplo,
        addr N,
        aWork.first,
        addr LDA,
        w.first,
        work.first,
        addr lwork,
        rwork.first,
        addr lrwork,
        iwork.first,
        addr liwork,
        addr info
    )

    doAssert info == 0

    var U: array[9, Complex[float64]]

    for i, x in aWork:
        U[i] = x.toComplex()

    result = (w, U)


proc eigh3r*(a: array[9, Complex[float64]]): (array[3, float64], array[9, Complex[float64]]) =
    var aWork: array[9, lapack_complex_double]

    # Compute the eigensystem of the complex Hermitian 3x3 matrix a, using
    # LAPACK's "zheevr" subroutine.

    for i, x in a:
        aWork[i] = x.toLCF()

    var
        jobz = cstring"V"
        range = cstring"A"
        uplo = cstring"L"
        N = 3.cint
        LDA = 3.cint
        VL = 0.cdouble
        VU = 0.cdouble
        IL = 0.cint
        IU = 0.cint
        ABSTOL = 0.cdouble
        M: cint
        w: array[3,cdouble]
        Z: array[9, lapack_complex_double]
        LDZ = 3.cint
        ISUPPZ: array[6, cint]
        work: array[18, lapack_complex_double]
        lwork = 18.cint
        rwork: array[72, cdouble]
        lrwork = 72.cint
        iwork: array[30,cint]
        liwork = 30.cint    
        info: cint

    zheevr(
        jobz,
        range,
        uplo,
        addr N,
        aWork.first,
        addr LDA,
        addr VL,
        addr VU,
        addr IL,
        addr IU,
        addr ABSTOL,
        addr M,
        w.first,
        Z.first,
        addr LDZ,
        ISUPPZ.first,
        work.first,
        addr lwork,
        rwork.first,
        addr lrwork,
        iwork.first,
        addr liwork,
        addr info
    )

    doAssert info == 0

    var U: array[9, Complex[float64]]

    for i, x in Z:
        U[i] = x.toComplex()

    result = (w, U)

