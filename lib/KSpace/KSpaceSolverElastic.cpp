


#include <KSpace/KSpaceSolver.h>

using std::ios;
/* Shortcut for Simulation dimensions. */
using SD = Parameters::SimulationDimension;
/* Shortcut for Matrix id in the container. */
using MI = MatrixContainer::MatrixIdx;
/* Shortcut for Output stream id in the container. */
using OI = OutputStreamContainer::OutputStreamIdx;

KSpaceSolverElastic::ComputeElasticImp KSpaceSolverElastic::sComputeElasticImp
        {
                /*2D cases: SD rho0 bOnA c0 s0 alpha */
                {
                        std::make_tuple(SD::k2D, false, false, false, false, false),
                        &KSpaceSolverElastic::computeElastic<SD::k2D,
                                false, false, false, false, false>
                }
        };


void
KSpaceSolverElastic::allocateMemory()
{
    // Add matrices into the container and create all matrices
    mMatrixContainer.init();
    mMatrixContainer.createMatrices();

    // Add output streams into container
    mOutputStreamContainer.init(mMatrixContainer);
}

void
KSpaceSolverElastic::loadInputData()
{
    // Load data from the input file
    mMatrixContainer.loadDataFromKConfig();

    Hdf5File &outputFile = mParameters.getOutputFile();

    auto & koutput = mParameters.getKOutput();
    auto filename = koutput.getOutputFileName();

    if (koutput.isOutputToFileFlag() && !outputFile.canAccess(filename))
        outputFile.create(filename);

    mOutputStreamContainer.createStreams();
}

void
KSpaceSolverElastic::compute()
{
    /* Initialize all used FFTW plans */
    initializeFftwPlans();

    preProcessing<SD::k2D>();

#if 0
    std::cout << mParameters.getRho0ScalarFlag() << std::endl;
    std::cout << mParameters.getBOnAScalarFlag() << std::endl;
    std::cout << mParameters.getC0ScalarFlag() << std::endl;
    std::cout << mParameters.getS0ScalarFlag() << std::endl;
    std::cout << mParameters.getAlphaPower() << std::endl;
#endif

    sComputeElasticImp[std::make_tuple(
            mParameters.getSimulationDimension(),
            mParameters.getRho0ScalarFlag(),
            mParameters.getBOnAScalarFlag(),
            mParameters.getC0ScalarFlag(),
            mParameters.getS0ScalarFlag(),
            mParameters.getAlphaCoeffScalarFlag())]
            (*this);

    /* Post processing phase */
    // mPostProcessingTime.start();

    postProcessing();


}

template<Parameters::SimulationDimension simulationDimension,
        bool rho0ScalarFlag,
        bool bOnAScalarFlag,
        bool c0ScalarFlag,
        bool s0ScalarFlag,
        bool alphaCoefScalarFlag>
void
KSpaceSolverElastic::computeElastic()
{
    auto &params = mParameters;
    mActPercent = 0;
    // Set the actual progress percentage to correspond the time index
    // after recovery
    if (params.getTimeIndex() > 0)
        mActPercent = (100 * params.getTimeIndex()) / params.getNt();

    // Progress header
    Logger::log(Logger::LogLevel::kBasic, kOutFmtSimulationHeader);


    mIterationTime.start();
    params.setTimeIndex(0);

    // Execute main loop
    while (params.getTimeIndex() < params.getNt()
           && !params.isTimeToCheckpoint(mTotalTime))
    {
        const size_t timeIndex = params.getTimeIndex();

        computePressureGradient<simulationDimension>();

        computeSplitVelocity<simulationDimension>();

        computeVelocity<simulationDimension>();

        // Compute gradient of velocity
        computeVelocityGradient<simulationDimension>();

        computeSplitPressure<simulationDimension>();

        // Calculate initial pressure
        if ((timeIndex == 0) && (mParameters.getInitialPressureSourceFlag() == 1))
            addInitialPressureSource<simulationDimension, rho0ScalarFlag,
                    c0ScalarFlag, s0ScalarFlag>();

        // Compute new pressure
        computePressure<simulationDimension,
                rho0ScalarFlag,
                bOnAScalarFlag,
                c0ScalarFlag,
                s0ScalarFlag,
                alphaCoefScalarFlag>();

        storeSensorData();
        printStatistics();

        mParameters.incrementTimeIndex();
    } // Time loop
}

template<Parameters::SimulationDimension simulationDimension>
void
KSpaceSolverElastic::preProcessing()
{
    // Get the correct sensor mask and recompute indices
    auto &params = mParameters;

    if (params.getSensorMaskType() == Parameters::SensorMaskType::kIndex)
        getIndexMatrix(MI::kSensorMaskIndex).recomputeIndicesToCPP();

    if (!params.getRho0ScalarFlag())
    {
        getRealMatrix(MI::kDtRho0Sgx).scalarDividedBy(params.getDt());
        getRealMatrix(MI::kDtRho0Sgy).scalarDividedBy(params.getDt());
    }

    // Generate shift variables
    generateDerivativeOperators();

    // Generate absorption variables and kappa
    switch (mParameters.getAbsorbingFlag())
    {

        case Parameters::AbsorptionType::kLossless:
            generateKappa();
            break;

        default:
            assert(false);
    }

    // Generate PML
    generatePml();

    generateLameConstant();
}

void
KSpaceSolverElastic::generateDerivativeOperators()
{
    const DimensionSizes &dimensionSizes =
            mParameters.getFullDimensionSizes();
    const DimensionSizes &reducedDimensionSizes =
            mParameters.getReducedDimensionSizes();

    constexpr FloatComplex imagUnit = FloatComplex(0, 1);
    constexpr FloatComplex posExp = FloatComplex(0, 1);
    constexpr FloatComplex negExp = FloatComplex(0, -1);
    constexpr double pi2 = 2.0f * double(M_PI);

    const double dx = mParameters.getDx();
    const double dy = mParameters.getDy();

    FloatComplex *ddxKShiftPos = getComplexData(MI::kDdxKShiftPosR);
    FloatComplex *ddyKShiftPos = getComplexData(MI::kDdyKShiftPos);

    FloatComplex *ddxKShiftNeg = getComplexData(MI::kDdxKShiftNegR);
    FloatComplex *ddyKShiftNeg = getComplexData(MI::kDdyKShiftNeg);

    // Calculate ifft shift
    auto iFftShift = [](ptrdiff_t i, ptrdiff_t size)
    {
        return (i + (size / 2)) % size - (size / 2);
    };// end of iFftShift

    for (size_t i = 0; i < reducedDimensionSizes.nx; i++)
    {
        const ptrdiff_t shift = iFftShift(i, dimensionSizes.nx);
        const double kx = (pi2 / dx) * (double(shift)
                                       / double(dimensionSizes.nx));
        const double exponent = kx * dx * 0.5f;

        ddxKShiftPos[i] = imagUnit * kx * std::exp(posExp * exponent);
        ddxKShiftNeg[i] = imagUnit * kx * std::exp(negExp * exponent);
    }

    // ddyKShiftPos, ddyKShiftPos
    for (size_t i = 0; i < dimensionSizes.ny; i++)
    {
        const ptrdiff_t shift = iFftShift(i, dimensionSizes.ny);
        const double ky = (pi2 / dy) * (double(shift)
                                       / double(dimensionSizes.ny));
        const double exponent = ky * dy * 0.5f;

        ddyKShiftPos[i] = imagUnit * ky * std::exp(posExp * exponent);
        ddyKShiftNeg[i] = imagUnit * ky * std::exp(negExp * exponent);
    }
}

template<Parameters::SimulationDimension simulationDimension,
        bool rho0ScalarFlag,
        bool c0ScalarFlag,
        bool s0ScalarFlag>
void
KSpaceSolverElastic::addInitialPressureSource()
{
    const size_t nElements = mParameters.getFullDimensionSizes().nElements();

    double *p0 = getRealData(MI::kInitialPressureSourceInput);

    double *tmp = getRealData(MI::kTmpReal1);

    for (size_t i = 0; i < nElements; i++)
    {
        tmp[i] = -p0[i] / 2.f;
    }

    getRealMatrix(MI::kSxxSplitX).copyData(getRealMatrix(MI::kTmpReal1));
    getRealMatrix(MI::kSxxSplitY).copyData(getRealMatrix(MI::kTmpReal1));
    getRealMatrix(MI::kSyySplitX).copyData(getRealMatrix(MI::kTmpReal1));
    getRealMatrix(MI::kSyySplitY).copyData(getRealMatrix(MI::kTmpReal1));

#if 0
    double* sxxSplitX = getRealData(MI::kSxxSplitX);
    double* sxxSplitY = getRealData(MI::kSxxSplitY);

    double* syySplitX = getRealData(MI::kSyySplitX);
    double* syySplitY = getRealData(MI::kSyySplitY);

    double* pmat = getRealData(MI::kSxxSplitY);
#endif
}

template<Parameters::SimulationDimension simulationDimension,
        bool rho0ScalarFlag,
        bool bOnAScalarFlag,
        bool c0ScalarFlag,
        bool s0ScalarFlag,
        bool alphaCoefScalarFlag>
void
KSpaceSolverElastic::computePressure()
{
    const auto &fullDimSizes = mParameters.getFullDimensionSizes();

    double *p = getRealData(MI::kP);

    const double *sxxSplitX = getRealData(MI::kSxxSplitX);
    const double *sxxSplitY = getRealData(MI::kSxxSplitY);

    const double *syySplitX = getRealData(MI::kSyySplitX);
    const double *syySplitY = getRealData(MI::kSyySplitY);

    /*
     p = -(sxx_split_x + sxx_split_y + syy_split_x + syy_split_y) / 2;
     */
    for (size_t i = 0; i < fullDimSizes.nElements(); ++i)
        p[i] = -(sxxSplitX[i] + sxxSplitY[i] + syySplitX[i] + syySplitY[i]) / 2;
}


template<Parameters::SimulationDimension simulationDimension>
void
KSpaceSolverElastic::computePressureGradient()
{
    const auto &params = mParameters;
    const auto &fullDimSizes = params.getFullDimensionSizes();

    const FloatComplex *ddxKShiftPos = getComplexData(MI::kDdxKShiftPosR);
    const FloatComplex *ddyKShiftPos = getComplexData(MI::kDdyKShiftPos);

    const FloatComplex *ddxKShiftNeg = getComplexData(MI::kDdxKShiftNegR);
    const FloatComplex *ddyKShiftNeg = getComplexData(MI::kDdyKShiftNeg);

    double *tmpReal1 = getRealData(MI::kTmpReal1);
    double *tmpReal2 = getRealData(MI::kTmpReal2);
    double *tmpReal3 = getRealData(MI::kTmpReal3);
    double *tmpReal4 = getRealData(MI::kTmpReal4);

    FloatComplex *ifftXXdx = getComplexData(MI::kTmpFftwXXdx);
    FloatComplex *ifftYYdy = getComplexData(MI::kTmpFftwYYdy);
    FloatComplex *ifftXYdx = getComplexData(MI::kTmpFftwXYdx);
    FloatComplex *ifftXYdy = getComplexData(MI::kTmpFftwXYdy);

    const double *sxxSplitX = getRealData(MI::kSxxSplitX);
    const double *sxxSplitY = getRealData(MI::kSxxSplitY);

    const double *syySplitX = getRealData(MI::kSyySplitX);
    const double *syySplitY = getRealData(MI::kSyySplitY);

    const double *sxySplitX = getRealData(MI::kSxySplitX);
    const double *sxySPlitY = getRealData(MI::kSxySplitY);

    auto nx = params.getFullDimensionSizes().nx;
    auto ny = params.getFullDimensionSizes().ny;

    const double dividerX = 1.0f / static_cast<double>(nx);
    const double dividerY = 1.0f / static_cast<double>(ny);

    /*
    dsxxdx = real( ifft( bsxfun(@times, ddx_k_shift_pos,
                   fft(sxx_split_x + sxx_split_y, [], 1)), [], 1) );
    dsyydy = real( ifft( bsxfun(@times, ddy_k_shift_pos,
                    fft(syy_split_x + syy_split_y, [], 2)), [], 2) );
    dsxydx = real( ifft( bsxfun(@times, ddx_k_shift_neg,
                    fft(sxy_split_x + sxy_split_y, [], 1)), [], 1) );
    dsxydy = real( ifft( bsxfun(@times, ddy_k_shift_neg,
                    fft(sxy_split_x + sxy_split_y, [], 2)), [], 2) );
    */

    /* tmpRet = sxx_split_x + sxx_split_y */
    for (size_t i = 0; i < fullDimSizes.nElements(); i++)
    {
        tmpReal1[i] = sxxSplitX[i] + sxxSplitY[i];
        tmpReal2[i] = syySplitX[i] + syySplitY[i];
        tmpReal3[i] = sxySplitX[i] + sxySPlitY[i];
        tmpReal4[i] = sxySplitX[i] + sxySPlitY[i];
    }

    /* tmpRet = fft(tmpRet, [], 1) */
    getTempFftwXXdx().computeR2CFft1DX(getRealMatrix(MI::kTmpReal1));
    getTempFftwYYdy().computeR2CFft1DY(getRealMatrix(MI::kTmpReal2));
    getTempFftwXYdx().computeR2CFft1DX(getRealMatrix(MI::kTmpReal3));
    getTempFftwXYdy().computeR2CFft1DY(getRealMatrix(MI::kTmpReal4));

    /* tmpRet = bsxfun(@times, ddx_k_shift_pos, tmpRet) */

    auto reducedDimXSizes = params.getReducedXDimensionSizes();

    for (size_t z = 0; z < reducedDimXSizes.nz; z++)
    {
        for (size_t y = 0; y < reducedDimXSizes.ny; y++)
        {
            for (size_t x = 0; x < reducedDimXSizes.nx; x++)
            {
                const size_t i = get1DIndex(z, y, x, reducedDimXSizes);

                ifftXXdx[i] = ifftXXdx[i] * ddxKShiftPos[x] * dividerX;
                ifftXYdx[i] = ifftXYdx[i] * ddxKShiftNeg[x] * dividerX;
            }
        }
    }

    auto reducedDimYSizes = params.getReducedYDimensionSizes();

    for (size_t z = 0; z < reducedDimYSizes.nz; z++)
    {
        for (size_t y = 0; y < reducedDimYSizes.ny; y++)
        {
            for (size_t x = 0; x < reducedDimYSizes.nx; x++)
            {
                const size_t i = get1DIndex(z, y, x, reducedDimYSizes);

                ifftYYdy[i] = ifftYYdy[i] * ddyKShiftPos[y] * dividerY;
                ifftXYdy[i] = ifftXYdy[i] * ddyKShiftNeg[y] * dividerY;
            }
        }
    }

    /* real(ifft(tmpRet)) */
    getTempFftwXXdx().computeC2RFft1DX(getRealMatrix(MI::kDSxxdx));
    getTempFftwYYdy().computeC2RFft1DY(getRealMatrix(MI::kDSyydy));
    getTempFftwXYdx().computeC2RFft1DX(getRealMatrix(MI::kDSxydx));
    getTempFftwXYdy().computeC2RFft1DY(getRealMatrix(MI::kDSxydy));

#if 0
    auto pmat1 = getRealData(MI::kDSxxdx);
    auto pmat2 = getRealData(MI::kDSyydy);
    auto pmat3 = getRealData(MI::kDSxydx);
    auto pmat4 = getRealData(MI::kDSxydy);
#endif

}

template<Parameters::SimulationDimension simulationDimension>
void
KSpaceSolverElastic::computeSplitVelocity()
{
    const auto &fullDimSizes = mParameters.getFullDimensionSizes();

    double *dsXXdx = getRealData(MI::kDSxxdx);
    double *dsYYdy = getRealData(MI::kDSyydy);

    double *dsXYdx = getRealData(MI::kDSxydx);
    double *dsXYdy = getRealData(MI::kDSxydy);

    double *mPmlX = getRealData(MI::kMPmlX);
    double *mPmlY = getRealData(MI::kMPmlY);

    double *mPmlXSgx = getRealData(MI::kMPmlXSgx);
    double *mPmlYSgy = getRealData(MI::kMPmlYSgy);

    double *pmlX = getRealData(MI::kPmlX);
    double *pmlY = getRealData(MI::kPmlY);

    double *pmlXSgx = getRealData(MI::kPmlXSgx);
    double *pmlYSgy = getRealData(MI::kPmlYSgy);

    double *uxSplitX = getRealData(MI::kUxSplitX);
    double *uxSplitY = getRealData(MI::kUxSplitY);
    double *uySplitX = getRealData(MI::kUySplitX);
    double *uySplitY = getRealData(MI::kUySplitY);

    double *dtRho0Sgx = getRealData(MI::kDtRho0Sgx);
    double *dtRho0Sgy = getRealData(MI::kDtRho0Sgy);

    /*
     ux_split_x = bsxfun(@times, mpml_y, ...
                    bsxfun(@times, pml_x_sgx, ...
                        bsxfun(@times, mpml_y, ...
                            bsxfun(@times, pml_x_sgx, ux_split_x) ...
                        ) + kgrid.dt .* rho0_sgx_inv .* dsxxdx ...
                     ) ...
                  );

    ux_split_y = bsxfun(@times, mpml_x_sgx, bsxfun(@times, pml_y, ...
                 bsxfun(@times, mpml_x_sgx, bsxfun(@times, pml_y, ux_split_y)) ...
                 + kgrid.dt .* rho0_sgx_inv .* dsxydy));

    uy_split_x = bsxfun(@times, mpml_y_sgy, bsxfun(@times, pml_x, ...
                 bsxfun(@times, mpml_y_sgy, bsxfun(@times, pml_x, uy_split_x)) ...
                 + kgrid.dt .* rho0_sgy_inv .* dsxydx));

    uy_split_y = bsxfun(@times, mpml_x,     bsxfun(@times, pml_y_sgy, ...
                 bsxfun(@times, mpml_x, bsxfun(@times, pml_y_sgy, uy_split_y)) ...
                 + kgrid.dt .* rho0_sgy_inv .* dsyydy));
    */

    for (size_t z = 0; z < fullDimSizes.nz; z++)
        for (size_t y = 0; y < fullDimSizes.ny; y++)
            for (size_t x = 0; x < fullDimSizes.nx; x++)
            {
                const size_t i = get1DIndex(z, y, x, fullDimSizes);

                auto coef = mPmlY[y] * pmlXSgx[x];
                auto expr = coef * uxSplitX[i] + dtRho0Sgx[i] * dsXXdx[i];
                uxSplitX[i] = coef * expr;

                coef = mPmlXSgx[x] * pmlY[y];
                expr = coef * uxSplitY[i] + dtRho0Sgx[i] * dsXYdy[i];
                uxSplitY[i] = coef * expr;

                coef = mPmlYSgy[y] * pmlX[x];
                expr = coef * uySplitX[i] + dtRho0Sgy[i] * dsXYdx[i];
                uySplitX[i] = coef * expr;

                coef = mPmlX[x] * pmlYSgy[y];
                expr = coef * uySplitY[i] + dtRho0Sgy[i] * dsYYdy[i];
                uySplitY[i] = coef * expr;
            }
}

template<Parameters::SimulationDimension simulationDimension>
void
KSpaceSolverElastic::computeVelocity()
{
    const auto &fullDimSizes = mParameters.getFullDimensionSizes();

    double *uxSplitX = getRealData(MI::kUxSplitX);
    double *uxSplitY = getRealData(MI::kUxSplitY);
    double *uySplitX = getRealData(MI::kUySplitX);
    double *uySplitY = getRealData(MI::kUySplitY);

    double *uxSgx = getRealData(MI::kUxSgx);
    double *uySgy = getRealData(MI::kUySgy);

    /* ux_sgx = ux_split_x + ux_split_y;
       uy_sgy = uy_split_x + uy_split_y;
    */

    for (size_t i = 0; i < fullDimSizes.nElements(); ++i)
    {
        uxSgx[i] = uxSplitX[i] + uxSplitY[i];
        uySgy[i] = uySplitX[i] + uySplitY[i];
    }
}

template<Parameters::SimulationDimension simulationDimension>
void
KSpaceSolverElastic::computeVelocityGradient()
{
    const auto &params = mParameters;

    auto nx = params.getFullDimensionSizes().nx;
    auto ny = params.getFullDimensionSizes().ny;

    const double dividerX = 1.0f / static_cast<double>(nx);
    const double dividerY = 1.0f / static_cast<double>(ny);


    const FloatComplex *ddxKShiftPos = getComplexData(MI::kDdxKShiftPosR);
    const FloatComplex *ddyKShiftPos = getComplexData(MI::kDdyKShiftPos);

    const FloatComplex *ddxKShiftNeg = getComplexData(MI::kDdxKShiftNegR);
    const FloatComplex *ddyKShiftNeg = getComplexData(MI::kDdyKShiftNeg);


    FloatComplex *ifftXXdx = getComplexData(MI::kTmpFftwXXdx);
    FloatComplex *ifftXXdy = getComplexData(MI::kTmpFftwXXdy);
    FloatComplex *ifftYYdx = getComplexData(MI::kTmpFftwYYdx);
    FloatComplex *ifftYYdy = getComplexData(MI::kTmpFftwYYdy);


    /*
    duxdx = real( ifft( bsxfun(@times, ddx_k_shift_neg,
                fft(ux_sgx, [], 1)), [], 1));
    duxdy = real( ifft( bsxfun(@times, ddy_k_shift_pos,
                fft(ux_sgx, [], 2)), [], 2));
    duydx = real( ifft( bsxfun(@times, ddx_k_shift_pos,
                fft(uy_sgy, [], 1)), [], 1));
    duydy = real( ifft( bsxfun(@times, ddy_k_shift_neg,
                fft(uy_sgy, [], 2)), [], 2));
    */

    getTempFftwXXdx().computeR2CFft1DX(getRealMatrix(MI::kUxSgx));
    getTempFftwXXdy().computeR2CFft1DY(getRealMatrix(MI::kUxSgx));

    getTempFftwYYdx().computeR2CFft1DX(getRealMatrix(MI::kUySgy));
    getTempFftwYYdy().computeR2CFft1DY(getRealMatrix(MI::kUySgy));

    auto reducedDimXSizes = params.getReducedXDimensionSizes();

    for (size_t z = 0; z < reducedDimXSizes.nz; z++)
        for (size_t y = 0; y < reducedDimXSizes.ny; y++)
            for (size_t x = 0; x < reducedDimXSizes.nx; x++)
            {
                const size_t i = get1DIndex(z, y, x, reducedDimXSizes);

                ifftXXdx[i] *= ddxKShiftNeg[x] * dividerX;
                ifftYYdx[i] *= ddxKShiftPos[x] * dividerX;
            }

    auto reducedDimYSizes = params.getReducedYDimensionSizes();

    for (size_t z = 0; z < reducedDimYSizes.nz; z++)
        for (size_t y = 0; y < reducedDimYSizes.ny; y++)
            for (size_t x = 0; x < reducedDimYSizes.nx; x++)
            {
                const size_t i = get1DIndex(z, y, x, reducedDimYSizes);

                ifftXXdy[i] *= ddyKShiftPos[y] * dividerY;
                ifftYYdy[i] *= ddyKShiftNeg[y] * dividerY;
            }

    getTempFftwXXdx().computeC2RFft1DX(getRealMatrix(MI::kDuxdx));
    getTempFftwXXdy().computeC2RFft1DY(getRealMatrix(MI::kDuxdy));

    getTempFftwYYdx().computeC2RFft1DX(getRealMatrix(MI::kDuydx));
    getTempFftwYYdy().computeC2RFft1DY(getRealMatrix(MI::kDuydy));

#if 0
    double *duxdx = getRealData(MI::kDuxdx);
    double *duxdy = getRealData(MI::kDuxdy);
    double *duydx = getRealData(MI::kDuydx);
    double *duydy = getRealData(MI::kDuydy);

    auto pmat1 = getRealData(MI::kDuxdx);
    auto pmat2 = getRealData(MI::kDuydx);
    auto pmat3 = getRealData(MI::kDuxdy);
    auto pmat4 = getRealData(MI::kDuydy);
#endif
}


template<Parameters::SimulationDimension simulationDimension>
void
KSpaceSolverElastic::computeSplitPressure()
{
    const auto &fullDimSizes = mParameters.getFullDimensionSizes();

    const double dt = mParameters.getDt();

    double *mu = getRealData(MI::kMu);
    double *lambda = getRealData(MI::kLambda);
    double *musgxy = getRealData(MI::kMuSgxy);

    double *sxxSplitX = getRealData(MI::kSxxSplitX);
    double *sxxSplitY = getRealData(MI::kSxxSplitY);

    double *syySplitX = getRealData(MI::kSyySplitX);
    double *syySplitY = getRealData(MI::kSyySplitY);

    double *sxySplitX = getRealData(MI::kSxySplitX);
    double *sxySplitY = getRealData(MI::kSxySplitY);

    double *mPmlX = getRealData(MI::kMPmlX);
    double *mPmlY = getRealData(MI::kMPmlY);

    double *mPmlXSgx = getRealData(MI::kMPmlXSgx);
    double *mPmlYSgy = getRealData(MI::kMPmlYSgy);

    double *pmlX = getRealData(MI::kPmlX);
    double *pmlY = getRealData(MI::kPmlY);

    double *pmlXSgx = getRealData(MI::kPmlXSgx);
    double *pmlYSgy = getRealData(MI::kPmlYSgy);

    double *duxdx = getRealData(MI::kDuxdx);
    double *duxdy = getRealData(MI::kDuxdy);
    double *duydx = getRealData(MI::kDuydx);
    double *duydy = getRealData(MI::kDuydy);

    /*
     sxx_split_x = bsxfun(@times, mpml_y, ...
                        bsxfun(@times, pml_x, ...
                            bsxfun(@times, mpml_y, bsxfun(@times, pml_x, sxx_split_x)) ...
                  + kgrid.dt .* (2 .* mu + lambda) .* duxdx ));

    sxx_split_y = bsxfun(@times, mpml_x, bsxfun(@times, pml_y, ...
                  bsxfun(@times, mpml_x, bsxfun(@times, pml_y, sxx_split_y)) ...
                  + kgrid.dt .* lambda .* duydy ));

    syy_split_x = bsxfun(@times, mpml_y, bsxfun(@times, pml_x, ...
                  bsxfun(@times, mpml_y, bsxfun(@times, pml_x, syy_split_x)) ...
                  + kgrid.dt .* lambda .* duxdx));
    syy_split_y = bsxfun(@times, mpml_x, bsxfun(@times, pml_y, ...
                  bsxfun(@times, mpml_x, bsxfun(@times, pml_y, syy_split_y)) ...
                  + kgrid.dt .* (2 .* mu + lambda) .* duydy ));

    sxy_split_x = bsxfun(@times, mpml_y_sgy, bsxfun(@times, pml_x_sgx, ...
                  bsxfun(@times, mpml_y_sgy, bsxfun(@times, pml_x_sgx, sxy_split_x)) ...
                  + kgrid.dt .* mu_sgxy .* duydx));
    sxy_split_y = bsxfun(@times, mpml_x_sgx, bsxfun(@times, pml_y_sgy, ...
                  bsxfun(@times, mpml_x_sgx, bsxfun(@times, pml_y_sgy, sxy_split_y)) ...
                  + kgrid.dt .* mu_sgxy .* duxdy));
    */


    for (size_t z = 0; z < fullDimSizes.nz; z++)
        for (size_t y = 0; y < fullDimSizes.ny; y++)
            for (size_t x = 0; x < fullDimSizes.nx; x++)
            {
                const size_t i = get1DIndex(z, y, x, fullDimSizes);
                auto moduli = (2.f * mu[i] + lambda[i]);

                auto coef = mPmlY[y] * pmlX[x];
                auto expr = coef * sxxSplitX[i] + dt * moduli * duxdx[i];
                sxxSplitX[i] = coef * expr;

                coef = mPmlX[x] * pmlY[y];
                expr = coef * sxxSplitY[i] + dt * lambda[i] * duydy[i];
                sxxSplitY[i] = coef * expr;

                coef = mPmlY[y] * pmlX[x];
                expr = coef * syySplitX[i] + dt * lambda[i] * duxdx[i];
                syySplitX[i] = coef * expr;

                coef = mPmlX[x] * pmlY[y];
                expr = coef * syySplitY[i] + dt * moduli * duydy[i];
                syySplitY[i] = coef * expr;

                coef = mPmlYSgy[y] * pmlXSgx[x];
                expr = coef * sxySplitX[i] + dt * musgxy[i] * duydx[i];
                sxySplitX[i] = coef * expr;

                coef = mPmlXSgx[x] * pmlYSgy[y];
                expr = coef * sxySplitY[i] + dt * musgxy[i] * duxdy[i];
                sxySplitY[i] = coef * expr;
            }
}

void
KSpaceSolverElastic::generateLameConstant()
{
    auto &params = mParameters;

    if (!params.getElasticFlag() || params.getC0ScalarFlag()
        || params.getS0ScalarFlag() || params.getRho0ScalarFlag())
        return;

    const size_t nElements = mParameters.getFullDimensionSizes().nElements();

    auto c2 = getRealData(MI::kC2);
    auto s2 = getRealData(MI::kS2);
    auto rho0 = getRealData(MI::kRho0);
    auto mu = getRealData(MI::kMu);
    auto lambda = getRealData(MI::kLambda);
    auto musgxy = getRealData(MI::kMuSgxy);

    for (size_t i = 0; i < nElements; i++)
    {
        mu[i] = s2[i] * s2[i] * rho0[i];
        lambda[i] = c2[i] * c2[i] * rho0[i] - 2 * mu[i];
        musgxy[i] = mu[i];
    }
}


void
KSpaceSolverElastic::generatePml()
{
    const DimensionSizes dimensionSizes = mParameters.getFullDimensionSizes();

    double pmlXAlpha = mParameters.getPmlXAlpha();
    double pmlYAlpha = mParameters.getPmlYAlpha();

    const size_t pmlXSize = mParameters.getPmlXSize();
    const size_t pmlYSize = mParameters.getPmlYSize();

    const double cRefDx = mParameters.getCRef() / mParameters.getDx();
    const double cRefDy = mParameters.getCRef() / mParameters.getDy();

    const double dt2 = mParameters.getDt() * 0.5f;

    double *pmlX = getRealData(MI::kPmlX);
    double *pmlY = getRealData(MI::kPmlY);

    double *pmlXSgx = getRealData(MI::kPmlXSgx);
    double *pmlYSgy = getRealData(MI::kPmlYSgy);

    auto mPmlX = getRealData(MI::kMPmlX);
    auto mPmlY = getRealData(MI::kMPmlY);

    auto mPmlXSgx = getRealData(MI::kMPmlXSgx);
    auto mPmlYSgy = getRealData(MI::kMPmlYSgy);

    auto multiAxialPmlRatio = mParameters.getMultiAxialPmlRatio();

    // Init arrays
    auto initPml = [](double *pml, double *pmlSg, size_t size)
    {
        for (size_t i = 0; i < size; i++)
        {
            pml[i] = 1.0f;
            pmlSg[i] = 1.0f;
        }
    };// end of initPml

    // Calculate left value of PML exponent,
    // for staggered use i + 0.5f, i shifted by -1 (Matlab indexing).
    auto pmlLeft = [dt2](double i, double cRef, double pmlAlpha, double pmlSize)
    {
        return exp(-dt2 * pmlAlpha * cRef * pow((i - pmlSize) / (-pmlSize), 4));
    };// end of pmlLeft.

    // Calculate right value of PML exponent,
    // for staggered use i + 0.5f, i shifted by +1 (Matlab indexing).
    auto pmlRight = [dt2](double i, double cRef, double pmlAlpha, double pmlSize)
    {
        return exp(-dt2 * pmlAlpha * cRef * pow((i + 1.0f) / pmlSize, 4));
    };// end of pmlRight.

    // PML in x dimension
    initPml(pmlX, pmlXSgx, dimensionSizes.nx);

    // Too difficult for SIMD
    for (size_t i = 0; i < pmlXSize; i++)
    {
        pmlX[i] = pmlLeft(double(i), cRefDx, pmlXAlpha, pmlXSize);
        pmlXSgx[i] = pmlLeft(double(i) + 0.5f, cRefDx, pmlXAlpha, pmlXSize);

        const size_t iR = dimensionSizes.nx - pmlXSize + i;

        pmlX[iR] = pmlRight(double(i), cRefDx, pmlXAlpha, pmlXSize);
        pmlXSgx[iR] = pmlRight(double(i) + 0.5f, cRefDx, pmlXAlpha, pmlXSize);
    }

    // PML in y dimension
    initPml(pmlY, pmlYSgy, dimensionSizes.ny);

    // Too difficult for SIMD
    for (size_t i = 0; i < pmlYSize; i++)
    {
        if (!mParameters.isSimulationAS())
        { // for axisymmetric code the PML is only on the outer side
            pmlY[i] = pmlLeft(double(i), cRefDy, pmlYAlpha, pmlYSize);
            pmlYSgy[i] = pmlLeft(double(i) + 0.5f, cRefDy, pmlYAlpha, pmlYSize);
        }

        const size_t iR = dimensionSizes.ny - pmlYSize + i;

        pmlY[iR] = pmlRight(double(i), cRefDy, pmlYAlpha, pmlYSize);
        pmlYSgy[iR] = pmlRight(double(i) + 0.5f, cRefDy, pmlYAlpha, pmlYSize);
    }

    pmlXAlpha *= multiAxialPmlRatio;
    pmlYAlpha *= multiAxialPmlRatio;

    initPml(mPmlX, mPmlXSgx, dimensionSizes.nx);

    // Too difficult for SIMD
    for (size_t i = 0; i < pmlXSize; i++)
    {
        mPmlX[i] = pmlLeft(double(i), cRefDx, pmlXAlpha, pmlXSize);
        mPmlXSgx[i] = pmlLeft(double(i) + 0.5f, cRefDx, pmlXAlpha, pmlXSize);

        const size_t iR = dimensionSizes.nx - pmlXSize + i;

        mPmlX[iR] = pmlRight(double(i), cRefDx, pmlXAlpha, pmlXSize);
        mPmlXSgx[iR] = pmlRight(double(i) + 0.5f, cRefDx, pmlXAlpha, pmlXSize);
    }

    // PML in y dimension
    initPml(mPmlY, mPmlYSgy, dimensionSizes.ny);

    // Too difficult for SIMD
    for (size_t i = 0; i < pmlYSize; i++)
    {
        if (!mParameters.isSimulationAS())
        { // for axisymmetric code the PML is only on the outer side
            mPmlY[i] = pmlLeft(double(i), cRefDy, pmlYAlpha, pmlYSize);
            mPmlYSgy[i] = pmlLeft(double(i) + 0.5f, cRefDy, pmlYAlpha, pmlYSize);
        }

        const size_t iR = dimensionSizes.ny - pmlYSize + i;

        mPmlY[iR] = pmlRight(double(i), cRefDy, pmlYAlpha, pmlYSize);
        mPmlYSgy[iR] = pmlRight(double(i) + 0.5f, cRefDy, pmlYAlpha, pmlYSize);
    }
}


void KSpaceSolverElastic::storeSensorInfo()
{

    // Unless the time for sampling has come, exit.
    if (mParameters.getTimeIndex() >= mParameters.getSamplingStartTimeIndex())
    {
        if (mParameters.getStoreVelocityNonStaggeredRawFlag())
        {
            if (mParameters.isSimulation3D())
            {
                computeShiftedVelocity<SD::k3D>();
            }
            else
            {
                computeShiftedVelocity<SD::k2D>();
            }
        }
        mOutputStreamContainer.sampleStreams();
    }
}

void KSpaceSolverElastic::postProcessing()
{

    auto & params = mParameters;
    auto koutput = params.getKOutput();

    if (koutput.isOutputToFileFlag()) {
        if (mParameters.getStorePressureFinalAllFlag())
        {
            getRealMatrix(MI::kP).writeData(mParameters.getOutputFile(),
                                            mOutputStreamContainer.getStreamHdf5Name(OI::kFinalPressure),
                                            mParameters.getCompressionLevel());
        }// p_final
        // Apply post-processing and close
        mOutputStreamContainer.postProcessStreams();
        mOutputStreamContainer.closeStreams();

        writeOutputDataInfo();
        params.getOutputFile().close();
    }

}


void
KSpaceSolverElastic::fftwVerify()
{
//    initializeFftwPlans();
//    generateDerivativeOperators();

#if 0
    FloatComplex *fftx = getComplexData(MI::kTmpFftwXXdx);
    double *x = getRealData(MI::kTmpReal1);

    x[49] = -5;

    x[49 + 128 - 3] = -5;
    x[49 + 128 - 2] = -5;
    x[49 + 128 - 1] = -5;
    x[49 + 128] = -5;
    x[49 + 128 + 1] = -5;
    x[49 + 128 + 2] = -5;
    x[49 + 128 + 3] = -5;

    getTempFftwXXdx().computeR2CFft1DX(getRealMatrix(MI::kTmpReal1));

    std::cout << fftx[0] << std::endl;
    std::cout << fftx[1] << std::endl;

    std::cout << " -------------------- " << std::endl;
    std::cout << fftx[62] << std::endl;
    std::cout << fftx[63] << std::endl;
    std::cout << fftx[64] << std::endl;
    std::cout << fftx[65] << std::endl;
    std::cout << fftx[66] << std::endl;

    std::cout << " -------------------- " << std::endl;

    std::cout << fftx[64] << std::endl;
    std::cout << fftx[65] << std::endl;

    std::cout << " -------------------- " << std::endl;
    std::cout << fftx[64+62] << std::endl;
    std::cout << fftx[64+63] << std::endl;
    std::cout << fftx[64+64] << std::endl;
    std::cout << fftx[64+65] << std::endl;
    std::cout << fftx[64+66] << std::endl;
#endif

#if 0
    FloatComplex *fftx = getComplexData(MI::kTmpFftwXXdx);
    const double divider = 1.0f / double(mParameters.getFullDimensionSizes().nElements());
    const auto & dimSizes = mParameters.getFullDimensionSizes();
    const size_t nElements = mParameters.getFullDimensionSizes().nElements();

    const DimensionSizes& reducedDimSizes= mParameters.getReducedDimensionSizes();

    double *p0 = getRealData(MI::kInitialPressureSourceInput);

    double* sxxSplitX = getRealData(MI::kSxxSplitX);
    double* sxxSplitY = getRealData(MI::kSxxSplitY);

    double* tmp = getRealData(MI::kTmpReal1);

    for (size_t i = 0; i < nElements; i++) {
        tmp[i] = -p0[i] / 2.f;
    }

    getRealMatrix(MI::kSxxSplitX).copyData(getRealMatrix(MI::kTmpReal1));


    getRealMatrix(MI::kSxxSplitY).copyData(getRealMatrix(MI::kTmpReal1));

    double *tmpReal1 = getRealData(MI::kTmpReal1);
    for (size_t i = 0; i < dimSizes.nElements(); i ++)
    {
        tmpReal1[i] = sxxSplitX[i] + sxxSplitY[i];
    }

    tmpReal1[49] = - 5;

    getTempFftwXXdx().computeR2CFft1DX(getRealMatrix(MI::kTmpReal1));
    for (size_t z = 0; z < reducedDimSizes.nz; z++)
    {
        for (size_t y = 0; y < reducedDimSizes.ny; y++)
        {
            for (size_t x = 0; x < reducedDimSizes.nx;  x++)
            {
                const size_t i = get1DIndex(z, y, x, reducedDimSizes);
                fftx[i] *= divider;
            }
        }
    }
    getTempFftwXXdx().computeC2RFft1DX(getRealMatrix(MI::kDSxxdx));
    double *dsxxdx = getRealData(MI::kDSxxdx);
    //std::cout << fftx[0 + 128] << std::endl;
    //std::cout << fftx[1 + 128] << std::endl;
    //std::cout << fftx[2 + 128] << std::endl;
#endif

#if 0
    FloatComplex *ffty = getComplexData(MI::kTmpFftwYYdy);
    const FloatComplex* ddyKShiftPos = getComplexData(MI::kDdyKShiftPos);

    const auto & dimSizes = mParameters.getFullDimensionSizes();

    const size_t nElements = mParameters.getFullDimensionSizes().nElements();

    const DimensionSizes& reducedDimSizes= mParameters.getReducedDimensionSizes();
    auto ny = dimSizes.ny;
    const double divider = 1.f / static_cast<double>(ny);

    double *p0 = getRealData(MI::kInitialPressureSourceInput);

    double* syySplitX = getRealData(MI::kSyySplitX);
    double* syySplitY = getRealData(MI::kSyySplitY);

    double* tmp = getRealData(MI::kTmpReal2);

    for (size_t i = 0; i < nElements; i++) {
        tmp[i] = -p0[i] / 2.f;
    }

    getRealMatrix(MI::kSyySplitX).copyData(getRealMatrix(MI::kTmpReal2));
    getRealMatrix(MI::kSyySplitY).copyData(getRealMatrix(MI::kTmpReal2));

    double *tmpReal2 = getRealData(MI::kTmpReal2);
    for (size_t i = 0; i < dimSizes.nElements(); i ++)
    {
        tmpReal2[i] = syySplitX[i] + syySplitY[i];
    }

    // tmpReal2[49 * 128 + 41] = - 5;

    getTempFftwYYdy().computeR2CFft1DY(getRealMatrix(MI::kTmpReal2));
    for (size_t z = 0; z < reducedDimSizes.nz; z++)
    {
        for (size_t y = 0; y < reducedDimSizes.ny; y++)
        {
            for (size_t x = 0; x < reducedDimSizes.nx;  x++)
            {
                const size_t i = get1DIndex(z, y, x, reducedDimSizes);
                ffty[i] *= ddyKShiftPos[y];
            }
        }
    }
    getTempFftwYYdy().computeC2RFft1DY(getRealMatrix(MI::kDSyydy));
    double *dsyydy = getRealData(MI::kDSyydy);
    for (size_t i = 0; i < nElements; ++ i)
        dsyydy[i] *= divider;
    //std::cout << fftx[0 + 128] << std::endl;
    //std::cout << fftx[1 + 128] << std::endl;
    //std::cout << fftx[2 + 128] << std::endl;
#endif
}


