


#include <KSpace/KSpaceSolver.h>

using std::ios;
/* Shortcut for Simulation dimensions. */
using SD = Parameters::SimulationDimension;
/* Shortcut for Matrix id in the container. */
using MI = MatrixContainer::MatrixIdx;
/* Shortcut for Output stream id in the container. */
using OI = OutputStreamContainer::OutputStreamIdx;

void
KSpaceSolverFluid::allocateMemory()
{
    // Add matrices into the container and create all matrices
    mMatrixContainer.init();
    mMatrixContainer.createMatrices();

    // Add output streams into container
    mOutputStreamContainer.init(mMatrixContainer);
}

void
KSpaceSolverFluid::loadInputData()
{
    // Load data from the input file
    mMatrixContainer.loadDataFromKConfig();

    Hdf5File &outputFile = mParameters.getOutputFile();

    auto & ksampler = mParameters.getKOutput();
    auto filename = ksampler.getOutputFileName();

    if (!filename.empty() && !outputFile.canAccess(filename))
        outputFile.create(filename);

    mOutputStreamContainer.createStreams();
}

void
KSpaceSolverFluid::compute()
{
    /* Initialize all used FFTW plans */
    initializeFftwPlans();

    preProcessing<SD::k2D>();

    sComputeMainLoopImp[std::make_tuple(mParameters.getSimulationDimension(),
                                        mParameters.getRho0ScalarFlag(),
                                        mParameters.getBOnAScalarFlag(),
                                        mParameters.getC0ScalarFlag(),
                                        mParameters.getAlphaCoeffScalarFlag())](*this);

    /* Post processing phase */
    postProcessing();

    writeOutputDataInfo();

    mParameters.getOutputFile().close();
}




