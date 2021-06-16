


#include <gtest/gtest.h>
#include <KSpace/KSpaceSolver.h>

TEST(KSpace, kFluid)
{
    /* version checking */
    const auto version = getVersion();
    std::cout << "The current version of the KSpaceSolver: " + version << std::endl;

    /* create special global configuration (default 2D)*/
    auto factory = AutoCreateInstance<KConfigFactory>();
    auto kcfg = factory->AutoCreateKConfig(kElastic);

    /* set fluid (default elastic) */
    kcfg->mMedium.mS0ScalarFlag = false;
    kcfg->mMedium.mS0 = KMatrix<double>::Zero(kcfg->mKGrid.mNx, kcfg->mKGrid.mNy);

    /* preprocessing */
    kcfg->preProcessing();

    std::cout << *kcfg << std::endl;

    /* set control parameters */
    auto koutput = AutoCreateInstance<KOutput>();

    std::string filename = "kspace_fluid_output.h5";
    koutput->setOutputFileName(filename);
    koutput->setOutputToFileFlag(true);

    /* initialize param object with preprocessed config */
    auto &params = Parameters::getInstance();
    params.init(*kcfg, *koutput);

    /* create first order PDE solver*/
    auto solver = AutoCreateInstance<KSpaceSolverElastic>();

    /* allocate memory for its container (such as matrix) */
    solver->allocateMemory();

    /* load data from kconfig */
    solver->loadInputData();

    /* solve the PDEs system*/
    solver->compute();

    // solver->freeMemory();

    /* results. note that : the matrix with column major order. */
    auto & result = koutput->mPressureRaw;
    EXPECT_EQ(result.rowSize(), kcfg->mKGrid.mNx);
    EXPECT_EQ(result.colSize(), kcfg->mKGrid.mNt);

    /* show some elements in the result matrix for validating */
//    std::cout << result[0][0] << std::endl;
//    std::cout << result[0][1] << std::endl;
//    std::cout << result[0][2] << std::endl;
//    std::cout << result[0][3] << std::endl;
}

TEST(KSpace, kElastic)
{
    /* version checking */
    const auto version = getVersion();
    std::cout << "The current version of the KSpaceSolver: " + version << std::endl;

    /* create special global configuration (default 2D)*/
    auto factory = AutoCreateInstance<KConfigFactory>();
    auto kcfg = factory->AutoCreateKConfig(kElastic);

    /* preprocessing */
    kcfg->preProcessing();

    std::cout << *kcfg << std::endl;

    /* set control parameters */
    auto koutput = AutoCreateInstance<KOutput>();

    std::string filename = "kspace_elastic_output.h5";
    koutput->setOutputFileName(filename);
    koutput->setOutputToFileFlag(false);

    /* initialize param object with preprocessed config */
    auto &params = Parameters::getInstance();
    params.init(*kcfg, *koutput);

    /* create first order PDE solver*/
    auto solver = AutoCreateInstance<KSpaceSolverElastic>();

    /* allocate memory for its container (such as matrix) */
    solver->allocateMemory();

    /* load data from kconfig */
    solver->loadInputData();

    /* solve the PDEs system*/
    solver->compute();

    // solver->freeMemory();

    /* results. note that : the matrix with column major order. */
    auto & result = koutput->mPressureRaw;
    EXPECT_EQ(result.rowSize(), kcfg->mKGrid.mNx);
    EXPECT_EQ(result.colSize(), kcfg->mKGrid.mNt);
}

TEST(KSpace, kMixture)
{
    /* version checking */
    const auto version = getVersion();
    std::cout << "The current version of the KSpaceSolver: " + version << std::endl;

    /* create special global configuration (default 2D)*/
    auto factory = AutoCreateInstance<KConfigFactory>();
    auto kcfg = factory->AutoCreateKConfig(kElastic);

    /* set mixture */
    auto Nx = 128;
    auto Ny = 128;

    kcfg->mKGrid.mNx = Nx;
    kcfg->mKGrid.mNy = Ny;

    /* wedge mask */
    auto r1_mask = KMatrix<double>::Rect(Nx, Ny, 0, 0, Ny - 1, Nx / 2 - 1);
    /* part mask */
    auto r2_mask = KMatrix<double>::Rect(Nx, Ny, 0, Nx / 2, Ny - 1, Nx - 1);

    /* wedge */
    auto c1_magnitude = 2700.f;
    auto s1_magnitude = 1100.f;
    auto rho1 = 1150.f;

    /* part */
    auto c2_magnitude = 5800.f;
    auto s2_magnitude = 3200.f;
    auto rho2 = 7800.f;

    /* medium */
    kcfg->mMedium.mC0 = c1_magnitude * r1_mask + c2_magnitude * r2_mask;
    kcfg->mMedium.mS0 = s1_magnitude * r1_mask + s2_magnitude * r2_mask;
    kcfg->mMedium.mRho0 = rho1 * r1_mask + rho2 * r2_mask;

    /* preprocessing */
    kcfg->preProcessing();

    std::cout << *kcfg << std::endl;

    /* set control parameters */
    auto koutput = AutoCreateInstance<KOutput>();

    std::string filename = "kspace_mixture_output.h5";
    koutput->setOutputFileName(filename);
    koutput->setOutputToFileFlag(false);

    /* initialize param object with preprocessed config */
    auto &params = Parameters::getInstance();
    params.init(*kcfg, *koutput);

    /* create first order PDE solver*/
    auto solver = AutoCreateInstance<KSpaceSolverElastic>();

    /* allocate memory for its container (such as matrix) */
    solver->allocateMemory();

    /* load data from kconfig */
    solver->loadInputData();

    /* solve the PDEs system*/
    solver->compute();

    // solver->freeMemory();

    /* results. note that : the matrix with column major order. */
    auto & result = koutput->mPressureRaw;
    EXPECT_EQ(result.rowSize(), kcfg->mKGrid.mNx);
    EXPECT_EQ(result.colSize(), kcfg->mKGrid.mNt);
}

TEST(KSpace, kMixtureWithFlaws)
{
    /* version checking */
    const auto version = getVersion();
    std::cout << "The current version of the KSpaceSolver: " + version << std::endl;

    /* create special global configuration (default 2D)*/
    auto factory = AutoCreateInstance<KConfigFactory>();
    auto kcfg = factory->AutoCreateKConfig(kElastic);

    /* set mixture */
    auto Nx = 128;
    auto Ny = 128;

    kcfg->mKGrid.mNx = Nx;
    kcfg->mKGrid.mNy = Ny;

    /* wedge mask */
    auto r1_mask = KMatrix<double>::Rect(Nx, Ny, 0, 0, Ny - 1, Nx / 2 - 1);
    /* part mask */
    auto r2_mask = KMatrix<double>::Rect(Nx, Ny, 0, Nx / 2, Ny - 1, Nx - 1);

    /* flaws definition */
    class PrdFillDisc
    {
    public:
        PrdFillDisc(size_t x, size_t y, size_t radius) : mX(x), mY(y), mRadius(radius), mVal(0) {}
        void setVal(double val) { mVal = val; }
        void operator()(double &e, size_t i, size_t j)
        {
            if ((mX - i) * (mX - i) + (mY - j) * (mY - j) < mRadius * mRadius)
                e = mVal;
        }

    private:
        size_t mX;
        size_t mY;
        double mVal;
        double mRadius;
    };

    /* wedge */
    auto c1_magnitude = 2700.;
    auto s1_magnitude = 1100.;
    auto rho1 = 1150.;

    /* part */
    auto c2_magnitude = 5800.;
    auto s2_magnitude = 3200.;
    auto rho2 = 7800.;

    /* flaws */
    auto x_pos = 89;
    auto y_pos = 63;
    auto radius = 10;
    auto flaw_c0_magnitude = 340.;
    auto flaw_s0_magnitude = 10.;
    auto flaw_rho0 = 1.29;

    auto fillDisc = PrdFillDisc(x_pos, y_pos, radius);

    /* medium */
    kcfg->mMedium.mC0 = c1_magnitude * r1_mask + c2_magnitude * r2_mask;
    fillDisc.setVal(flaw_c0_magnitude);
    KMatrix<double>::Fill(kcfg->mMedium.mC0, fillDisc);

    kcfg->mMedium.mS0 = s1_magnitude * r1_mask + s2_magnitude * r2_mask;
    fillDisc.setVal(flaw_s0_magnitude);
    KMatrix<double>::Fill(kcfg->mMedium.mS0, fillDisc);

    kcfg->mMedium.mRho0 = rho1 * r1_mask + rho2 * r2_mask;
    fillDisc.setVal(flaw_rho0);
    KMatrix<double>::Fill(kcfg->mMedium.mRho0, fillDisc);

    /* preprocessing */
    kcfg->preProcessing();

    std::cout << *kcfg << std::endl;

    /* set control parameters */
    auto koutput = AutoCreateInstance<KOutput>();

    /* save output data to hdf5file */
    std::string filename = "kspace_mixture_with_flaws.h5";
    koutput->setOutputFileName(filename);
    koutput->setOutputToFileFlag(true);

    /* initialize param object with preprocessed config */
    auto &params = Parameters::getInstance();
    params.init(*kcfg, *koutput);

    /* create first order PDE solver*/
    auto solver = AutoCreateInstance<KSpaceSolverElastic>();

    /* allocate memory for its container (such as matrix) */
    solver->allocateMemory();

    /* load data from kconfig */
    solver->loadInputData();

    /* solve the PDEs system*/
    solver->compute();

    // solver->freeMemory();

    /* results. note that : the matrix with column major order. */
    auto & result = koutput->mPressureRaw;
    EXPECT_EQ(result.rowSize(), kcfg->mKGrid.mNx);
    EXPECT_EQ(result.colSize(), kcfg->mKGrid.mNt);
}


TEST(KSpace, kMiscVerify)
{
    /* version checking */
    const auto version = getVersion();
    std::cout << "The current version of the KSpaceSolver: " + version << std::endl;

    /* create special global configuration (default 2D)*/
    auto factory = AutoCreateInstance<KConfigFactory>();
    auto kcfg = factory->AutoCreateKConfig(kElastic);

    /* set mixture */
    auto Nx = 8;
    auto Ny = 8;

    kcfg->mKGrid.mNx = Nx;
    kcfg->mKGrid.mNy = Ny;

    /* wedge mask */
    auto r1_mask = KMatrix<double>::Rect(Nx, Ny, 0, 0, Ny - 1, Nx / 2 - 1);
    /* part mask */
    auto r2_mask = KMatrix<double>::Rect(Nx, Ny, 0, Nx / 2, Ny - 1, Nx - 1);

    /* wedge */
    auto c1_magnitude = 2700.f;
    auto s1_magnitude = 1100.f;
    auto rho1 = 1150.f;

    /* part */
    auto c2_magnitude = 5800.f;
    auto s2_magnitude = 3200.f;
    auto rho2 = 7800.f;

    /* medium */
    kcfg->mMedium.mC0 = c1_magnitude * r1_mask + c2_magnitude * r2_mask;
    kcfg->mMedium.mS0 = s1_magnitude * r1_mask + s2_magnitude * r2_mask;
    kcfg->mMedium.mRho0 = rho1 * r1_mask + rho2 * r2_mask;

    auto line_magnitude = 5;
    auto line_y = 3;
    auto line_start = 0;
    auto line_end = Nx - 1;

    auto initPressure = line_magnitude * KMatrix<double>::Line(Nx, Ny, line_y, line_start, line_end);
    kcfg->mSource.mInitialPressureSourceInput = initPressure;
    /* preprocessing */
    kcfg->preProcessing();

    std::cout << *kcfg << std::endl;

    /* set control parameters */
    auto koutput = AutoCreateInstance<KOutput>();

    /* initialize param object with preprocessed config */
    auto &params = Parameters::getInstance();
    params.init(*kcfg, *koutput);

    /* create first order PDE solver*/
    auto solver = AutoCreateInstance<KSpaceSolverElastic>();

    /* allocate memory for its container (such as matrix) */
    solver->allocateMemory();

    /* load data from kconfig */
    solver->loadInputData();


    solver->miscVerify();
}

