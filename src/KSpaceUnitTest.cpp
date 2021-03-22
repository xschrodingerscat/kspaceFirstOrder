



#include <gtest/gtest.h>
#include <KSpace/KSpaceSolver.h>


TEST(KSpace, kFluid)
{
	/* create special global configuration (default 2D)*/
	auto factory = AutoCreateInstance<KConfigFactory>();
	auto kcfg = factory->AutoCreateKConfig(kFluid);

	auto kcmds = AutoCreateInstance<KCmds>();
	/* preprocessing */
	kcfg->preProcessing();

	/* initialize praram object with preprocessed config */
	auto& params = Parameters::getInstance();
	params.init(*kcfg, *kcmds);

	/* create first order PDE solver*/
	auto solver = AutoCreateInstance<KSpaceSolverFluid>();

	/* allocate memory for its container (such as matrix) */
	solver->allocateMemory();

	/* load data from kconfig */
	solver->loadInputData();

	/* solve the pdes system*/
	solver->compute();

}

TEST(KSpace, kElastic)
{
    /* create special global configuration (default 2D)*/
    auto factory = AutoCreateInstance<KConfigFactory>();
    auto kcfg = factory->AutoCreateKConfig(kElastic);

    /* preprocessing */
    kcfg->preProcessing();

    auto kcmds = AutoCreateInstance<KCmds>();

    /* initialize praram object with preprocessed config */
    auto& params = Parameters::getInstance();
    params.init(*kcfg, *kcmds);

    /* create first order PDE solver*/
    auto solver = AutoCreateInstance<KSpaceSolverElastic>();

    /* allocate memory for its container (such as matrix) */
    solver->allocateMemory();

    /* load data from kconfig */
    solver->loadInputData();

    /* solve the pdes system*/
    solver->compute();

}

TEST(KSpace, kFftXY)
{
	/* create special global configuration (default 2D)*/
	auto factory = AutoCreateInstance<KConfigFactory>();
	auto kcfg = factory->AutoCreateKConfig(kElastic);

	/* preprocessing */
	kcfg->preProcessing();

	auto kcmds = AutoCreateInstance<KCmds>();

	/* initialize praram object with preprocessed config */
	auto& params = Parameters::getInstance();
	params.init(*kcfg, *kcmds);

	/* create first order PDE solver*/
	auto solver = AutoCreateInstance<KSpaceSolverElastic>();

	/* allocate memory for its container (such as matrix) */
	solver->allocateMemory();

	solver->loadInputData();

	solver->fftTest();
}



