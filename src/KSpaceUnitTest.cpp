


#include <gtest/gtest.h>
#include <KSpace/KSpaceSolver.h>


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

    /* set control parameters */
    auto koutput = AutoCreateInstance<KOutput>();

    std::string filename = "kspace_elastic_output.h5";
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
}


