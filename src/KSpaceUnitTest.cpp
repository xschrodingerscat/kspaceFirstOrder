



#include <gtest/gtest.h>

#include <KSpace/KSpaceSolver.h>


TEST(KSpace, kFluid)
{
	auto factory = AutoCreateInstance<KConfigFactory>();
	auto kcfg = factory->CreateKConfig(kFluid);

	kcfg->preProcessing();

	auto& params = Parameters::getInstance();
	params.init(*kcfg);


	KSpaceFirstOrderSolver kSpaceSolver;

	kSpaceSolver.allocateMemory();
}


TEST(KSpace, kElastic)
{

}




