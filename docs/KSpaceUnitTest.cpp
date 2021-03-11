
#include <exception>

#include <KSpaceSolver/KSpaceFirstOrderSolver.h>
#include <Logger/Logger.h>
#include <gtest/gtest.h>

using std::string;

TEST(KSpace, Trivial01)
{
	// Create k-Space solver
	KSpaceFirstOrderSolver kSpaceSolver;

	// Print header
	Logger::log(Logger::LogLevel::kBasic, kOutFmtFirstSeparator);
	Logger::log(Logger::LogLevel::kBasic, kOutFmtCodeName, kSpaceSolver.getCodeName().c_str());
	Logger::log(Logger::LogLevel::kBasic, kOutFmtSeparator);

	// Create parameters and parse command line
	Parameters& parameters = Parameters::getInstance();

	//-------------------------------------------- Initialize simulation -----------------------------------------------//
	int argc = 1;
	char *argv[1] = {"KSpace"};
	try
	{
		// Initialize parameters by parsing the command line and reading input file scalars
		parameters.init(argc, argv);

		// When we know the GPU, we can print out the code version
		if (parameters.isPrintVersionOnly())
		{
			kSpaceSolver.printFullCodeNameAndLicense();
			exit(EXIT_SUCCESS);
		}
	}
	catch (const std::exception& e)
	{
		Logger::log(Logger::LogLevel::kBasic, kOutFmtFailed);
		Logger::log(Logger::LogLevel::kBasic, kOutFmtLastSeparator);
		Logger::errorAndTerminate(Logger::wordWrapString(e.what(),kErrFmtPathDelimiters, 9));
	}

	// Print simulation setup
	parameters.printSimulatoinSetup();

	Logger::log(Logger::LogLevel::kBasic, kOutFmtInitializationHeader);

	//------------------------------------------------ Allocate memory -------------------------------------------------//
	try
	{
		kSpaceSolver.allocateMemory();
	}
	catch (const std::bad_alloc& e)
	{
		Logger::log(Logger::LogLevel::kBasic, kOutFmtFailed);
		Logger::log(Logger::LogLevel::kBasic, kOutFmtLastSeparator);
		// 9 = Indentation of Error:
		Logger::errorAndTerminate(Logger::wordWrapString(kErrFmtOutOfMemory," ", 9));
	}
	catch (const std::exception& e)
	{
		Logger::log(Logger::LogLevel::kBasic, kOutFmtFailed);
		Logger::log(Logger::LogLevel::kBasic, kOutFmtLastSeparator);

		const string errorMessage = string(kErrFmtUnknownError) + e.what();
		Logger::errorAndTerminate(Logger::wordWrapString(errorMessage, kErrFmtPathDelimiters, 9));
	}

	//------------------------------------------------ Load input data -------------------------------------------------//
	try
	{
		kSpaceSolver.loadInputData();
	}
	catch (const std::ios::failure& e)
	{
		Logger::log(Logger::LogLevel::kBasic, kOutFmtFailed);
		Logger::log(Logger::LogLevel::kBasic, kOutFmtLastSeparator);
		// 9 = Indentation of Error:
		Logger::errorAndTerminate(Logger::wordWrapString(e.what(), kErrFmtPathDelimiters, 9));
	}
	catch (const std::exception& e)
	{
		Logger::log(Logger::LogLevel::kBasic, kOutFmtFailed);
		Logger::log(Logger::LogLevel::kBasic, kOutFmtLastSeparator);

		const string errorMessage = string(kErrFmtUnknownError) + e.what();
		Logger::errorAndTerminate(Logger::wordWrapString(errorMessage, kErrFmtPathDelimiters, 9));
	}

	// Did we recover from checkpoint?
	if (parameters.getTimeIndex() > 0)
	{
		Logger::log(Logger::LogLevel::kBasic, kOutFmtSeparator);
		Logger::log(Logger::LogLevel::kBasic, kOutFmtRecoveredFrom, parameters.getTimeIndex());
	}

	//-------------------------------------------------- Simulation ----------------------------------------------------//
	Logger::log(Logger::LogLevel::kBasic, kOutFmtSeparator);
	// Exception are caught inside due to different log formats
	kSpaceSolver.compute();

	//----------------------------------------------------- Summary ----------------------------------------------------//
	Logger::log(Logger::LogLevel::kBasic, kOutFmtSummaryHeader);
	Logger::log(Logger::LogLevel::kBasic, kOutFmtMemoryUsage, kSpaceSolver.getMemoryUsage());

	Logger::log(Logger::LogLevel::kBasic, kOutFmtSeparator);

	// Elapsed Time time
	if (kSpaceSolver.getCumulatedTotalTime() != kSpaceSolver.getTotalTime())
	{
		Logger::log(Logger::LogLevel::kBasic, kOutFmtLegExecutionTime, kSpaceSolver.getTotalTime());
	}
	Logger::log(Logger::LogLevel::kBasic, kOutFmtTotalExecutionTime, kSpaceSolver.getCumulatedTotalTime());

	Logger::log(Logger::LogLevel::kBasic, kOutFmtEndOfSimulation);

}// end of main

//----------------------------------------------------------------------------------------------------------------------




