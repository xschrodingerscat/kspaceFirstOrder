


#ifndef __KSPACESOLVER_INCLUDE_H__
#define __KSPACESOLVER_INCLUDE_H__

#include <memory>
#include <exception>
#include <functional>

#include <Logger/Logger.h>

#include <KSpace/KConfig.h>
#include <Parameters/Parameters.h>

#include <Containers/MatrixContainer.h>
#include <Containers/OutputStreamContainer.h>

#include <MatrixClasses/RealMatrix.h>
#include <MatrixClasses/ComplexMatrix.h>
#include <MatrixClasses/IndexMatrix.h>
#include <MatrixClasses/FftwComplexMatrix.h>
#include <MatrixClasses/FftwRealMatrix.h>

#include <OutputStreams/BaseOutputStream.h>
#include <Utils/TimeMeasure.h>

#include <KSpace/KOutput.h>
#include <KSpaceSolver/KSpaceFirstOrderSolver.h>
#include <KSpace/KVersion.h>

using KSpaceSolver = KSpaceFirstOrderSolver;

class KSpaceSolverFluid : public KSpaceSolver
{

public:
    KSpaceSolverFluid() = default;

    ~KSpaceSolverFluid() override = default;

    void allocateMemory() override;

    void loadInputData() override;

    void freeMemory() override { KSpaceSolver::freeMemory(); }

    void compute() override;

};

class KSpaceSolverElastic : public KSpaceSolver
{

public:
    KSpaceSolverElastic() = default;

    ~KSpaceSolverElastic() override = default;

    void fftwVerify();

    void allocateMemory() override;

    void loadInputData() override;

//    void freeMemory() override { KSpaceSolver::freeMemory(); }

    void compute() override;

    template<Parameters::SimulationDimension simulationDimension>
    void preProcessing();

    void generateDerivativeOperators();

    template<Parameters::SimulationDimension simulationDimension,
            bool rho0ScalarFlag,
            bool bOnAScalarFlag,
            bool c0ScalarFlag,
            bool s0ScalarFlag,
            bool alphaCoefScalarFlag>
    void computeElastic();

    template<Parameters::SimulationDimension simulationDimension>
    void computePressureGradient();

    template<Parameters::SimulationDimension simulationDimension>
    void computeSplitVelocity();

    template<Parameters::SimulationDimension simulationDimension>
    void computeVelocity();

    template<Parameters::SimulationDimension simulationDimension>
    void computeVelocityGradient();

    template<Parameters::SimulationDimension simulationDimension>
    void computeSplitPressure();

    template<Parameters::SimulationDimension simulationDimension,
            bool rho0ScalarFlag,
            bool c0ScalarFlag,
            bool s0ScalarFlag>
    void addInitialPressureSource();

    template<Parameters::SimulationDimension simulationDimension,
            bool rho0ScalarFlag,
            bool bOnAScalarFlag,
            bool c0ScalarFlag,
            bool s0ScalarFlag,
            bool alphaCoefScalarFlag>
    void computePressure();

    /* local function member */
    /// Post processing and closing the output streams.
    void postProcessing();

    /// Store sensor data.
    void storeSensorInfo();

private:
    using ComputeElasticFunc = std::function<void(KSpaceSolverElastic &)>;

    using ComputeElasticImp = std::map<
            std::tuple<Parameters::SimulationDimension,
                    bool, bool, bool, bool, bool>, ComputeElasticFunc>;

    /// Map with different implementations of ComputeMainLoop.
    static ComputeElasticImp sComputeElasticImp;

private:
    void generateLameConstant();

    void generatePml() override;

};

#endif /* ifndef __KSPACESOLVER_INCLUDE_H__ */




