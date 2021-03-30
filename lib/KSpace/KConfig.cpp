


#include <KSpace/KConfig.h>


KGrid::KGrid()
{
    mNx = 0;
    mNy = 0;
    mNz = 0;

    mDx = 0.f;
    mDy = 0.f;
    mDz = 0.f;

    mNt = 0;
    mDt = 0.f;

    mTimeEnd = 0.f;

//    mDim = 0;
    mNonUniformGridFlag = 0;
}

std::ostream &operator<<(std::ostream &os, const KGrid &grid)
{
    os << "Non uniform grid flag: " << std::boolalpha << grid.mNonUniformGridFlag << std::endl;
    os << "Nx: " << grid.mNx << std::endl;
    os << "Ny: " << grid.mNy << std::endl;
    os << "Nz: " << grid.mNz << std::endl;
    os << "Dx: " << grid.mDx << std::endl;
    os << "Dy: " << grid.mDy << std::endl;
    os << "Dz: " << grid.mDz << std::endl;
    os << "Dt: " << grid.mDt << std::endl;
    os << "Nt: " << grid.mNt << std::endl;
    return os;
}

Medium::Medium()
{
    mRho0ScalarFlag = false;
    mRho0Scalar = 0.f;

    mC0ScalarFlag = false;
    mC0Scalar = 0.f;

    mS0ScalarFlag = false;
    mS0Scalar = 0.f;

    mRho0SgxScalar = 0.f;
    mRho0SgyScalar = 0.f;
//    mRho0SgzScalar = 0.f;

    mCref = 0;
}

std::ostream &operator<<(std::ostream &os, const Medium &medium)
{
    os << "Rho0 (density): ";
    medium.mRho0ScalarFlag ? (os << "scalar " << medium.mRho0Scalar)
                           : (os << "matrix " << medium.mRho0.rowSize() << "*" << medium.mRho0.colSize());
    os << std::endl;

    os << "C0 (compression speed): ";
    medium.mC0ScalarFlag ? (os << "scalar " << medium.mC0Scalar)
                         : (os << "matrix " << medium.mC0.rowSize() << "*" << medium.mC0.colSize());
    os << std::endl;

    os << "S0 (shear speed): ";
    medium.mS0ScalarFlag ? (os << "scalar " << medium.mS0Scalar)
                         : (os << "matrix " << medium.mS0.rowSize() << "*" << medium.mS0.colSize());
    os << std::endl;

    return os;
}

Sensor::Sensor()
{
    mSensorMaskType = SensorMaskType::kIndex;

}

std::ostream &operator<<(std::ostream &os, const Sensor &sensor)
{
    using SMT = Sensor::SensorMaskType;

    os << "Sensor mask type: ";
    os << ((sensor.mSensorMaskType == SMT::kIndex) ? "Index" : "Unknown") << std::endl;
    os << "Mask matrix: ";
    os << sensor.mMask.rowSize() << "*" << sensor.mMask.colSize() << std::endl;
    os << "Sensor mask index matrix: ";
    os << sensor.mSensorMaskIndex.rowSize() << "*" << sensor.mSensorMaskIndex.colSize() << std::endl;

    return os;
}

Source::Source()
{
    mPressureSourceFlag = 0;
    mInitialPressureSourceFlag = 0;
    mTransducerSourceFlag = 0;

    mVelocityXSourceFlag = 0;
    mVelocityYSourceFlag = 0;
//    mVelocityZSourceFlag = 0;
}

std::ostream &operator<<(std::ostream &os, const Source &source)
{
    os << "Pressure source flag: " << std::boolalpha << source.mPressureSourceFlag << std::endl;
    os << "Initial pressure source flag: " << std::boolalpha << source.mInitialPressureSourceFlag << std::endl;
    os << "Initial pressure source input matrix: "
       << source.mInitialPressureSourceInput.rowSize() << "*"
       << source.mInitialPressureSourceInput.colSize() << std::endl;

    return os;
}

KPml::KPml()
{
    mPmlXSize = 0;
    mPmlYSize = 0;
//    mPmlZSize = 0;

    mPmlXAlpha = 0.f;
    mPmlYAlpha = 0.f;
//    mPmlZAlpha = 0.f;

    mMultiAxialPmlRatio = 0.f;
    mAbsorbingFlag = AbsorptionType::kLossless;
}

std::ostream &operator<<(std::ostream &os, const KPml &pml)
{
    using AT = KPml::AbsorptionType;

    os << "Pml x size: " << pml.mPmlXSize << std::endl;
    os << "Pml y size: " << pml.mPmlYSize << std::endl;
    os << "Pml x alpha: " << pml.mPmlXAlpha << std::endl;
    os << "Pml y Alpha: " << pml.mPmlYAlpha << std::endl;
    os << "Multi axial pml ratio: " << pml.mMultiAxialPmlRatio << std::endl;
    os << "Absorbing flag: " << (pml.mAbsorbingFlag == AT::kLossless ? "kLossless" : "Unknown") << std::endl;

    return os;
}


KConfig::KConfig() :
        mKGrid(),
        mSource(),
        mMedium(),
        mSensor(),
        mKPml()
{
    mElasticFlag = 0;
    mNonLinearFlag = 0;
    mAxisymmetricFlag = false;

}

std::ostream &operator<<(std::ostream &os, const KConfig &config)
{
    os << "---------------------------------------------------------" << std::endl;
    os << "Elastic flag: " << std::boolalpha << config.mElasticFlag << std::endl;
    os << "NonLinear flag: " << std::boolalpha << config.mNonLinearFlag << std::endl;
    os << "Axisymmetric flag: " << std::boolalpha << config.mAxisymmetricFlag << std::endl;
    os << "---------------------------------------------------------" << std::endl;
    os << "Grid: " << std::endl;
    os << config.mKGrid << std::endl;
    os << "---------------------------------------------------------" << std::endl;
    os << "Source: " << std::endl;
    os << config.mSource << std::endl;
    os << "---------------------------------------------------------" << std::endl;
    os << "Medium: " << std::endl;
    os << config.mMedium << std::endl;
    os << "---------------------------------------------------------" << std::endl;
    os << "Sensor: " << std::endl;
    os << config.mSensor << std::endl;
    os << "---------------------------------------------------------" << std::endl;
    os << "KPml: " << std::endl;
    os << config.mKPml << std::endl;
    os << "---------------------------------------------------------" << std::endl;

    return os;
}

void
KFluidConfig::init()
{
    /* Global */
    mElasticFlag = 0;
    mNonLinearFlag = 0;
    mAxisymmetricFlag = false;

    /* Grid */
//    mKGrid.mDim = 2;
    mKGrid.mNonUniformGridFlag = 0;

    mKGrid.mNx = 128;
    mKGrid.mNy = 256;
    mKGrid.mNz = 1;

    mKGrid.mDx = 1e-04;
    mKGrid.mDy = 1e-04;

    mKGrid.mTimeEnd = 1e-5;
    /* temp */
    auto Nx = mKGrid.mNx;
    auto Ny = mKGrid.mNy;

    /* Medium */
    auto density = 1000;
    mMedium.mRho0ScalarFlag = false;
    mMedium.mRho0 = density * KMatrix<float>::Ones(Nx, Ny);

    auto sound_speed = 1500;
    mMedium.mC0ScalarFlag = false;
    mMedium.mC0 = sound_speed * KMatrix<float>::Ones(Nx, Ny);

    /* Sensor */
    mSensor.mSensorMaskType = Sensor::SensorMaskType::kIndex;

    auto row = 20;
    auto right = Ny - 1;
    mSensor.mMask = KMatrix<size_t>::Line(Nx, Ny, row, 0, right);

    /* Source */
    mSource.mPressureSourceFlag = 0;
    mSource.mTransducerSourceFlag = 0;

    mSource.mInitialPressureSourceFlag = 1;

    auto disc_magnitude = 5;
    auto disc_x_pos = 50 - 1;
    auto disc_y_pos = 50 - 1;
    auto disc_radius = 8;

    mSource.mInitialPressureSourceInput = disc_magnitude
                                          * KMatrix<float>::Disc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

    mSource.mVelocityXSourceFlag = 0;
    mSource.mVelocityYSourceFlag = 0;

    /* Pml */
    mKPml.mPmlXSize = 20;
    mKPml.mPmlYSize = 20;

    mKPml.mPmlXAlpha = 2.0;
    mKPml.mPmlYAlpha = 2.0;

    mKPml.mAbsorbingFlag = KPml::AbsorptionType::kLossless;
}


bool
KFluidConfig::preProcessing()
{
    /**
     *  time:		mNt,
     *  density:		mRho0Sgx, mRho0Sgy
     *  speed:		mCref
     *	Sensor:		mSensorMaskIndex
     */


    /* preprocessing mCref */
    auto v = KMatrix<float>::Min(mMedium.mC0);
    v = KMatrix<float>::Max(mMedium.mC0);
    mMedium.mCref = *std::max_element(v.begin(), v.end());

    /* preprocessing mDt and mNt */
    auto min_grid = std::min(mKGrid.mDx, mKGrid.mDy);

    const auto KSPACE_CFL = 0.3f;
    mKGrid.mDt = KSPACE_CFL * min_grid / mMedium.mCref;
    mKGrid.mNt = std::floor(mKGrid.mTimeEnd / mKGrid.mDt) + 1;

    /* preprocessing mSensorMaskIndex */
    /* column major                    */
    auto &mask = mSensor.mMask;
    auto &maskIndex = mSensor.mSensorMaskIndex;

    auto base = KMatrixType<size_t>(1, std::vector<size_t>(0));

    auto col = mask.colSize();
    auto row = mask.rowSize();

    for (size_t i = 0; i < row; ++i)
        for (size_t j = 0; j < col; ++j)
            if (mask[i][j] == 1)
            {
                size_t index = j * row + i + 1;
                base[0].push_back(index);
            }

    maskIndex = KMatrix<size_t>(base);

    /* preprocessing  density: mRho0Sgx, mRho0Sgy  */
    mMedium.mRho0Sgx = mMedium.mRho0;
    mMedium.mRho0Sgy = mMedium.mRho0;

    return true;
}

void
KElasticConfig::init()
{
    /* Global */
    mElasticFlag = 1;
    mNonLinearFlag = 0;
    mAxisymmetricFlag = false;

    /* Grid */
//    mKGrid.mDim = 2;
    mKGrid.mNonUniformGridFlag = 0;

    mKGrid.mNx = 128;
    mKGrid.mNy = 128;
    mKGrid.mNz = 1;

    mKGrid.mDx = 1e-04;
    mKGrid.mDy = 1e-04;

    mKGrid.mTimeEnd = 1e-5;

    /* temp */
    auto Nx = mKGrid.mNx;
    auto Ny = mKGrid.mNy;

    /* Medium */
    auto density = 1000;
    mMedium.mRho0ScalarFlag = false;
    mMedium.mRho0 = density * KMatrix<float>::Ones(Nx, Ny);

    auto cs_sound_speed = 1500;
    mMedium.mC0ScalarFlag = false;
    mMedium.mC0 = cs_sound_speed * KMatrix<float>::Ones(Nx, Ny);

    auto ss_sound_speed = 500;
    mMedium.mS0ScalarFlag = false;
    mMedium.mS0 = ss_sound_speed * KMatrix<float>::Ones(Nx, Ny);

    /* Sensor */
    mSensor.mSensorMaskType = Sensor::SensorMaskType::kIndex;

    auto row = 20;
    auto right = Ny - 1;
    mSensor.mMask = KMatrix<size_t>::Line(Nx, Ny, row, 0, right);

    /* Source */
    mSource.mPressureSourceFlag = 0;
    mSource.mTransducerSourceFlag = 0;

    mSource.mInitialPressureSourceFlag = 1;

    auto line_magnitude = 5;
    auto line_y = 29;
    auto line_start = 29;
    auto line_end = 97;

    auto initPressure = line_magnitude * KMatrix<float>::Line(Nx, Ny, line_y, line_start, line_end);
    mSource.mInitialPressureSourceInput = initPressure;

    mSource.mVelocityXSourceFlag = 0;
    mSource.mVelocityYSourceFlag = 0;

    /* Pml */
    mKPml.mPmlXSize = 20;
    mKPml.mPmlYSize = 20;

    mKPml.mPmlXAlpha = 2.0;
    mKPml.mPmlYAlpha = 2.0;

    mKPml.mMultiAxialPmlRatio = 0.1f;

    mKPml.mAbsorbingFlag = KPml::AbsorptionType::kLossless;

}

bool
KElasticConfig::preProcessing()
{
    /**
     *  time:		mNt,
     *  density:		mRho0Sgx, mRho0Sgy
     *  speed:		mCref
     *	Sensor:		mSensorMaskIndex
     */

    /* preprocessing mCref */
    auto v = KMatrix<float>::Min(mMedium.mC0);
    v = KMatrix<float>::Max(mMedium.mC0);
    mMedium.mCref = *std::max_element(v.begin(), v.end());

    /* preprocessing mDt and mNt */
    auto min_grid = std::min(mKGrid.mDx, mKGrid.mDy);

    const auto KSPACE_CFL = 0.3f;
    mKGrid.mDt = KSPACE_CFL * min_grid / mMedium.mCref;
    mKGrid.mNt = std::floor(mKGrid.mTimeEnd / mKGrid.mDt) + 1;

    /* preprocessing mSensorMaskIndex */
    /* column major                    */
    auto &mask = mSensor.mMask;
    auto &maskIndex = mSensor.mSensorMaskIndex;

    auto base = KMatrixType<size_t>(1, std::vector<size_t>(0));

    auto col = mask.colSize();
    auto row = mask.rowSize();

    for (size_t i = 0; i < row; ++i)
        for (size_t j = 0; j < col; ++j)
            if (mask[i][j] == 1)
            {
                size_t index = j * row + i + 1;
                base[0].push_back(index);
            }

    maskIndex = KMatrix<size_t>(base);

    /* preprocessing  density: mRho0Sgx, mRho0Sgy  */
    mMedium.mRho0Sgx = mMedium.mRho0;
    mMedium.mRho0Sgy = mMedium.mRho0;

    return true;
}






