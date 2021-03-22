


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

//    mDim = 0;
    mNonUniformGridFlag = 0;
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

Sensor::Sensor()
{
    mSensorMaskType = SensorMaskType::kIndex;
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

    mKGrid.mDt = 2e-08;

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


void
KFluidConfig::preProcessing()
{
    /**
     *  time:		mNt,
     *  density:		mRho0Sgx, mRho0Sgy
     *  speed:		mCref
     *	Sensor:		mSensorMaskIndex
     */

    /* preprocessing mNt */
    auto v = KMatrix<float>::Min(mMedium.mC0);
    auto c_min = *std::min_element(v.begin(), v.end());

    auto x_size = mKGrid.mDx * mKGrid.mNx;
    auto y_size = mKGrid.mDy * mKGrid.mNy;

    auto t_end = sqrt(x_size * x_size + y_size * y_size) / c_min;
    mKGrid.mNt = std::floor(t_end / mKGrid.mDt) + 1;

    /* preprocessing mCref */
    v = KMatrix<float>::Max(mMedium.mC0);
    mMedium.mCref = *std::max_element(v.begin(), v.end());

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
    mKGrid.mNy = 256;
    mKGrid.mNz = 1;

    mKGrid.mDx = 1e-04;
    mKGrid.mDy = 1e-04;

    mKGrid.mDt = 2e-08;

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

    mKPml.mMultiAxialPmlRatio = 0.1f;

    mKPml.mAbsorbingFlag = KPml::AbsorptionType::kLossless;
}

void
KElasticConfig::preProcessing()
{
    /**
     *  time:		mNt,
     *  density:		mRho0Sgx, mRho0Sgy
     *  speed:		mCref
     *	Sensor:		mSensorMaskIndex
     */

    /* preprocessing mNt */
    auto v = KMatrix<float>::Min(mMedium.mC0);
    auto c_min = *std::min_element(v.begin(), v.end());

    auto x_size = mKGrid.mDx * mKGrid.mNx;
    auto y_size = mKGrid.mDy * mKGrid.mNy;

    auto t_end = sqrt(x_size * x_size + y_size * y_size) / c_min;
    mKGrid.mNt = std::floor(t_end / mKGrid.mDt) + 1;

    /* preprocessing mCref */
    v = KMatrix<float>::Max(mMedium.mC0);
    mMedium.mCref = *std::max_element(v.begin(), v.end());

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
}






