


#ifndef __KCONFIG_INCLUDE_H__
#define __KCONFIG_INCLUDE_H__

#include <cmath>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <functional>

#include <KSpace/KMatrix.h>


struct KGrid
{

    // size_t			mDim;
    size_t mNonUniformGridFlag;

    size_t mNx;
    size_t mNy;
    size_t mNz;

    float mDx;
    float mDy;
    float mDz;

    float mDt;

    /* derivative */
    size_t mNt;

    KGrid();
};

struct Medium
{

    bool mRho0ScalarFlag;
    float mRho0Scalar;            /* kg / m^3 */

    KMatrix<float> mRho0;

    bool mC0ScalarFlag;
    bool mS0ScalarFlag;

    float mC0Scalar;                /* m/s */
    float mS0Scalar;                /* m/s */

    KMatrix<float> mC0;
    KMatrix<float> mS0;

    /* derivative */
    float mRho0SgxScalar;
    float mRho0SgyScalar;
    // float			mRho0SgzScalar;

    KMatrix<float> mRho0Sgx;
    KMatrix<float> mRho0Sgy;
    // KMatrix<float> mRho0Sgz;

    float mCref;

    Medium();
};


struct Sensor
{

    enum class SensorMaskType
    {
        kIndex = 0,
        // kCorners		= 1
    };

    SensorMaskType mSensorMaskType;
    KMatrix<size_t> mMask;

    /* derivative */
    KMatrix<size_t> mSensorMaskIndex;

    Sensor();
};


struct Source
{

    size_t mPressureSourceFlag;
    size_t mInitialPressureSourceFlag;
    KMatrix<float> mInitialPressureSourceInput;
    size_t mTransducerSourceFlag;

    /* derivative */
    size_t mVelocityXSourceFlag;
    size_t mVelocityYSourceFlag;
    // size_t			mVelocityZSourceFlag;

    Source();
};


struct KPml
{

    enum class AbsorptionType : size_t
    {
        kLossless = 0,
        // kPowerLaw		= 1,
        // kStokes			= 2
    };

    size_t mPmlXSize;
    size_t mPmlYSize;
    // size_t		mPmlZSize;

    float mPmlXAlpha;
    float mPmlYAlpha;
    // float			mPmlZAlpha;

    float mMultiAxialPmlRatio;
    AbsorptionType mAbsorbingFlag;

    KPml();
};

class KConfig
{
public:
    KConfig();

    virtual ~KConfig() = default;

    KConfig(const KConfig &) = default;

    KConfig &operator=(const KConfig &) = default;

    virtual bool preProcessing() = 0;

    virtual KConfig *clone() const = 0;

public:
//	enum class SimulationDimension {
//		k2D,
//		k3D
//	};


    size_t mElasticFlag;
    size_t mNonLinearFlag;
    bool mAxisymmetricFlag;


    KGrid mKGrid;
    Source mSource;
    Medium mMedium;
    Sensor mSensor;

    KPml mKPml;

};

class KFluidConfig : public KConfig
{
public:
    KFluidConfig(const KFluidConfig &config) = default;;

    KFluidConfig &operator=(const KFluidConfig &) = default;

    KFluidConfig() { init(); };

    ~KFluidConfig() override = default;

    void init();

    bool preProcessing() override;

    KConfig *clone() const override { return CreateInstance<KFluidConfig>(*this); };
};

class KElasticConfig : public KConfig
{
public:
    KElasticConfig(const KElasticConfig &config) = default;

    KElasticConfig &operator=(const KElasticConfig &) = default;

    KElasticConfig() { init(); };

    ~KElasticConfig() override = default;

    void init();

    bool preProcessing() override;

    KConfig *clone() const override { return CreateInstance<KElasticConfig>(*this); };
};

enum KConfigType : size_t
{
    kFluid,
    kElastic,
//	kFluid3D,
//	kElastic3D,
};


class KConfigFactory
{
public:
    KConfigFactory()
    {
        mContainer[kFluid] = AutoCreateInstance<KFluidConfig>();
        mContainer[kElastic] = AutoCreateInstance<KElasticConfig>();
    };

    std::shared_ptr<KConfig> AutoCreateKConfig(KConfigType type)
    {
        KConfig *cfg = mContainer[type]->clone();
        return std::shared_ptr<KConfig>(cfg);
    };

private:
    using ContainerType = std::unordered_map<KConfigType,
            std::shared_ptr<KConfig>, std::hash<size_t>>;

    ContainerType mContainer;
};


#endif /* ifndef __CONFIG_INCLUDE_H__ */





