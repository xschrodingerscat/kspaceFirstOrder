


#ifndef __KCONFIG_INCLUDE_H__
#define __KCONFIG_INCLUDE_H__ 

#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <functional>

#include <KSpace/KMatrix.h>

template<typename T,  typename... Args>
T* CreateInstance(Args... args)
{
	return new T(args...);
}

template<typename T,  typename... Args>
std::shared_ptr<T> AutoCreateInstance(Args... args)
{
	return std::make_shared<T>(args...);
}


struct KGrid {

	size_t			mDim;
    size_t			mNonUniformGridFlag;

	size_t			mNx;
	size_t			mNy;
	size_t			mNz;

	double			mDx;
	double			mDy;
	double			mDz;

	double			mDt;

	/* derivate */
	size_t			mNt;

	KGrid();
};

struct Medium {

    bool			mRho0ScalarFlag;
    float			mRho0Scalar;			/* kg / m^3 */

	KMatrix<double> mRho0;

    bool			mC0ScalarFlag;
    float			mC0Scalar;				/* m/s */

	KMatrix<double> mC0;

	/* derivate */
    float			mRho0SgxScalar;
    float			mRho0SgyScalar;
    float			mRho0SgzScalar;

	KMatrix<double> mRho0Sgx;
	KMatrix<double> mRho0Sgy;
	KMatrix<double> mRho0Sgz;
	
	double			mCref;

	Medium();
};


struct Sensor {

	enum class SensorMaskType {
		kIndex			= 0,
		kCorners		= 1
	};

    SensorMaskType	mSensorMaskType;
	KMatrix<size_t> mMask;

	/* derivate */
	std::vector<size_t> mSensorMaskIndex;

	Sensor();
};


struct Source {

    size_t			mPressureSourceFlag;
    size_t			mInitialPressureSourceFlag;
	KMatrix<double> mInitialPressureSourceInput;
    size_t			mTransducerSourceFlag;
	
	/* derivate */
	size_t			mVelocityXSourceFlag;
	size_t			mVelocityYSourceFlag;
	size_t			mVelocityZSourceFlag;

	Source();
};


struct KPml {

	enum class AbsorptionType: size_t {
		kLossless		= 0,
		kPowerLaw		= 1,
		kStokes			= 2
	};

	size_t			mPmlXSize;
	size_t			mPmlYSize;
	size_t			mPmlZSize;
	double			mPmlXAlpha;
	double			mPmlYAlpha;
	double			mPmlZAlpha;

    AbsorptionType	mAbsorbingFlag;

	KPml();
};

class KConfig {
public:
	KConfig();
	virtual ~KConfig() {};

	KConfig(const KConfig &) = default;
	KConfig& operator=(const KConfig &) = default;

	virtual void preProcessing() = 0;
	virtual KConfig *clone() const = 0;

public:
	enum class SimulationDimension {
		k2D,
		k3D
	};

	size_t			mElasticFlag;
    size_t			mNonLinearFlag;
    bool			mAxisymmetricFlag;


	KGrid			mKGrid;
	Source			mSource;
	Medium			mMedium;
	Sensor			mSensor;

	KPml			mKPml;

};

class KFluidConfig : public KConfig {
public:
	KFluidConfig(const KFluidConfig &) = default; 
	KFluidConfig() { init(); };

	void init();

	void preProcessing()  override;
	KConfig *clone() const override { return CreateInstance<KFluidConfig>(*this); };
};

class KElasticConfig: public KConfig {
public:
	KElasticConfig(const KElasticConfig &) = default; 
	KElasticConfig() {init(); };

	void init();

	void preProcessing()  override;
	KConfig *clone() const override { return CreateInstance<KElasticConfig>(*this); };
};

enum KConfigType : size_t {
	kFluid,
	kElastic,
	kFluid3D,
	kElastic3D,
};


class KConfigFactory {
public:
	KConfigFactory() {
		 mContainer[kFluid] = AutoCreateInstance<KFluidConfig>();
		 mContainer[kElastic] = AutoCreateInstance<KElasticConfig>();
	};

	KConfig *CreateKConfig(KConfigType type) {
		return mContainer[type]->clone();
	};

private:
	using ContainerType = std::unordered_map<KConfigType, 
		                  std::shared_ptr<KConfig>, std::hash<size_t>>;

	ContainerType mContainer;
};


#endif /* ifndef __CONFIG_INCLUDE_H__ */





