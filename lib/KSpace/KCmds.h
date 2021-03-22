



#ifndef __KCMDS_INCLUDE_H__
#define __KCMDS_INCLUDE_H__ 

#include <string>

struct KCmds {

	/// Output file name.
    std::string mOutputFileName;

    /// Store raw time-series of pressure over the sensor mask?
    bool mStorePressureRawFlag;
    /// Store RMS of pressure over the the sensor mask?
    bool mStorePressureRmsFlag;
    /// Store maximum of pressure over the sensor mask?
    bool mStorePressureMaxFlag;
    /// Store minimum of pressure over the sensor mask?
    bool mStorePressureMinFlag;
    /// Store maximum of pressure over the whole domain?
    bool mStorePressureMaxAllFlag;
    /// Store minimum of pressure over the whole domain?
    bool mStorePressureMinAllFlag;
    /// Store pressure in the final time step over the whole domain?
    bool mStorePressureFinalAllFlag;

	KCmds() {
		/// Output file name.
		mOutputFileName = "";
		/// Store raw time-series of pressure over the sensor mask?
		mStorePressureRawFlag = true;
		/// Store RMS of pressure over the the sensor mask?
		mStorePressureRmsFlag = false;
		/// Store maximum of pressure over the sensor mask?
		mStorePressureMaxFlag = false;
		/// Store minimum of pressure over the sensor mask?
		mStorePressureMinFlag = false;
		/// Store maximum of pressure over the whole domain?
		mStorePressureMaxAllFlag = false;
		/// Store minimum of pressure over the whole domain?
		mStorePressureMinAllFlag = false;
		/// Store pressure in the final time step over the whole domain?
		mStorePressureFinalAllFlag = false;
	}

};

#endif /* ifndef __KCMDS_INCLUDE_H__ */




