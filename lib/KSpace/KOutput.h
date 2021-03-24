//
// Created by shaune on 2021/3/24.
//

#ifndef KSPACESOLVER_KOUTPUT_H
#define KSPACESOLVER_KOUTPUT_H


#include <string>
#include <Containers/OutputStreamContainer.h>

class KOutput
{
public:
    KOutput()
    {
        /// Output file name.
        mOutputFileName = "";
        /// Store raw time-series of pressure over the sensor mask?
        mStorePressureRawFlag = true;

        mOutputToFileFlag = false;
        /*
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
         */
    }

    KOutput(const KOutput& ) = default;

    ~KOutput() = default;

    bool isOutputToFileFlag() const
    {
        return mOutputToFileFlag;
    }

    const std::string &getOutputFileName() const {
        return mOutputFileName;
    }

    const KMatrix<float> &getPressureRaw() const
    {
        return mPressureRaw;
    }

    void setOutputFileName(const std::string &filename)
    {
        mOutputFileName = filename;
    }

    void setOutputToFileFlag(bool flag)
    {
        mOutputToFileFlag = flag;
    }

    /*
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
     */

    /// Output file name.
    std::string mOutputFileName;

    /// Store raw time-series of pressure over the sensor mask?
    bool mStorePressureRawFlag;

    bool mOutputToFileFlag;


    KMatrix<float> mPressureRaw;
};


#endif //KSPACESOLVER_KOUTPUT_H
