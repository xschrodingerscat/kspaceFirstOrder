//
// Created by shaune on 2021/3/24.
//

#include "KIndexOutputStream.h"

#include <algorithm>
#include <Parameters/Parameters.h>

/**
 * Constructor - links the HDF5 dataset, SourceMatrix, SensorMask and reduction operator together.
 */
KIndexOutputStream::KIndexOutputStream(Hdf5File&            file,
                                     const MatrixName&    datasetName,
                                     const RealMatrix&    sourceMatrix,
                                     const IndexMatrix&   sensorMask,
                                     const ReduceOperator reduceOp,
                                     float*               bufferToReuse)
        : BaseOutputStream(file, datasetName, sourceMatrix, reduceOp, bufferToReuse),
          mSensorMask(sensorMask) { }

/**
 * Destructor.
 */
KIndexOutputStream::~KIndexOutputStream() { }// end of ~IndexOutputStream
//----------------------------------------------------------------------------------------------------------------------

/**
 * Create a HDF5 stream, create a dataset, and allocate data for it.
 */
void KIndexOutputStream::create()
{
    auto & params = Parameters::getInstance();
    auto & koutput = params.getKOutput();

    KMatrixType<float> elem;
    for (size_t i = 0; i < mBufferSize; ++ i)
    {
        elem.push_back(std::vector<float>());
    }

    koutput.mPressureRaw = KMatrix<float>(elem);
}// end of create
//----------------------------------------------------------------------------------------------------------------------

/**
 * Reopen the output stream after restart and reload data.
 */
void KIndexOutputStream::reopen() { }// end of reopen
//----------------------------------------------------------------------------------------------------------------------

/**
 * Sample grid points, line them up in the buffer an flush to the disk unless a reduction operator is applied.
 */
void KIndexOutputStream::sample()
{
    const float*  sourceData = mSourceMatrix.getData();
    const size_t* sensorData = mSensorMask.getData();

    auto & params = Parameters::getInstance();

    switch (mReduceOp)
    {
        case ReduceOperator::kNone:
        {
            auto & koutput = params.getKOutput();
            for (size_t i = 0; i < mBufferSize; i++)
            {
                auto data = sourceData[sensorData[i]];
                koutput.mPressureRaw[i].push_back(data);
            }
            break;
        }// case kNone

        default:
            break;
    }// switch
}// end of sample
//----------------------------------------------------------------------------------------------------------------------

/**
 * Apply post-processing on the buffer and flush it to the file.
 */
void KIndexOutputStream::postProcess() { }// end of postProcess
//----------------------------------------------------------------------------------------------------------------------

/**
 * Checkpoint the stream.
 */
void KIndexOutputStream::checkpoint() { }
//----------------------------------------------------------------------------------------------------------------------

/**
 * Close stream (apply post-processing if necessary, flush data and close).
 */
void KIndexOutputStream::close() { }// end of close
//----------------------------------------------------------------------------------------------------------------------
