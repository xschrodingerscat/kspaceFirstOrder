//
// Created by shaune on 2021/3/24.
//

#ifndef KSPACESOLVER_KINDEXOUTPUTSTREAM_H
#define KSPACESOLVER_KINDEXOUTPUTSTREAM_H


#include <OutputStreams/BaseOutputStream.h>

/**
 * @class   IndexOutputStream
 * @brief   Output stream for quantities sampled by an index sensor mask.
 * @details Output stream for quantities sampled by an index sensor mask. This class writes data to a single dataset in a
 *          root group of the HDF5 file (time-series as well as aggregations).
 */
class KIndexOutputStream : public BaseOutputStream
{
public:
    /// Default constructor not allowed.
    KIndexOutputStream() = delete;
    /**
     * @brief Constructor.
     *
     * @details The Constructor links the HDF5 dataset, source (sampled matrix), sensor mask and the reduction operator
     *          together. The constructor DOES NOT allocate memory because the size of the sensor mask is not known at
     *          the time the instance of the class is being created.
     *
     * @param [in] file          - HDF5 file to write the output to.
     * @param [in] datasetName   - HDF5 dataset name. Index based sensor data is stored in a single dataset.
     * @param [in] sourceMatrix  - Source real matrix to be sampled.
     * @param [in] sensorMask    - Index sensor mask.
     * @param [in] reduceOp      - Reduction operator.
     * @param [in] bufferToReuse - An external buffer can be used to line up the grid points.
     */
    KIndexOutputStream(Hdf5File&            file,
                      const MatrixName&    datasetName,
                      const RealMatrix&    sourceMatrix,
                      const IndexMatrix&   sensorMask,
                      const ReduceOperator reduceOp,
                      double*               bufferToReuse = nullptr);

    /// Copy constructor not allowed.
    KIndexOutputStream(const KIndexOutputStream&) = delete;

    /// Destructor.
    virtual ~KIndexOutputStream() override;

    /// Operator = not allowed.
    KIndexOutputStream& operator=(const KIndexOutputStream&) = delete;

    /// Create a HDF5 stream, allocate data for it and open the dataset.
    virtual void create()      override;

    /// Reopen the output stream after restart and reload data.
    virtual void reopen()      override;

    /// Sample grid points, line them up in the buffer, if necessary a reduce operator is applied.
    virtual void sample()      override;

    /// Apply post-processing on the buffer and flush it to the file.
    virtual void postProcess() override;

    /// Checkpoint the stream.
    virtual void checkpoint()  override;

    /// Close stream.
    virtual void close()       override;

protected:
    /// Sensor mask to sample data.
    const IndexMatrix& mSensorMask;
};// end of IndexOutputStream
//----------------------------------------------------------------------------------------------------------------------

#endif //KSPACESOLVER_KINDEXOUTPUTSTREAM_H


