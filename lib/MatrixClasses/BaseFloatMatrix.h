/**
 * @file      BaseFloatMatrix.h
 *
 * @author    Jiri Jaros \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            jarosjir@fit.vutbr.cz
 *
 * @brief     The header file containing the base class for single precisions floating point numbers (floats).
 *
 * @version   kspaceFirstOrder 2.17
 *
 * @date      11 July      2011, 12:13 (created) \n
 *            11 February  2020, 14:43 (revised)
 *
 * @copyright Copyright (C) 2011 - 2020 SC\@FIT Research Group, Brno University of Technology, Brno, CZ.
 *
 * This file is part of the C++ extension of the [k-Wave Toolbox](http://www.k-wave.org).
 *
 * k-Wave is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with k-Wave.
 * If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).
 */

#ifndef BASE_FLOAT_MATRIX_H
#define BASE_FLOAT_MATRIX_H

#include <MatrixClasses/BaseMatrix.h>
#include <Utils/DimensionSizes.h>

/**
 * @class   BaseFloatMatrix
 * @brief   Abstract base class for double based matrices defining basic interface. Higher dimensional
 *          matrices stored as 1D arrays, row-major order.
 *
 * @details Abstract base class for double based matrices defining basic interface. Higher dimensional matrices stored
 *          as 1D arrays, row-major order. The I/O is done via HDF5 files.
 */
class BaseFloatMatrix : public BaseMatrix
{
  public:
    /// Default constructor.
    BaseFloatMatrix();
    /// Copy constructor not allowed.
    BaseFloatMatrix(const BaseFloatMatrix&) = delete;
    /// Destructor.
    virtual ~BaseFloatMatrix() override = default;

    /// Operator = not allowed.
    BaseFloatMatrix& operator=(const BaseFloatMatrix&) = delete;

    /**
     * @brief  Get dimension sizes of the matrix.
     * @return Dimension sizes of the matrix.
     */
    virtual const DimensionSizes& getDimensionSizes() const override { return mDimensionSizes; };
    /**
     * @brief  Size of the matrix.
     * @return Number of elements.
     */
    virtual size_t size()                             const override { return mSize; };
    /**
     * @brief  The capacity of the matrix (this may differ from size due to padding, etc.).
     * @return Capacity of the currently allocated storage.
     */
    virtual size_t capacity()                         const override { return mCapacity; };

    /**
     * @brief  Get matrix data out of the class (for direct kernel access).
     * @return Pointer to mutable matrix data.
     */
    virtual double*       getData()       { return mData; };
    /**
     * @brief Get matrix data out of the class (for direct kernel access), const version.
     * @return Pointer to immutable matrix data.
     */
    virtual const double* getData() const { return mData; };

    /**
     * @brief Copy data from other matrix with the same size.
     * @param [in] src - Matrix to copy data in.
     */
    virtual void copyData(const BaseFloatMatrix& src);

    /// Zero all elements of the matrix (NUMA first touch).
    virtual void zeroMatrix();
    /**
     * @brief Calculate matrix = scalar / matrix.
     * @param [in] scalar - Scalar constant
     */
    virtual void scalarDividedBy(const double scalar);

  protected:
   /**
    * @brief Aligned memory allocation with zeroing.
    * @throw std::bad_alloc - If there's not enough memory.
    */
    virtual void allocateMemory();
    /// Memory deallocation.
    virtual void freeMemory();

    /// Dimension sizes.
    DimensionSizes mDimensionSizes;

    /// Total number of used elements.
    size_t mSize;
    /// Total number of allocated elements (in terms of floats).
    size_t mCapacity;

    /// Raw matrix data.
    double* mData;

  private:

};//end of BaseFloatMatrix
//----------------------------------------------------------------------------------------------------------------------

#endif /* BASE_FLOAT_MATRIX_H */
