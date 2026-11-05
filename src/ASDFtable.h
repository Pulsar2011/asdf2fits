//
//  ASDFtable.h
//  asdf2fits
//
//  Private: conversion between ASDF tables and FITStable.
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFtable.h
 * @brief Conversion between ASDF tables and FITStable.
 * @note Private header: it isn't installed with the library.
 */
#ifndef _DSL_ASDFtable_
#define _DSL_ASDFtable_

#include <memory>
#include <vector>

#include <DSTfits/FITStable.h>

#include "ASDFhandle.h"

namespace DSL::asdfdetail
{
    /**
     * @brief Get the column layout of a structured ndarray, from metadata only.
     * @param ndarray Structured ndarray; the returned specs point to it.
     * @return One spec per field, in the order of the fields.
     * @throws FITSexception BAD_TFORM if a field is itself structured or has no size.
     */
    std::vector<ColumnSpec> structuredColumns(asdf_ndarray_t* ndarray);

    /**
     * @brief Get the column layout of a value tagged .../table-*, from metadata only.
     * @details A column is either a bare ndarray, named COLn, or a mapping with @c data and optional
     * @c name and @c unit. Dimensions after the first of an ndarray give the elements per row.
     * @param table  Value holding the table.
     * @param arrays Receives the column ndarrays, which the returned specs point to.
     * @return The columns. Those whose data isn't a plain, non-structured ndarray (Time, Quantity, ...) are left out.
     * @throws FITSexception NOT_TABLE if the value isn't a mapping with a list of columns.
     */
    std::vector<ColumnSpec> columnTable(asdf_value_t* table, std::vector<NdarrayPtr>& arrays);

    /**
     * @brief Load a table.
     * @details A column without a unit takes TUNITn from the header if TTYPEn there matches its name.
     * @param columns    Layout of the columns; their ndarrays must still be alive.
     * @param rows       Number of rows.
     * @param structured True if the columns share one structured ndarray, false if each has its own.
     * @param header     Header giving the keywords of the table.
     * @return The new table.
     * @throws FITSexception BAD_DIMEN if a source is strided, READ_ERROR if its data can't be read,
     * BAD_TFORM if a column has a type FITStable can't hold.
     */
    std::shared_ptr<FITStable> readTable(const std::vector<ColumnSpec>& columns, uint64_t rows, bool structured, const FITShdu& header);

    /**
     * @brief Copy a table into a new structured ndarray value of the file.
     * @details TSCALn, TZEROn, TNULLn and TDIMn are dropped from the header, and the other column
     * keywords are rebuilt from the columns. Strings are as wide as the longest value, padded with NUL.
     * @param file   File the value is created in.
     * @param table  Table to copy.
     * @param header Header of the output HDU; receives the keywords describing the table.
     * @return The new value.
     * @throws FITSexception BAD_DIMEN if a column's row count differs from the table's, BAD_TFORM for
     * variable-length or empty cells or an unknown data type, MEMORY_ALLOCATION if the row buffer can't
     * be allocated, WRITE_ERROR if the value can't be created.
     */
    asdf_value_t* writeTable(asdf_file_t* file, const FITStable& table, FITShdu& header);
}

#endif
