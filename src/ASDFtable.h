//
//  ASDFtable.h
//  asdf2fits
//
//  Private: conversion between ASDF tables and FITStable.
//
//  Licensed under PolyForm Noncommercial 1.0.0
//  You may use, share and adapt this code for noncommercial purposes.
//  Commercial use requires a separate licence from the author.
//  Licence text: https://polyformproject.org/licenses/noncommercial/1.0.0/

#ifndef _DSL_ASDFtable_
#define _DSL_ASDFtable_

#include <memory>
#include <vector>

#include <DSTfits/FITStable.h>

#include "ASDFhandle.h"

namespace DSL::asdfdetail
{
    /** @brief Column layout of a structured ndarray, read from metadata only. */
    std::vector<ColumnSpec> structuredColumns(asdf_ndarray_t* ndarray);

    /**
     * @brief Column layout of a value tagged .../table-*, read from metadata only.
     * @param arrays Receives the column ndarrays, which the returned specs point to.
     * @return Columns whose data isn't a plain ndarray (Time, Quantity, ...) are left out.
     */
    std::vector<ColumnSpec> columnTable(asdf_value_t* table, std::vector<NdarrayPtr>& arrays);

    /** @brief Load a table. A structured ndarray holds all columns in one source. */
    std::shared_ptr<FITStable> readTable(const std::vector<ColumnSpec>& columns, uint64_t rows, bool structured, const FITShdu& header);

    /**
     * @brief Copy a table into a new structured ndarray value of the file.
     * @param header Header of the output HDU; receives the keywords describing the table.
     */
    asdf_value_t* writeTable(asdf_file_t* file, const FITStable& table, FITShdu& header);
}

#endif
