//
//  ASDFhandle.h
//  asdf2fits
//
//  Private: RAII wrappers around libasdf handles and the column layout shared by the table code.
//
//  Licensed under PolyForm Noncommercial 1.0.0
//  You may use, share and adapt this code for noncommercial purposes.
//  Commercial use requires a separate licence from the author.
//  Licence text: https://polyformproject.org/licenses/noncommercial/1.0.0/

#ifndef _DSL_ASDFhandle_
#define _DSL_ASDFhandle_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <asdf.h>

#include <DSTfits/FITSexception.h>

namespace DSL::asdfdetail
{
    struct FileCloser    { void operator()(asdf_file_t* f)    const { asdf_close(f); } };
    struct ValueDeleter  { void operator()(asdf_value_t* v)   const { asdf_value_destroy(v); } };
    struct NdarrayDeleter{ void operator()(asdf_ndarray_t* n) const { asdf_ndarray_destroy(n); } };

    using FilePtr    = std::unique_ptr<asdf_file_t, FileCloser>;
    using ValuePtr   = std::unique_ptr<asdf_value_t, ValueDeleter>;
    using NdarrayPtr = std::unique_ptr<asdf_ndarray_t, NdarrayDeleter>;

    /** @brief Text of the last libasdf error on a file, never null. */
    inline std::string lastError(asdf_file_t* file)
    {
        const char* text = asdf_error(file);
        return text ? std::string(text) : std::string("unknown libasdf error");
    }

    /** @brief True if the value carries a tag containing the given text. */
    inline bool hasTag(asdf_value_t* value, const char* text)
    {
        const char* tag = asdf_value_tag(value);
        return tag && std::string(tag).find(text) != std::string::npos;
    }

    /** @brief Layout of one table column, inside a row or in its own ndarray. */
    struct ColumnSpec
    {
        std::string name;
        std::string unit;
        asdf_scalar_datatype_t type = ASDF_DATATYPE_UNKNOWN;
        bool     swap      = false;     //!< Stored in the other byte order than this machine
        uint64_t elemBytes = 0;         //!< One element (for a string: all its characters)
        uint64_t repeat    = 1;         //!< Elements per row
        uint64_t offset    = 0;         //!< Start of the column inside a structured row
        asdf_ndarray_t* source = nullptr;   //!< ndarray holding the bytes; owned by the caller
    };
}

#endif
