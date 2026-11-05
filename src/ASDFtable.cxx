//
//  ASDFtable.cxx
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

#include "ASDFtable.h"

#include <algorithm>
#include <any>
#include <bit>
#include <cstring>
#include <limits>
#include <regex>
#include <typeindex>
#include <utility>

#include "ASDFheader.h"

namespace DSL::asdfdetail
{
    namespace
    {
        const char* const kWho = "asdf2fits";

        const asdf_byteorder_t kHostOrder = (std::endian::native == std::endian::little) ? ASDF_BYTEORDER_LITTLE : ASDF_BYTEORDER_BIG;

        void reverseBytes(void* data, size_t size)
        {
            auto* bytes = static_cast<unsigned char*>(data);
            std::reverse(bytes, bytes + size);
        }

        // One stored element converted to its in-memory form
        template<typename T> struct Element
        {
            static constexpr size_t bytes = sizeof(T);
            static T read(const uint8_t* p, bool swap)
            {
                T value;
                std::memcpy(&value, p, sizeof(T));
                if(swap)
                    reverseBytes(&value, sizeof(T));
                return value;
            }
        };

        template<> struct Element<bool>
        {
            static constexpr size_t bytes = 1;
            static bool read(const uint8_t* p, bool) { return *p != 0; }
        };

        template<typename R> struct Element<std::pair<R, R>>
        {
            static constexpr size_t bytes = 2 * sizeof(R);
            static std::pair<R, R> read(const uint8_t* p, bool swap)
            {
                return { Element<R>::read(p, swap), Element<R>::read(p + sizeof(R), swap) };
            }
        };

        void appendUtf8(std::string& out, uint32_t c)
        {
            if(c < 0x80)
                out += static_cast<char>(c);
            else if(c < 0x800)
            {
                out += static_cast<char>(0xC0 | (c >> 6));
                out += static_cast<char>(0x80 | (c & 0x3F));
            }
            else if(c < 0x10000)
            {
                out += static_cast<char>(0xE0 | (c >> 12));
                out += static_cast<char>(0x80 | ((c >> 6) & 0x3F));
                out += static_cast<char>(0x80 | (c & 0x3F));
            }
            else
            {
                out += static_cast<char>(0xF0 | (c >> 18));
                out += static_cast<char>(0x80 | ((c >> 12) & 0x3F));
                out += static_cast<char>(0x80 | ((c >> 6) & 0x3F));
                out += static_cast<char>(0x80 | (c & 0x3F));
            }
        }

        std::string readString(const ColumnSpec& spec, const uint8_t* p)
        {
            std::string text;
            if(spec.type == ASDF_DATATYPE_ASCII)
            {
                text.assign(reinterpret_cast<const char*>(p), spec.elemBytes);
            }
            else
            {
                for(uint64_t i = 0; i + 4 <= spec.elemBytes; i += 4)
                    appendUtf8(text, Element<uint32_t>::read(p + i, spec.swap));
            }

            // Strings are padded with NUL up to the width of the field
            const size_t last = text.find_last_not_of('\0');
            text.erase(last == std::string::npos ? 0 : last + 1);
            return text;
        }

        // Fills a column, skipping the per-row metadata update: only the first value goes through it
        template<typename T>
        std::shared_ptr<FITSform> finish(std::shared_ptr<FITScolumn<T>> column, std::vector<T>&& values)
        {
            if(!values.empty())
            {
                T first = values[0];
                column->push_back(std::any(first));
                column->template values<T>() = std::move(values);
            }
            return column;
        }

        template<typename T>
        std::shared_ptr<FITSform> readNumeric(const ColumnSpec& spec, const uint8_t* base, uint64_t rowBytes, uint64_t rows, dtype type, size_t position)
        {
            if(spec.repeat == 1)
            {
                std::vector<T> values(rows);
                for(uint64_t r = 0; r < rows; r++)
                    values[r] = Element<T>::read(base + r * rowBytes, spec.swap);
                auto column = std::make_shared<FITScolumn<T>>(spec.name, type, spec.unit, position);
                return finish<T>(column, std::move(values));
            }

            std::vector<std::vector<T>> values(rows, std::vector<T>(spec.repeat));
            for(uint64_t r = 0; r < rows; r++)
                for(uint64_t i = 0; i < spec.repeat; i++)
                    values[r][i] = Element<T>::read(base + r * rowBytes + i * Element<T>::bytes, spec.swap);

            auto column = std::make_shared<FITScolumn<std::vector<T>>>(spec.name, type, static_cast<int64_t>(spec.repeat), static_cast<int64_t>(Element<T>::bytes), spec.unit, position);
            return finish<std::vector<T>>(column, std::move(values));
        }

        std::shared_ptr<FITSform> readStrings(const ColumnSpec& spec, const uint8_t* base, uint64_t rowBytes, uint64_t rows, size_t position)
        {
            const int64_t width = static_cast<int64_t>(spec.type == ASDF_DATATYPE_UCS4 ? spec.elemBytes / 4 : spec.elemBytes);
            if(spec.repeat == 1)
            {
                std::vector<std::string> values(rows);
                for(uint64_t r = 0; r < rows; r++)
                    values[r] = readString(spec, base + r * rowBytes);
                auto column = std::make_shared<FITScolumn<std::string>>(spec.name, tstring, int64_t(1), width, spec.unit, position);
                return finish<std::string>(column, std::move(values));
            }

            std::vector<std::vector<std::string>> values(rows, std::vector<std::string>(spec.repeat));
            for(uint64_t r = 0; r < rows; r++)
                for(uint64_t i = 0; i < spec.repeat; i++)
                    values[r][i] = readString(spec, base + r * rowBytes + i * spec.elemBytes);
            auto column = std::make_shared<FITScolumn<std::vector<std::string>>>(spec.name, tstring, static_cast<int64_t>(spec.repeat), width, spec.unit, position);
            return finish<std::vector<std::string>>(column, std::move(values));
        }

        std::shared_ptr<FITSform> readColumn(const ColumnSpec& spec, const uint8_t* base, uint64_t rowBytes, uint64_t rows, size_t position)
        {
            switch(spec.type)
            {
                case ASDF_DATATYPE_INT8:       return readNumeric<int8_t>  (spec, base, rowBytes, rows, tsbyte, position);
                case ASDF_DATATYPE_UINT8:      return readNumeric<uint8_t> (spec, base, rowBytes, rows, tbyte, position);
                case ASDF_DATATYPE_INT16:      return readNumeric<int16_t> (spec, base, rowBytes, rows, tshort, position);
                case ASDF_DATATYPE_UINT16:     return readNumeric<uint16_t>(spec, base, rowBytes, rows, tushort, position);
                case ASDF_DATATYPE_INT32:      return readNumeric<int32_t> (spec, base, rowBytes, rows, tint, position);
                case ASDF_DATATYPE_UINT32:     return readNumeric<uint32_t>(spec, base, rowBytes, rows, tuint, position);
                case ASDF_DATATYPE_INT64:      return readNumeric<int64_t> (spec, base, rowBytes, rows, tlonglong, position);
                case ASDF_DATATYPE_UINT64:     return readNumeric<uint64_t>(spec, base, rowBytes, rows, tulonglong, position);
                case ASDF_DATATYPE_FLOAT32:    return readNumeric<float>   (spec, base, rowBytes, rows, tfloat, position);
                case ASDF_DATATYPE_FLOAT64:    return readNumeric<double>  (spec, base, rowBytes, rows, tdouble, position);
                case ASDF_DATATYPE_COMPLEX64:  return readNumeric<FITSform::complex>(spec, base, rowBytes, rows, tcplx, position);
                case ASDF_DATATYPE_COMPLEX128: return readNumeric<FITSform::dblcomplex>(spec, base, rowBytes, rows, tdbcplx, position);
                case ASDF_DATATYPE_BOOL8:      return readNumeric<bool>    (spec, base, rowBytes, rows, tlogical, position);
                case ASDF_DATATYPE_ASCII:
                case ASDF_DATATYPE_UCS4:       return readStrings(spec, base, rowBytes, rows, position);
                default:
                    throw FITSexception(BAD_TFORM, kWho, "readTable", "column '"+spec.name+"' has a datatype that FITStable can't hold");
            }
        }

        void requireContiguous(const asdf_ndarray_t* ndarray, uint64_t elementBytes, const std::string& what)
        {
            if(!ndarray->strides)
                return;

            int64_t expected = static_cast<int64_t>(elementBytes);
            for(uint32_t i = ndarray->ndim; i > 0; i--)
            {
                if(ndarray->strides[i - 1] != expected)
                    throw FITSexception(BAD_DIMEN, kWho, "readTable", what+" is strided, which isn't supported for tables");
                expected *= static_cast<int64_t>(ndarray->shape[i - 1]);
            }
        }

        // Size in bytes of one element of a scalar or string datatype of a field or column
        uint64_t elementBytes(const asdf_datatype_t& datatype, uint64_t count)
        {
            if(datatype.type == ASDF_DATATYPE_ASCII || datatype.type == ASDF_DATATYPE_UCS4)
                return count > 0 ? datatype.size / count : datatype.size;
            return asdf_scalar_datatype_size(datatype.type);
        }

        // Writing --------------------------------------------------------------------------------

        struct OutColumn
        {
            const FITSform* form = nullptr;
            ColumnSpec spec;
            std::vector<uint64_t> fieldShape;
        };

        template<typename T> void writeElement(uint8_t* p, const T& value) { std::memcpy(p, &value, sizeof(T)); }
        inline void writeElement(uint8_t* p, bool value) { *p = value ? 1 : 0; }
        template<typename R> void writeElement(uint8_t* p, const std::pair<R, R>& value)
        {
            std::memcpy(p, &value.first, sizeof(R));
            std::memcpy(p + sizeof(R), &value.second, sizeof(R));
        }
        inline void writeString(uint8_t* p, const std::string& value, uint64_t width)
        {
            std::memset(p, 0, width);
            std::memcpy(p, value.data(), std::min<uint64_t>(value.size(), width));
        }

        template<typename T> asdf_scalar_datatype_t outType()
        {
            if constexpr (std::is_same_v<T, bool>)                  return ASDF_DATATYPE_BOOL8;
            if constexpr (std::is_same_v<T, int8_t>)                return ASDF_DATATYPE_INT8;
            if constexpr (std::is_same_v<T, uint8_t>)               return ASDF_DATATYPE_UINT8;
            if constexpr (std::is_same_v<T, int16_t>)               return ASDF_DATATYPE_INT16;
            if constexpr (std::is_same_v<T, uint16_t>)              return ASDF_DATATYPE_UINT16;
            if constexpr (std::is_same_v<T, int32_t>)               return ASDF_DATATYPE_INT32;
            if constexpr (std::is_same_v<T, uint32_t>)              return ASDF_DATATYPE_UINT32;
            if constexpr (std::is_same_v<T, int64_t>)               return ASDF_DATATYPE_INT64;
            if constexpr (std::is_same_v<T, uint64_t>)              return ASDF_DATATYPE_UINT64;
            if constexpr (std::is_same_v<T, float>)                 return ASDF_DATATYPE_FLOAT32;
            if constexpr (std::is_same_v<T, double>)                return ASDF_DATATYPE_FLOAT64;
            if constexpr (std::is_same_v<T, FITSform::complex>)     return ASDF_DATATYPE_COMPLEX64;
            if constexpr (std::is_same_v<T, FITSform::dblcomplex>)  return ASDF_DATATYPE_COMPLEX128;
        }

        // Each method recognises one kind of payload and fills in the column layout for it
        struct LayoutVisitor
        {
            OutColumn& column;

            template<typename T> bool scalar()
            {
                if(column.form->payloadType() != std::type_index(typeid(T)))
                    return false;
                column.spec.type = outType<T>();
                column.spec.elemBytes = Element<T>::bytes;
                column.spec.repeat = 1;
                return true;
            }

            template<typename T> bool vector()
            {
                if(column.form->payloadType() != std::type_index(typeid(std::vector<T>)))
                    return false;
                const auto& values = column.form->template values<std::vector<T>>();
                uint64_t repeat = 0;
                checkCells(values, repeat);
                column.spec.type = outType<T>();
                column.spec.elemBytes = Element<T>::bytes;
                column.spec.repeat = repeat;
                return true;
            }

            bool strings()
            {
                if(column.form->payloadType() == std::type_index(typeid(std::string)))
                {
                    const auto& values = column.form->template values<std::string>();
                    uint64_t width = 1;
                    for(const std::string& v : values)
                        width = std::max<uint64_t>(width, v.size());
                    column.spec.type = ASDF_DATATYPE_ASCII;
                    column.spec.elemBytes = width;
                    column.spec.repeat = 1;
                    return true;
                }
                if(column.form->payloadType() == std::type_index(typeid(std::vector<std::string>)))
                {
                    const auto& values = column.form->template values<std::vector<std::string>>();
                    uint64_t repeat = 0, width = 1;
                    checkCells(values, repeat);
                    for(const auto& cell : values)
                        for(const std::string& v : cell)
                            width = std::max<uint64_t>(width, v.size());
                    column.spec.type = ASDF_DATATYPE_ASCII;
                    column.spec.elemBytes = width;
                    column.spec.repeat = repeat;
                    return true;
                }
                return false;
            }

            template<typename V> void checkCells(const std::vector<V>& values, uint64_t& repeat) const
            {
                if(values.empty())
                {
                    repeat = std::max<int64_t>(column.form->getNelem(), 1);
                    return;
                }
                repeat = values.front().size();
                for(const V& cell : values)
                    if(cell.size() != repeat)
                        throw FITSexception(BAD_TFORM, kWho, "writeTable", "column '"+column.form->getName()+"' has variable-length cells, which have no ASDF form");
                if(repeat == 0)
                    throw FITSexception(BAD_TFORM, kWho, "writeTable", "column '"+column.form->getName()+"' has empty cells, which have no ASDF form");
            }
        };

        // Without any row the storage type isn't known: rely on the declared FITS type
        void layoutFromType(OutColumn& column)
        {
            ColumnSpec& spec = column.spec;
            spec.repeat = std::max<int64_t>(column.form->getNelem(), 1);
            switch(column.form->getType())
            {
                case tsbyte:     spec.type = ASDF_DATATYPE_INT8;       break;
                case tbyte:      spec.type = ASDF_DATATYPE_UINT8;      break;
                case tshort:     spec.type = ASDF_DATATYPE_INT16;      break;
                case tushort:    spec.type = ASDF_DATATYPE_UINT16;     break;
                case tint:
                case tlong:      spec.type = ASDF_DATATYPE_INT32;      break;
                case tuint:
                case tulong:     spec.type = ASDF_DATATYPE_UINT32;     break;
                case tlonglong:  spec.type = ASDF_DATATYPE_INT64;      break;
                case tulonglong: spec.type = ASDF_DATATYPE_UINT64;     break;
                case tfloat:     spec.type = ASDF_DATATYPE_FLOAT32;    break;
                case tdouble:    spec.type = ASDF_DATATYPE_FLOAT64;    break;
                case tcplx:      spec.type = ASDF_DATATYPE_COMPLEX64;  break;
                case tdbcplx:    spec.type = ASDF_DATATYPE_COMPLEX128; break;
                case tlogical:
                case tbit:       spec.type = ASDF_DATATYPE_BOOL8;      break;
                default:         spec.type = ASDF_DATATYPE_ASCII;      break;
            }
            spec.elemBytes = (spec.type == ASDF_DATATYPE_ASCII)
                ? static_cast<uint64_t>(std::max<int64_t>(column.form->getWidth(), 1))
                : asdf_scalar_datatype_size(spec.type);
        }

        template<typename F> bool forEachElementType(F&& f)
        {
            return f.template operator()<bool>()    || f.template operator()<int8_t>()   || f.template operator()<uint8_t>()
                || f.template operator()<int16_t>() || f.template operator()<uint16_t>() || f.template operator()<int32_t>()
                || f.template operator()<uint32_t>()|| f.template operator()<int64_t>()  || f.template operator()<uint64_t>()
                || f.template operator()<float>()   || f.template operator()<double>()
                || f.template operator()<FITSform::complex>() || f.template operator()<FITSform::dblcomplex>();
        }

        void computeLayout(OutColumn& column, uint64_t rows)
        {
            LayoutVisitor visitor{column};
            const bool known = visitor.strings() || forEachElementType([&]<typename T>() {
                return visitor.template scalar<T>() || visitor.template vector<T>();
            });
            if(!known)
            {
                if(rows > 0)
                    throw FITSexception(BAD_TFORM, kWho, "writeTable", "column '"+column.form->getName()+"' holds data of an unknown type");
                layoutFromType(column);
            }
        }

        void copyColumn(const OutColumn& column, uint8_t* base, uint64_t rowBytes, uint64_t rows)
        {
            const FITSform& form = *column.form;
            const ColumnSpec& spec = column.spec;
            if(rows == 0)
                return;

            const bool copied = forEachElementType([&]<typename T>() {
                if(form.payloadType() == std::type_index(typeid(T)))
                {
                    const auto& values = form.template values<T>();
                    for(uint64_t r = 0; r < rows; r++)
                    {
                        const T value = values[r];
                        writeElement(base + r * rowBytes, value);
                    }
                    return true;
                }
                if(form.payloadType() == std::type_index(typeid(std::vector<T>)))
                {
                    const auto& values = form.template values<std::vector<T>>();
                    for(uint64_t r = 0; r < rows; r++)
                        for(uint64_t i = 0; i < spec.repeat; i++)
                        {
                            const T value = values[r][i];
                            writeElement(base + r * rowBytes + i * spec.elemBytes, value);
                        }
                    return true;
                }
                return false;
            });
            if(copied)
                return;

            if(form.payloadType() == std::type_index(typeid(std::string)))
            {
                const auto& values = form.template values<std::string>();
                for(uint64_t r = 0; r < rows; r++)
                    writeString(base + r * rowBytes, values[r], spec.elemBytes);
            }
            else
            {
                const auto& values = form.template values<std::vector<std::string>>();
                for(uint64_t r = 0; r < rows; r++)
                    for(uint64_t i = 0; i < spec.repeat; i++)
                        writeString(base + r * rowBytes + i * spec.elemBytes, values[r][i], spec.elemBytes);
            }
        }
    }

    std::vector<ColumnSpec> structuredColumns(asdf_ndarray_t* ndarray)
    {
        std::vector<ColumnSpec> columns;
        uint64_t offset = 0;
        for(uint32_t i = 0; i < ndarray->datatype.nfields; i++)
        {
            const asdf_datatype_t& field = ndarray->datatype.fields[i];

            uint64_t repeat = 1;
            for(uint32_t d = 0; d < field.ndim; d++)
                repeat *= field.shape[d];

            ColumnSpec spec;
            spec.name = field.name ? field.name : "COL" + std::to_string(i + 1);
            spec.type = field.type;
            spec.repeat = repeat;
            spec.source = ndarray;
            spec.offset = offset;
            spec.swap = (field.byteorder != ASDF_BYTEORDER_DEFAULT ? field.byteorder : ndarray->byteorder) != kHostOrder
                     && (field.byteorder != ASDF_BYTEORDER_DEFAULT || ndarray->byteorder != ASDF_BYTEORDER_DEFAULT);
            spec.elemBytes = elementBytes(field, repeat);

            if(field.type == ASDF_DATATYPE_STRUCTURED || spec.elemBytes == 0)
                throw FITSexception(BAD_TFORM, kWho, "structuredColumns", "field '"+spec.name+"' has a datatype that FITStable can't hold");

            offset += spec.elemBytes * repeat;
            columns.push_back(spec);
        }
        return columns;
    }

    std::vector<ColumnSpec> columnTable(asdf_value_t* table, std::vector<NdarrayPtr>& arrays)
    {
        std::vector<ColumnSpec> columns;

        asdf_mapping_t* mapping = nullptr;
        if(asdf_value_as_mapping(table, &mapping) != ASDF_VALUE_OK)
            throw FITSexception(NOT_TABLE, kWho, "columnTable", "the table isn't a mapping");

        ValuePtr list(asdf_mapping_get(mapping, "columns"));
        asdf_sequence_t* sequence = nullptr;
        if(!list || asdf_value_as_sequence(list.get(), &sequence) != ASDF_VALUE_OK)
            throw FITSexception(NOT_TABLE, kWho, "columnTable", "the table has no list of columns");

        asdf_sequence_iter_t* it = asdf_sequence_iter_init(sequence);
        while(asdf_sequence_iter_next(&it))
        {
            ValuePtr dataHolder;
            asdf_value_t* data = it->value;
            std::string name = "COL" + std::to_string(it->index + 1);
            std::string unit;

            asdf_mapping_t* column = nullptr;
            if(!hasTag(it->value, "core/ndarray-") && asdf_value_as_mapping(it->value, &column) == ASDF_VALUE_OK)
            {
                ValuePtr nameValue(asdf_mapping_get(column, "name"));
                ValuePtr unitValue(asdf_mapping_get(column, "unit"));
                const char* text = nullptr;
                if(nameValue && asdf_value_as_scalar0(nameValue.get(), &text) == ASDF_VALUE_OK && text)
                    name = text;
                if(unitValue && asdf_value_as_scalar0(unitValue.get(), &text) == ASDF_VALUE_OK && text)
                    unit = text;

                dataHolder.reset(asdf_mapping_get(column, "data"));
                data = dataHolder.get();
            }

            asdf_ndarray_t* raw = nullptr;
            if(!data || !hasTag(data, "core/ndarray-") || asdf_value_as_ndarray(data, &raw) != ASDF_VALUE_OK)
                continue;   // Time, Quantity, SkyCoord... have no FITStable column

            NdarrayPtr array(raw);
            if(array->ndim < 1 || array->datatype.type == ASDF_DATATYPE_STRUCTURED)
                continue;

            uint64_t repeat = 1;
            for(uint32_t d = 1; d < array->ndim; d++)
                repeat *= array->shape[d];

            ColumnSpec spec;
            spec.name = name;
            spec.unit = unit;
            spec.type = array->datatype.type;
            spec.repeat = repeat;
            spec.source = array.get();
            spec.swap = array->byteorder != ASDF_BYTEORDER_DEFAULT && array->byteorder != kHostOrder;
            spec.elemBytes = elementBytes(array->datatype, 1);
            if(spec.elemBytes == 0)
                continue;

            arrays.push_back(std::move(array));
            columns.push_back(spec);
        }
        return columns;
    }

    std::shared_ptr<FITStable> readTable(const std::vector<ColumnSpec>& columns, uint64_t rows, bool structured, const FITShdu& header)
    {
        auto table = std::make_shared<FITStable>();

        uint64_t structuredRow = 0;
        for(const ColumnSpec& spec : columns)
            structuredRow += spec.elemBytes * spec.repeat;

        size_t position = 1;
        for(ColumnSpec spec : columns)
        {
            const uint64_t rowBytes = structured ? structuredRow : spec.elemBytes * spec.repeat;
            requireContiguous(spec.source, rowBytes, "column '"+spec.name+"'");

            if(spec.unit.empty())
            {
                const std::string n = std::to_string(position);
                if(header.Exists("TUNIT" + n) && header.Exists("TTYPE" + n) && header.GetValueForKey("TTYPE" + n) == spec.name)
                    spec.unit = header.GetValueForKey("TUNIT" + n);
            }

            size_t size = 0;
            const uint8_t* base = nullptr;
            if(rows > 0)
            {
                base = static_cast<const uint8_t*>(asdf_ndarray_data(spec.source, &size));
                if(!base || size < rows * rowBytes)
                    throw FITSexception(READ_ERROR, kWho, "readTable", "can't read the data of column '"+spec.name+"'");
            }

            table->InsertColumn(readColumn(spec, base ? base + spec.offset : nullptr, rowBytes, rows, position));
            position++;
        }

        table->HDU().Import(header, {"TFIELDS"});
        for(const char* key : {"COMMENT", "HISTORY", "EXTNAME"})
            if(header.Exists(key))
                table->HDU().ValueForKey(key, header.GetEntry(key)->second.value(), fChar);
        return table;
    }

    asdf_value_t* writeTable(asdf_file_t* file, const FITStable& table, FITShdu& header)
    {
        const uint64_t rows = table.nrows();

        std::vector<OutColumn> columns(table.ncols());
        std::vector<ColumnSpec> specs;
        uint64_t rowBytes = 0;
        for(size_t i = 0; i < columns.size(); i++)
        {
            OutColumn& column = columns[i];
            column.form = table.getColumn(i + 1).get();
            column.spec.name = column.form->getName();
            column.spec.unit = column.form->getUnit();
            if(column.form->size() != rows)
                throw FITSexception(BAD_DIMEN, kWho, "writeTable", "column '"+column.spec.name+"' doesn't have as many rows as the table");
            computeLayout(column, rows);
            column.spec.offset = rowBytes;
            rowBytes += column.spec.elemBytes * column.spec.repeat;
            specs.push_back(column.spec);
        }

        // Scaling and null values don't apply to the physical values written here, and the column
        // keywords are rebuilt from the columns themselves.
        static const std::regex stale("(TSCAL|TZERO|TNULL|TTYPE|TFORM|TUNIT|TDIM)[0-9]+");
        std::vector<std::string> remove;
        for(auto it = header.begin(); it != header.end(); ++it)
            if(std::regex_match(it->first, stale))
                remove.push_back(it->first);
        for(const std::string& key : remove)
            header.DeleteKey(key);
        setTableStructure(header, specs, rows);

        // Fields point into these until the datatype is serialized
        std::vector<asdf_datatype_t> fields(columns.size());
        for(size_t i = 0; i < columns.size(); i++)
        {
            OutColumn& column = columns[i];
            asdf_datatype_t& field = fields[i];
            field = asdf_datatype_t{};
            field.type = column.spec.type;
            field.name = column.spec.name.c_str();
            field.byteorder = kHostOrder;
            // libasdf writes this as the length of a string field even when it also has a shape
            field.size = column.spec.elemBytes;
            if(column.spec.repeat > 1)
            {
                column.fieldShape = {column.spec.repeat};
                field.ndim = 1;
                field.shape = column.fieldShape.data();
            }
        }

        const uint64_t shape[1] = {rows};
        asdf_ndarray_t ndarray{};
        ndarray.ndim = 1;
        ndarray.shape = shape;
        ndarray.datatype.type = ASDF_DATATYPE_STRUCTURED;
        ndarray.datatype.nfields = static_cast<uint32_t>(fields.size());
        ndarray.datatype.fields = fields.data();
        ndarray.datatype.size = rowBytes;
        ndarray.byteorder = kHostOrder;

        auto* buffer = static_cast<uint8_t*>(asdf_ndarray_data_alloc(&ndarray));
        if(!buffer && rows > 0)
            throw FITSexception(MEMORY_ALLOCATION, kWho, "writeTable", lastError(file));

        for(const OutColumn& column : columns)
            copyColumn(column, buffer ? buffer + column.spec.offset : nullptr, rowBytes, rows);

        asdf_value_t* out = asdf_value_of_ndarray(file, &ndarray);
        if(!out)
        {
            asdf_ndarray_data_dealloc(&ndarray);
            throw FITSexception(WRITE_ERROR, kWho, "writeTable", lastError(file));
        }
        return out;
    }
}
