# asdf2fits

Load [ASDF](https://www.asdf-format.org/) files as [DeepSkyFits](https://gitlab.in2p3.fr/deepskyproject/deepskyfits)
objects, and save those objects back to ASDF. Built on the
[libasdf](https://github.com/asdf-format/libasdf) C library.

```c++
#include <asdf2fits/ASDFmanager.h>
#include <asdf2fits/fits2asdf.h>
using namespace DSL;

ASDFmanager file("observation.asdf");

// Metadata only: no pixel or table data is read
std::shared_ptr<FITShdu> header = file.GetHeader("SCI");
double exptime = header->GetDoubleValueForKey("EXPTIME");

// Data, by name or by index (the first HDU is 1, as in FITSmanager)
std::shared_ptr<FITScube>  image = file.GetImage("SCI");
std::shared_ptr<FITStable> table = file.GetTable("CATALOG");
auto typed = std::dynamic_pointer_cast<FITSimg<float>>(image);  // FITSimg<T>, like any FITScube

image->Write("sci.fits");                                       // the DeepSkyFits API applies as usual

// Saving
fits2asdf(*image, "image.asdf");                                // as the primary HDU
fits2asdf(*table, "table.asdf", /*replace=*/true);              // after an empty primary HDU

FITSmanager fits("observation.fits");
fits2asdf(fits, "observation.asdf");                            // every HDU of a FITS file
```

Failures throw `FITSexception`.

## What is read

| In the ASDF file | Becomes |
|---|---|
| Value tagged `.../fits/fits-*` (astropy FITS layout) | One HDU per list entry, with its exact header |
| ndarray of a numeric type | `FITSimg<T>` (`float16` as `float`, `bool8` as `uint8_t`) |
| Structured ndarray, or value tagged `.../table-*` | `FITStable` |
| Scalars next to an ndarray | Keywords of its HDU, nested keys as `META.INSTRUMENT.NAME` |
| Scalars at the root of the tree | Keywords of the primary HDU |

- Without the FITS layout, an HDU is named after its tree path. The one called `PRIMARY`, or else the first found,
  is the primary HDU.
- Axes are reversed to follow FITS: an array of shape `[ny, nx]` has `NAXIS1 = nx`.
- Values are physical, so `BSCALE`/`BZERO` found in a header are not applied.
- NaN (real images) and the `BLANK` value (integer images) mark masked pixels.
- The FITS WCS is loaded from the header keywords.

Not converted:

- GWCS objects (Roman, JWST): such files load without a WCS.
- Table columns that aren't plain ndarrays (`Time`, `Quantity`, `SkyCoord`) are left out.
- Complex and string images, and ndarray `mask` properties, which libasdf doesn't read.
- Strided tables, which throw.

## What is written

`fits2asdf` uses the layout of [asdf-astropy](https://github.com/astropy/asdf-astropy)
(`tag:astropy.org:astropy/fits/fits-1.2.0`): a list of HDUs with `header` as `[key, value, comment]` cards, and
`data` as an ndarray. Tables are structured ndarrays, which that schema allows. In Python:

```python
import asdf
with asdf.open("observation.asdf") as f:
    hdul = f["fits"]    # astropy.io.fits.HDUList
```

- Data are written as physical values, so `TSCALn`/`TZEROn` are not carried over.
- Masked pixels become NaN, or a `BLANK` value for integers (the lowest value of a signed type, the highest of an
  unsigned one).
- Variable-length columns (`P`/`Q`) have no ASDF form and throw.

Two limits of astropy's `BinTableHDU` affect tables opened through astropy; plain `asdf` and this library read them
correctly:

- `int8` columns are read as booleans.
- `TUNITn` is ignored, so column units don't show up in `hdul[i].columns`.

## Build

Requires CMake 3.21 or newer, a C++20 compiler, autoconf, automake, libtool, pkg-config, lz4, zlib and bzip2.
On macOS: `brew install cmake autoconf automake libtool pkgconf lz4`.
DeepSkyFits (`DSTfits` 1.5 or newer) must already be installed in `$DeepSkyROOT` (by default, the parent directory).

```sh
git clone --recurse-submodules https://gitlab.in2p3.fr/deepskyproject/asdf2fits.git
cd asdf2fits
cmake -S . -B build
cmake --build build -j
(cd build && ctest)
cmake --install build      # to $DeepSkyROOT unless CMAKE_INSTALL_PREFIX is set
```

libasdf and its dependencies libfyaml and libmd are built from the submodules. libasdf is built from a copy with [a fix for an endless loop](cmake/patches/libasdf-document-end-marker.patch) applied, until upstream has it; the submodule itself stays unmodified. libfyaml and libmd are linked into
`libasdf`, so `libasdf2fits`, `libasdf` and their headers are all that gets installed. `-DWITH_ASAN=ON` builds with
AddressSanitizer.

To use the library from CMake:

```cmake
find_package(DSTfits 1.5 REQUIRED)
find_package(asdf2fits REQUIRED)
target_link_libraries(myapp asdf2fits::asdf2fits DSTfits::DSTfits)
```

## License

PolyForm Noncommercial 1.0.0, as DeepSkyFits. libasdf, libfyaml and libmd keep their own licences.
