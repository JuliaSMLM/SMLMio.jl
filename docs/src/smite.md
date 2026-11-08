# SMITE files

SMITE is the lab's MATLAB package for single-molecule localization microscopy. Its results files
(`*_Results.mat`, MATLAB v7.3) hold the localizations in a structure called `SMD`. SMLMio reads
an SMD into a `SmiteSMLD` and writes any SMLD back as an SMD.

## Units

SMITE and SMLMCore place positions differently:

| | SMITE's SMD | SMLMCore |
|---|---|---|
| `X`, `Y` | pixels; (1, 1) is the center of the top-left pixel | microns from the camera's corner |
| `X_SE`, `Y_SE` | pixels | microns |
| `Z`, `Z_SE` | microns | microns |

Loading converts `x = (X - 0.5) * PixelSize` and `σ_x = X_SE * PixelSize` (the same for `y`);
saving converts back. The camera is an `IdealCamera` of `XSize` by `YSize` pixels of
`PixelSize` microns, so the center of SMITE's pixel (1, 1) is the center of the camera's first
pixel.

## Types

```julia
# Names an SMD in a .mat file
mutable struct SmiteSMD
    filepath::String    # directory holding the .mat file
    filename::String    # name of the .mat file
    varname::String     # variable name in the file (default "SMD")
end
SmiteSMD(filepath::String, filename::String) = SmiteSMD(filepath, filename, "SMD")

# An SMLD loaded from SMITE; metadata keeps the SMD's other fields
struct SmiteSMLD{T,E<:AbstractEmitter} <: AbstractSMLD
    emitters::Vector{E}
    camera::AbstractCamera
    n_frames::Int
    n_datasets::Int
    metadata::Dict{String,Any}
end
```

## Loading

```julia
using SMLMio

smd = SmiteSMD("path/to/Results/ResultsStructs", "Cell_01_Label_01_Results.mat")
smld = load_smite_2d(smd)      # Emitter2DFit{Float64}
smld3 = load_smite_3d(smd)     # Emitter3DFit{Float64}, needs Z and Z_SE

smld.metadata["pixel_size"]    # microns per pixel, from the file's PixelSize
smld.metadata["data_size"]     # [YSize, XSize]
smld.metadata["DriftX"]        # every other SMD field, as SMITE wrote it (drift in pixels)
```

What the loaders do with fields SMITE leaves empty or writes in another shape:

- **No pixel size.** When the file has no `PixelSize`, or an empty one, pass it:
  `load_smite_2d(smd; pixel_size=0.0978)`. Without it the loader throws; it never assumes one.
  When the file has a pixel size, a different `pixel_size` is an error.
- **Optional fields.** `Photons_SE` and `Bg_SE` load as `NaN`, `DatasetNum` as 1 and `ConnectID`
  as 0 when the file leaves them empty (SMITE leaves `ConnectID` empty before frame connection)
  or holds them with one value per dataset (some older SMITE versions). Their raw values stay in
  the metadata, `metadata["smite_defaulted_fields"]` lists them, and a non-empty field of the
  wrong length warns.
- **Required fields.** `X`, `Y`, `Photons`, `Bg`, `X_SE`, `Y_SE` and `FrameNum` (and `Z`, `Z_SE`
  in 3D) must hold one value per localization; otherwise the loader throws and names the field.
- **No localizations.** SMITE leaves `XSize` and `YSize` empty in a result without
  localizations, so no camera can be built and the loader throws, naming the missing image size.
- **3D.** SMITE's SMD has no `ZSize`; when a file has one, it is appended to `data_size`.
- Rows or columns: an SMD built by hand in MATLAB with row vectors loads like SMITE's columns.
- Localizations with a non-zero imaginary part in any position, photon or background field are
  dropped with a warning.

## Saving

```julia
save_smite(smld, "output/directory", "processed_Results.mat")
```

`save_smite` takes any SMLD (a `SmiteSMLD`, a `BasicSMLD`, ...) and writes MATLAB v7.3 that
`load_smite_2d` and `load_smite_3d` read back to the same emitters, camera, frames and datasets:

- `X`, `Y`, `X_SE`, `Y_SE` in SMITE's pixels; `Z` and `Z_SE` in microns, for 3D data
- `XSize`, `YSize` and `PixelSize` from the camera, and `NDims`. SMITE has one pixel size, so a
  camera whose pixels are not uniform and square is an error.
- every metadata entry whose name is not an emitter column, except the loaders' own
  `data_size`, `pixel_size` and `smite_defaulted_fields`

Example workflow:

```julia
using SMLMio
using SMLMCore: @filter, filter_roi
using Dates

smld = load_smite_2d(SmiteSMD("raw_data", "experiment1_Results.mat"))

# Positions are microns: keep bright localizations in a 10 µm square
bright = @filter(smld, photons > 1000)
roi = filter_roi(bright, (5.0, 15.0), (5.0, 15.0))

roi.metadata["analysis_date"] = string(Dates.now())
save_smite(roi, "processed_data", "experiment1_processed_Results.mat")
```
