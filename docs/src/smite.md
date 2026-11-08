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
smld.metadata["PValue"]        # one value per row of the file; emitter e is row e.id
```

Each emitter's `id` is its row in the file. The fields SMITE stacks per localization (its
catSMD list: `PValue`, `LogLikelihood`, `PSFSigma`, `ThreshFlag`, `XBoxCorner`, ...) that the
file holds with one value per row are listed in `metadata["smite_per_localization_fields"]`;
index them by `id`, for example `smld.metadata["PValue"][e.id]`. `metadata["smite_row_frames"]`
and `metadata["smite_row_photons"]` record each row's frame and photons. Some older SMITE files hold such a field with one value
per dataset; in a file with as many localizations as datasets that cannot be told apart.

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
  `load_smite_2d` on a file with `NDims` 3 warns that it leaves `Z` out (kept in metadata).
- Rows or columns: an SMD built by hand in MATLAB with row vectors loads like SMITE's columns.
- Localizations with a non-zero imaginary part in any position, photon or background field are
  dropped with a warning.

## Saving

```julia
save_smite(smld, "output/directory", "processed_Results.mat")
```

`save_smite` takes any SMLD (a `SmiteSMLD`, a `BasicSMLD`, ...) and writes MATLAB v7.3 that
`load_smite_2d`, `load_smite_3d` and SMITE read. It writes:

- `X`, `Y`, `X_SE`, `Y_SE` in SMITE's pixels, and `Z`, `Z_SE` in microns for 3D data
- `XSize`, `YSize` and `PixelSize` from the camera, `NFrames`, `NDatasets` and `NDims`; counts
  and indices as doubles, as SMITE's results files hold them. SMITE has one pixel size, so a
  camera whose pixels are not uniform and square is an error.
- every metadata entry that is not an emitter column, except the loaders' own: `original_file`,
  `data_size`, `pixel_size`, `smite_defaulted_fields`, `smite_per_localization_fields`,
  `smite_row_frames`, `smite_row_photons`, `complex_fields_removed`, `complex_fields`,
  `original_emitter_count` and `removed_emitter_count`.
- SMITE's per-localization fields a loader kept, for the saved emitters' rows, so after a filter
  `PValue` still belongs to its localization. They are written only when every emitter is
  still the file row its `id` names (the same frame and photons as when loaded); for combined
  localizations (frame connection's combined output), emitters joined from another file, or
  changed frames or photons, they are dropped with a warning. A field replaced since loading
  that no longer has one value per file row is dropped with a warning too.
- An optional field the file left empty, and nobody has set since loading, back empty (an
  empty `ConnectID` stays empty).

Loading the saved file gives back the emitters' positions, uncertainties, photons, background,
frames, datasets and track ids, and `n_frames` and `n_datasets`. It does not give back:

- `σ_xy`: SMITE has no such field; it loads as 0.
- `id`: it loads as the row in the new file.
- the camera type: an sCMOS camera's calibration maps are not saved, and the file loads with an
  `IdealCamera` of the same pixels.
- a camera origin other than 0: SMITE has none, so `save_smite` warns, and the file loads with
  the camera at 0 and the positions moved by the same offset.

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
