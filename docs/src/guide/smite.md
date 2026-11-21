# SMITE files

SMITE is the lab's MATLAB package for single-molecule localization microscopy. Its results files
(`*_Results.mat`, MATLAB v7.3) hold the localizations in a structure called `SMD`. SMLMio reads
an SMD into a [`SmiteSMLD`](@ref) and writes SMLDs of fit emitters back as an SMD. This page
holds the rules for both directions; the docstrings of [`load_smite_2d`](@ref),
[`load_smite_3d`](@ref) and [`save_smite`](@ref) link here.

## Units

SMITE and SMLMCore place positions differently:

| | SMITE's SMD | SMLMCore |
|---|---|---|
| `X`, `Y` | pixels; (1, 1) is the center of the top-left pixel | µm from the camera's corner |
| `X_SE`, `Y_SE` | pixels | µm |
| `Z`, `Z_SE` | µm | µm |
| `Photons` | photons | photons |
| `Bg` | photons per pixel | photons per pixel |

Loading converts `x = (X - 0.5) * PixelSize` and `σ_x = X_SE * PixelSize` (the same for `y`);
saving converts back. The camera is an `IdealCamera` of `XSize` by `YSize` pixels of
`PixelSize` µm, so the center of SMITE's pixel (1, 1) is the center of the camera's first
pixel. `X` is the column and `Y` the row in both. Other SMD fields stay in SMITE's units in the
metadata (`PSFSigma` and `DriftX` are in pixels). SMLMCore's
[conventions](https://JuliaSMLM.github.io/SMLMCore.jl/stable/guide/conventions/) describe its
side.

## Loading

```julia
using SMLMio

smd = SmiteSMD("path/to/results", "localizations_Results.mat")
smld = load_smite_2d(smd)      # Emitter2DFit{Float64}
smld3 = load_smite_3d(smd)     # Emitter3DFit{Float64}, needs Z and Z_SE

smld.metadata["pixel_size"]    # µm per pixel, from the file's PixelSize
smld.metadata["data_size"]     # [YSize, XSize]
smld.metadata["DriftX"]        # every other SMD field, as SMITE wrote it (drift in pixels)
smld.metadata["PValue"]        # if listed (below): one value per file row; emitter e is row e.id
```

Each emitter's `id` is its row in the file. The fields SMITE stacks per localization (its
catSMD list: `PValue`, `LogLikelihood`, `PSFSigma`, `ThreshFlag`, `XBoxCorner`, ...) that the
file holds with one value per row are listed in `metadata["smite_per_localization_fields"]`;
index them by `id`, for example `smld.metadata["PValue"][e.id]`.
`metadata["smite_row_frames"]` and `metadata["smite_row_photons"]` record each row's frame and
photons. To edit or replace a listed field, keep one value per file row, in file-row order:
`save_smite` writes it by row. Some older SMITE files hold such a field with one value per
dataset; in a file with as many localizations as datasets that cannot be told apart.

`n_frames` is SMITE's `NFrames`, the number of frames per dataset, and in SMITE's convention
`frame` runs from 1 to `n_frames` in each dataset; the loaders do not check frames against it.
`n_datasets` is `NDatasets`. An `n_frames` of 0 means "unknown" in the lab's packages.

What the loaders do with fields SMITE leaves empty or writes in another shape:

- **No pixel size.** When the file has no `PixelSize`, or an empty one, pass it:
  `load_smite_2d(smd; pixel_size=0.0978)`. Without it the loader throws; it never assumes one.
  When the file has a pixel size, a different `pixel_size` is an error.
- **Pixel size type.** A MATLAB single `PixelSize` gives an `IdealCamera{Float32}`, as in
  SMLMData 0.7.1; any other pixel size becomes a `Float64` and gives an `IdealCamera{Float64}`.
  A `pixel_size` passed to the loader sets the type only when the file has no `PixelSize`
  (a `Float32` one gives a `Float32` camera). The emitters are `Float64` either way. A pixel
  size that is complex (with a non-zero imaginary part), not finite or not positive is an
  error, and so is a `Float32` one whose camera edges overflow or round onto each other.
- **Sizes and counts.** `XSize` and `YSize` must each be one whole number of at least 1, and
  `NFrames` and `NDatasets` one whole number of at least 0, in `Int`'s range and compared as
  stored, so a count above 2^53 is not rounded. A complex value with a non-zero imaginary part,
  text, several values or a value below the bound is an error naming the field (a complex value
  whose imaginary part is zero is its real part).
- **Optional fields.** `Photons_SE` and `Bg_SE` load as `NaN`, `DatasetNum` as 1 and `ConnectID`
  as 0 when the file leaves them empty (SMITE leaves `ConnectID` empty before frame connection)
  or holds them with one value per dataset (some older SMITE versions). Their raw values stay in
  the metadata, `metadata["smite_defaulted_fields"]` lists them, and a non-empty field of the
  wrong length warns.
- **Required fields.** `X`, `Y`, `Photons`, `Bg`, `X_SE`, `Y_SE` and `FrameNum` (and `Z`, `Z_SE`
  in 3D) must hold one value per localization; otherwise the loader throws and names the field.
  `FrameNum`, `DatasetNum` and `ConnectID` must hold whole numbers in `Int`'s range, or the
  loader throws and names the field.
- **No localizations.** SMITE leaves `XSize` and `YSize` empty in a result without
  localizations, so no camera can be built and the loader throws, naming the missing image size.
- **One value each.** SMITE's `catSMD` stores a list for a field that differs between the
  datasets it joins, such as `NFrames` when they differ in length. The loaders throw for a
  `PixelSize`, `XSize`, `YSize`, `NFrames` or `NDatasets` of several values, saying how many the
  file holds; a `pixel_size` passed to the loader does not stand in for such a `PixelSize`.
- **3D.** SMITE's SMD has no `ZSize`. Every `ZSize` stays in the metadata, so `save_smite` writes
  it back, and none is an error: when it is one whole number of at least 0 that fits an `Int`,
  `load_smite_3d` also appends it to `data_size`; any other `ZSize` (several values, text, -1)
  is only kept in the metadata.
  `load_smite_2d` on a file whose `Z` and `Z_SE` hold one value per localization warns that it leaves `Z`
  out (kept in the metadata); SMITE writes `NDims` 2 for its own 3D fits, so `NDims` does not
  say whether a file is 3D.
- **Rows or columns.** An SMD built by hand in MATLAB with row vectors loads like SMITE's
  columns.
- **Cells of numbers.** A per-localization field stored as a MATLAB cell of numbers (single
  and double values together, say) loads like a numeric column, each value checked and
  converted as stored (an integer in a cell is not rounded through single precision). A
  per-localization field that holds anything else (text, an empty cell element) is an error
  naming it, for the optional fields too; an empty one is a field left empty.
- **Complex values.** Localizations with a non-zero imaginary part in any position, photon or
  background field are dropped with a warning.
- **Number types.** Positions, photons, background and their uncertainties load as `Float64`,
  whatever the file stores; an integer above 2^53 in one of these fields (in a hand-built file,
  or one 0.1.0 wrote from integer data) can load rounded. Frames, datasets and track ids load as
  `Int`, exactly.

## Saving

```julia
mkpath("output/directory")
save_smite(smld, "output/directory", "processed_Results.mat")
```

`save_smite` takes an SMLD whose emitters have `x`, `y`, `photons`, `bg`, `σ_x`, `σ_y`,
`σ_photons`, `σ_bg`, `frame`, `dataset` and `track_id`, as `Emitter2DFit` and `Emitter3DFit` do
(plus `z` and `σ_z` for 3D), and whose camera has `pixel_edges_x` and `pixel_edges_y`. It can be
a [`SmiteSMLD`](@ref), a `BasicSMLD`, or your own type with `emitters`, `camera`, `n_frames`,
`n_datasets` and `metadata`. `Emitter2D` and `Emitter3D` cannot be saved. The folder must exist.
A missing folder, a folder (or anything but a file) of the file's name, a symlink to nothing, a
name that is no file, an existing file you may not write, too little free space to rewrite one
(where the system reports it), a camera whose pixels are not uniform and square (its edges,
relative to its first, are not within a thousandth of a pixel of those the loaders rebuild from
`XSize`, `YSize` and the written `PixelSize`, the differences of edges taken exactly), whose
edges are not all finite, or that
the loaders could not build back (no pixels along an axis, a first pixel of zero or negative
width, a pixel size that is not finite as written, or `Float32` edges that round onto each
other), and
metadata keys that are not MATLAB names (an ASCII letter, then up to 62 letters, digits or
underscores; Dict values and arrays of Dicts need such string keys) are each an
`ArgumentError` naming what is wrong, and nothing is written.

The file holds an `SMD` variable in SMITE's layout (MATLAB v7.3), which the loaders read. It
holds no `SMF`, so SMITE tools that also load an `SMF` need one from the original results file.
The file is written in full to a temporary file in the system's temporary folder (`tempdir()`,
which `TMPDIR` sets on Linux and macOS; it must exist and have room for one copy) and only then
copied into place. So a save that fails while writing (MAT
refusing a value, the file held open by HDF5 in the same session) leaves an existing file as it
was and creates none. An existing file keeps its owner, permissions and links, as with a direct
write (for a symlink, the file it points to is rewritten). An error during the final copy (a full
disk, a lost network share, an interruption) can leave the file partly written, as with a direct
write; the error then says the file is incomplete. The data are unchanged in memory, so save
again once the cause is fixed. It contains:

- `X`, `Y`, `X_SE`, `Y_SE` in SMITE's pixels, and `Photons`, `Bg`, `Photons_SE`, `Bg_SE`. Each
  column is converted by its values: all `Float32` stays MATLAB single, real values are
  written as double, and a column with a complex value (one with a non-zero imaginary part;
  otherwise it is its real part) as complex double, which the loaders screen as they do
  SMITE's complex fits. So `Float32` and `Float64` emitters together are written as double,
  never as the cell MAT would write for a mixed vector (which SMITE cannot compute with), and
  integers or `Float16` values (which MAT cannot write) become doubles. A value that is not a
  number, or a number that a double would round (a `BigFloat`, an integer above 2^53), is an
  `ArgumentError` naming the field and the emitter, and nothing is written. Positions are
  computed in `Float64` from the written pixel size (the camera's origin subtracted exactly),
  then stored as single for `Float32` data; a finite one whose value in pixels would load back
  as infinite, or a non-zero uncertainty that would load back as 0, is an `ArgumentError` too.
- `Z`, `Z_SE` in µm for `Emitter3DFit` emitters, or for non-empty 3D data whose every emitter
  has properties `z` and `σ_z`, converted as above; a `z` or `σ_z` a double would round is an
  `ArgumentError`, wherever it is in the column. 3D data whose `z` or `σ_z` are not all numbers
  (`missing`), and mixed 2D/3D data, save without `Z`, with a warning: an SMD has one `Z` column
  for all localizations, so save 2D and 3D data separately to keep it. The metadata's `Z` and
  `Z_SE` are then not saved either, so they never stand in for the emitters' `z`.
- `XSize`, `YSize` and `PixelSize` from the camera (`PixelSize` as single for a `Float32`
  camera, so the file loads back with the same camera), `NFrames`, `NDatasets` and `NDims`; counts
  and indices as doubles, as most SMITE results files hold them. They follow the loaders' rule,
  so every file `save_smite` writes loads: `n_frames` and `n_datasets` must be whole numbers of
  at least 0, and each emitter's `frame`, `dataset` and `track_id` whole numbers, each in
  `Int`'s range and held exactly by a double (every whole number up to 2^53, and some above).
  Any other value is an `ArgumentError` naming it, and nothing is written. SMITE has one pixel
  size, so a camera whose pixels are not uniform and square is an error.
- every metadata entry that is not an emitter column, except the loaders' own: `original_file`,
  `data_size`, `pixel_size`, `smite_defaulted_fields`, `smite_per_localization_fields`,
  `smite_row_frames`, `smite_row_photons`, `complex_fields_removed`, `complex_fields`,
  `original_emitter_count` and `removed_emitter_count`. Metadata keys must be MATLAB names: an
  ASCII letter, then up to 62 letters, digits or underscores. A Dict value needs string keys,
  each a MATLAB name. Otherwise `save_smite` throws an `ArgumentError` naming every offending
  key, before it writes anything.
- SMITE's per-localization fields a loader kept, for the saved emitters' rows, so after a filter
  `PValue` still belongs to its localization. Each keeps its MATLAB class (SMITE's `IndSMD` is
  a cell, `ThreshFlag` uint32), except `Z` and `Z_SE`, which `load_smite_2d` keeps from a 3D
  file: on 2D emitters they are written as `z` and `σ_z` are, above (a number a double would
  round is an `ArgumentError` naming the field and file row; a value that is not a number drops
  both, with a warning), and the file says `NDims` 3 when both are written and it has
  localizations, the loaded file's `Z` had values (SMITE writes `NDims` 2 for its own 3D fits)
  or the metadata's `NDims` is 3. `Z` and `Z_SE` in the metadata of data not loaded by SMLMio
  (built by hand, or loaded by SMLMData 0.7.1, which kept a 3D file's `Z` that way) are written
  the same way when they hold one value per emitter, in the order given, as 0.1.0 wrote them,
  so a save and a save of its load agree; with another count they are ordinary metadata. Such
  data carry no record of which localization each value belongs to (an `id` is not that: a
  join, a filter or frame connection keeps or reuses ids), so `save_smite` cannot check that
  order: it warns, naming every SMITE per-localization field (`Z`, `PValue`, ...) it writes in
  the order given. Sorting, filtering or joining such data before saving misplaces those
  values; load the file with SMLMio first, whose row record lets `save_smite` check them. An
  emitter of data SMLMio loaded is taken to be the row its `id` names
  when its frame and photons equal that row's. The fields are written only when every emitter is
  such a row; for combined localizations (frame connection's combined output), changed frames or
  photons, or emitters joined from another file, they are dropped with a warning. Emitters from
  another file whose frame and photons equal a row of this file (constant simulated photons, for
  example) are not detected and get this file's values. A field replaced since loading that no
  longer has one value per file row is dropped with a warning too. These fields are written by
  row because `metadata["smite_per_localization_fields"]` lists them; deleting that list makes
  `save_smite` write them as ordinary metadata, and after a filter they then no longer match the
  localizations; with one value per emitter they are then written in the order given, with the
  warning above.
- a field SMITE stacks per localization that is not listed in
  `metadata["smite_per_localization_fields"]` (the file held it with other than one value per
  row, often one per dataset: `PValue`, `ThreshFlag`, `XBoxCorner`, ... in some older files; or
  it was added after loading) as it is, unless the saved localizations are as many as its
  values: then a load would take it for one value per localization, and it is dropped with a
  warning. To save such a field by row, list it there with one value per file row. This applies
  to data loaded from SMITE; metadata on other data is written as given.
- an optional field the file left empty, which nobody has set since loading, back empty (an
  empty `ConnectID` stays empty). A `Photons_SE` or `Bg_SE` the file held with one value per
  dataset is saved as the `NaN` it loaded as; its per-dataset values are not saved. A
  `DatasetNum` or `ConnectID` the file held with one value per dataset is saved as the 1s or 0s
  it loaded as.

Loading the saved file gives back the emitters' positions, uncertainties, photons, background,
frames, datasets and track ids, and `n_frames` and `n_datasets`. It does not give back:

- `σ_xy`, and in 3D `σ_xz` and `σ_yz`: SMITE has no covariance fields; they load as 0.
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

# Positions are µm: keep bright localizations in a 10 µm square
bright = @filter(smld, photons > 1000)
roi = filter_roi(bright, (5.0, 15.0), (5.0, 15.0))

roi.metadata["analysis_date"] = string(Dates.now())
mkpath("processed_data")
save_smite(roi, "processed_data", "experiment1_processed_Results.mat")
```
