# Changelog

## 0.1.1 (unreleased)

Fixed in 0.1.0's SMITE I/O; each fix has a test that fails without it:

- Scalar fields were converted to `Float64` before they were checked. A complex `NFrames`,
  `NDatasets`, `XSize`, `YSize` or `PixelSize` loaded as its real part, counts above 2^53 were
  rounded, and a MATLAB single `PixelSize` gave an `IdealCamera{Float64}` (SMLMData 0.7.1 gave
  an `IdealCamera{Float32}`). Each is now checked as stored, and counts load exactly. A single
  `PixelSize` gives an `IdealCamera{Float32}` again (any other pixel size a `Float64` camera),
  and `save_smite` writes a `Float32` camera's `PixelSize` as single, so the file loads back
  with the same camera; a `Float32` pixel size whose camera edges overflow or round onto each
  other, or a pixel size that is not finite, is an `ArgumentError`. `XSize`
  and `YSize` must be whole numbers of at least 1 (0 or below loaded a camera with no pixels),
  `NFrames` and `NDatasets` whole numbers of at least 0; a value with a non-zero imaginary part,
  text or several values is an `ArgumentError` naming the field, and a `pixel_size` passed to
  the loader no longer stands in for a `PixelSize` of several values or text. A `FrameNum`,
  `DatasetNum` or `ConnectID` value too large for an `Int` is an `ArgumentError`, not an
  `InexactError`. In 3D, `ZSize` (which SMITE does not write) now stays in the metadata, so a
  save and load keeps it, and is never an error: one whole number of at least 0 that fits an
  `Int` is appended to `data_size`, anything else is only kept (0.1.0 refused a `ZSize` that was
  not whole).
- `save_smite` wrote `n_frames`, `n_datasets` and each emitter's `frame`, `dataset` and
  `track_id` as doubles without checking them, so a value a double cannot hold exactly was saved
  rounded, and values the loaders refuse (a negative `n_frames`, a frame of 2.5) were written.
  They now follow the loaders' rule: whole numbers in `Int`'s range that a double holds exactly,
  with `n_frames` and `n_datasets` at least 0. Any other value is an `ArgumentError` naming it,
  and nothing is written.
- The loaders refused a per-localization field stored as a MATLAB cell of numbers (MAT reads it
  as an array of `Any`), which SMLMData 0.7.1 loaded: "Photons is not numeric". Such a field
  now loads as a numeric column, each value checked and converted as stored, so an integer
  frame or track id in a cell is never rounded. An optional field (`Photons_SE`, `Bg_SE`,
  `DatasetNum`, `ConnectID`) holding something other than numbers, text included, is an
  `ArgumentError` naming it, where it failed later with a `MethodError` or, for text, loaded
  as a column of defaults. Empty text, which MAT 0.10 and 0.11 read for an empty MATLAB char
  field (`''`), is made an empty array in SMITE's own fields as the file is read, so it is an
  empty field everywhere, as an empty array is (other fields keep their text); an empty-text
  `ConnectID` set by hand is written back empty.
- `save_smite` of emitters of `Float32` and `Float64` together (a `Union` element type) wrote
  `Photons`, `Bg`, their errors and `Z`, `Z_SE` as MATLAB cells, which MAT makes of a vector of
  mixed types, so `load_smite_3d` refused its own file; a 3D file with such cells, loaded in 2D
  and saved, kept its cell `Z`. Every emitter column is now converted by its values: all
  `Float32` stays single, real values are written as double when a double holds every value
  exactly (integers, `Float16`; a complex value with a zero imaginary part is its real part),
  complex values as complex double (as 0.1.0 wrote them; the loaders screen them), and a value
  that is not a number or a number a double would round is an `ArgumentError` naming the field
  and the emitter (for `z` and `σ_z` that are not numbers, such as `missing`, the file is saved
  without `Z`, with a warning, as before, and so is the metadata's). Positions in SMITE's pixels are
  computed in `Float64` from the written pixel size, where a `Float32` camera rounded them, and
  each finite one must load back finite (a non-zero uncertainty also not as 0); `save_smite`
  refuses a camera whose edges are not, within a thousandth of a pixel, those its loaders rebuild
  from the written pixel size (the check rounded in the camera's own type, so a `Float32` camera
  far from the origin or a `Float16` one passed and loaded back with other edges), or that they
  could not build back (no pixels along an axis, a pixel size that is not finite as written,
  `Float32` edges that round onto each other). 3D
  data whose `z` mixes `Float32` and `Float64`, or is `Float16`, save their `Z` instead of
  dropping it. `Z` and `Z_SE` kept in the metadata by `load_smite_2d` are written the same way
  for 2D emitters, as are metadata `Z` and `Z_SE` with one value per emitter on data not loaded
  by SMLMio (built by hand, or loaded by SMLMData 0.7.1), so a save and a save of its load agree;
  the file says `NDims` 3 when both are written and it has localizations, its `Z` had values
  (SMITE writes `NDims` 2 for its 3D fits) or the metadata's `NDims` is 3. Every other
  per-localization field keeps its MATLAB class (SMITE's `IndSMD` stays a cell).

## 0.1.0 (2026-11-16)

- Created from SMLMData 0.7.1's SMITE I/O: `SmiteSMD`, `SmiteSMLD`, `load_smite_2d`,
  `load_smite_3d`, and `save_smite`, with git history kept.
- Depends on SMLMCore and MAT; requires Julia 1.10 or later (SMLMCore's floor).
- Weak dependency on SMLMData >= 0.7.2 guards against two copies of each type
  beside pre-split SMLMData: an environment that holds SMLMData 0.7.1 or earlier cannot add
  SMLMio.

Fixed relative to SMLMData 0.7.1. The positions, optional-field, empty-result, save-and-load
and filter-then-save fixes were reproduced on real SMITE results files first; the others were
found by reading SMITE's code or by review, and each has a test that fails without it:

- Positions: SMITE stores `X`, `Y`, `X_SE`, `Y_SE` in pixels with (1, 1) the center of the
  top-left pixel. The loaders kept them as microns (every file loaded 1/PixelSize times too large,
  outside its camera); they now load `x = (X - 0.5) * PixelSize`, `σ_x = X_SE * PixelSize`, and
  `save_smite` writes the inverse. `Z` and `Z_SE` are microns in both and unchanged.
- Files whose `Photons_SE` and `Bg_SE` hold one value per dataset or are empty, or whose
  `ConnectID` is empty, crashed with a `BoundsError`. Those fields now load
  as `NaN` (and `ConnectID` as 0, `DatasetNum` as 1), raw values kept in the metadata and listed
  in `metadata["smite_defaulted_fields"]`; a non-empty field of the wrong length warns. A
  required field of the wrong length is an `ArgumentError` naming it. Row-vector SMDs load every
  localization, not only the first.
- Pixel size: a file without `PixelSize` silently loaded at 0.1 µm per pixel. The loaders now
  take a `pixel_size` keyword for such files and throw without one; a `pixel_size` that differs
  from the file's is an error.
- Results without localizations (`XSize`, `YSize` empty) failed with a `MethodError`; they now
  throw an `ArgumentError` naming the missing image size.
- `load_smite_3d` required a `ZSize` field that SMITE never writes.
- `save_smite` files could not be read back (no `XSize`, `YSize`, `PixelSize`; the loaders'
  `data_size` and `pixel_size` written as SMD fields). It now writes those from the camera, and
  `NDims`, and skips the loaders' own metadata; a camera whose pixels are not uniform and square
  is an `ArgumentError`. `save_smite` takes any `AbstractSMLD`, not only a `SmiteSMLD`.
- After a filter, or complex screening at load, `save_smite` wrote the SMD fields kept in
  metadata (`PValue`, `LogLikelihood`, ...) unchanged next to fewer localizations, so they
  belonged to other localizations. The loaders list SMITE's per-localization fields (catSMD's
  list) in `metadata["smite_per_localization_fields"]` and record each row's frame and photons;
  `save_smite` writes those fields by each emitter's file row (`id`) when every emitter is still
  its row (its frame and photons equal that row's), and otherwise (combined localizations such
  as frame connection's combined output, changed frames or photons, another file's emitters)
  drops them with a warning. Emitters from another file whose frame and photons equal a row of
  this file (constant simulated photons, for example) are not detected and get this file's
  values.
- Some older SMITE files hold `PValue`, `ThreshFlag`, `XBoxCorner` or `YBoxCorner` with one value
  per dataset. Saved with as many localizations as values, `save_smite` wrote them where a load
  (and SMITE) took them for one value per localization; such a field is now dropped with a
  warning in that case.
- `save_smite` wrote the loaders' own metadata (`original_file`, `complex_*`) into the file, and
  loading it back reported the earlier file's name and complex screening. The loaders' entries
  are now never saved and never taken from a file.
- `save_smite` wrote counts and indices (`FrameNum`, `DatasetNum`, `ConnectID`, `NFrames`,
  `NDatasets`) as MATLAB int64; most SMITE results files hold doubles, and MATLAB integer
  arithmetic rounds. They are now doubles. Optional fields the file left empty and nobody has set
  since loading are saved empty (an empty `ConnectID` stays empty).
- The loaders left the file open when reading it failed, and a file without the variable gave a
  bare `KeyError`; the file is now closed on any error and a missing variable is an
  `ArgumentError` naming it.
- `load_smite_2d` on a 3D file now warns that it leaves out `Z`; a camera whose first edge is
  not at 0 now warns on save, since SMITE has no origin.
- `save_smite` threw from inside MAT for a metadata key that is not a MATLAB name or a Dict value
  with non-string keys, after MAT had truncated the target, so saving over an existing file
  destroyed it. It now checks every key (every nested Dict's and struct array's too) first and
  throws one `ArgumentError` naming all offending keys. The file is written in full to a
  temporary file before the target is touched and then copied into place, so a save that fails
  while writing leaves an existing file as it was (it keeps its owner, permissions and links);
  only an error during the final copy, such as a full disk, can leave it partly written, and the
  error then says the file is incomplete and to save again. A
  missing folder, a folder at the target, a symlink to nothing, a name that is no file, an
  existing file that may not be written and too little free space to rewrite one are
  `ArgumentError`s before anything is written.
- A `FrameNum`, `DatasetNum` or `ConnectID` that is not a whole number threw `InexactError`; it
  is now an `ArgumentError` naming the field (`DatasetNum` and `ConnectID` are checked on the
  loaded rows only). So are `XSize`, `YSize`, `NFrames`, `NDatasets` and `ZSize`. An SMD whose
  `NFrames` or `NDatasets` holds a list of values was reported as having none; the error now says
  how many values the file holds.
- `SmiteSMD` and `save_smite` take any string type for paths (a `SubString`, for example), not
  only `String`. `SmiteSMD` shows its whole folder path, so its compact form reads as the
  constructor call that rebuilds it.
- Documentation: `n_frames` is frames per dataset; `save_smite` needs fit emitters
  (`Emitter2DFit` or `Emitter3DFit`); `σ_xz` and `σ_yz` are not saved either; the file holds an
  `SMD` in SMITE's layout, which the loaders read, and no `SMF`, so SMITE tools that also load
  `SMF` need one from the original results file.
