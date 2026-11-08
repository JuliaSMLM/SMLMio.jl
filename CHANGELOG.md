# Changelog

## 0.1.0 (unreleased)

- Created from SMLMData 0.7.1's SMITE I/O: `SmiteSMD`, `SmiteSMLD`, `load_smite_2d`,
  `load_smite_3d`, and `save_smite`, with git history kept.
- Depends on SMLMCore and MAT.
- Weak dependency on SMLMData >= 0.7.2 guards against two copies of each type
  beside pre-split SMLMData.

Fixed relative to SMLMData 0.7.1:

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
  belonged to other localizations. The loaders list those fields in
  `metadata["smite_per_localization_fields"]` and `save_smite` writes them by each emitter's file
  row (`id`); for combined data sets it drops them with a warning.
- `save_smite` wrote the loaders' own metadata (`original_file`, `complex_*`) into the file, and
  loading it back reported the earlier file's name and complex screening. The loaders' entries
  are now never saved and never taken from a file.
- `save_smite` wrote counts and indices (`FrameNum`, `DatasetNum`, `ConnectID`, `NFrames`,
  `NDatasets`) as MATLAB int64; SMITE writes doubles, and MATLAB integer arithmetic rounds. They
  are now doubles. Defaulted fields that nobody has set since loading are saved as the file had
  them (an empty `ConnectID` stays empty).
- `load_smite_2d` on a 3D file now warns that it leaves out `Z`; a camera whose first edge is
  not at 0 now warns on save, since SMITE has no origin.
