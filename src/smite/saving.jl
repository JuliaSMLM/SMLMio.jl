"""
    save_smite(smld::AbstractSMLD, filepath::String, filename::String)

Save SMLD data to SMITE's SMD .mat format (MATLAB v7.3), readable by `load_smite_2d`,
`load_smite_3d` and SMITE.

# Arguments
- `smld::AbstractSMLD`: SMLD object to save (a `SmiteSMLD`, a `BasicSMLD` or any other SMLD
  with `emitters`, `camera`, `n_frames`, `n_datasets` and `metadata`)
- `filepath::String`: Directory path where to save the file
- `filename::String`: Name of the output .mat file

# What it writes
- `X`, `Y`, `X_SE`, `Y_SE` in SMITE's pixels, (1, 1) the center of the camera's top-left pixel:
  `X = (x - x_edge) / pixel_size + 0.5`, `X_SE = σ_x / pixel_size`, with `x_edge` the camera's
  first pixel edge, in the emitters' float type. SMITE has no origin: a camera whose first edge
  is not at 0 warns, and the file loads back with the camera at 0 and the positions moved by
  that offset.
- `Photons`, `Bg`, their standard errors, `FrameNum`, `DatasetNum` and `ConnectID` (the track
  id). Counts and indices are written as doubles, as SMITE's results files hold them (SMITE
  creates an empty `ConnectID` as uint32; MATLAB integer arithmetic would round).
- `Z` and `Z_SE` (microns) for `Emitter3DFit` emitters, or for non-empty 3D data
  (`emitter_ndims`) whose every emitter has properties `z` and `σ_z` and whose `z` and `σ_z`
  columns are each all `Float32` or all `Float64` (MATLAB single or double; MAT.jl cannot write
  `Float16`, `Rational`, `BigFloat` or `missing`, and writes a mixed column as a cell). Non-empty
  3D or mixed 2D/3D data without such columns save without `Z`, with a warning: an SMD has one
  `Z` column for all localizations, so save 2D and 3D data separately to keep it.
- `XSize`, `YSize` and `PixelSize` from the camera, `NFrames`, `NDatasets`, and `NDims` (3 with
  `Z`, else 2). SMITE has one pixel size, so a camera whose pixels are not uniform and square is
  an `ArgumentError`. An sCMOS camera's calibration maps are not saved; the file loads back
  with an `IdealCamera`.
- Every metadata entry that is not an emitter column, except the loaders' own entries:
  `original_file`, `data_size`, `pixel_size`, `smite_defaulted_fields`,
  `smite_per_localization_fields`, `smite_row_frames`, `smite_row_photons`,
  `complex_fields_removed`, `complex_fields`, `original_emitter_count` and
  `removed_emitter_count`.
- SMITE's per-localization fields a loader kept (`metadata["smite_per_localization_fields"]`:
  `PValue`, `LogLikelihood`, `PSFSigma`, ...) for the saved emitters' file rows, the row being
  each emitter's `id`, when every emitter is still the row its `id` names: the same frame and
  photons the loader recorded for that row (`metadata["smite_row_frames"]`,
  `metadata["smite_row_photons"]`). Otherwise (combined localizations such as frame
  connection's combined output, emitters joined from another file, changed frames or photons,
  ids that are not rows, no row record) those fields are dropped with a warning; so is a field
  replaced since loading that no longer has one value per file row.
- An optional field the file left empty (`ConnectID` before frame connection, ...), which the
  loader defaulted, back empty while it still holds only the default.

Not saved: `σ_xy` (SMITE has no such field; it loads back as 0) and the emitter `id`, which
loads back as the row in the new file. A type whose `propertynames` lists a property that
`getproperty` refuses breaks Julia's property interface and is not supported.
"""
function save_smite(smld::AbstractSMLD, filepath::String, filename::String)
    # Create SMD structure
    s = Dict{String,Any}()

    # Z columns for an Emitter3DFit element type, or for non-empty 3D data whose z and σ_z columns
    # MAT writes as MATLAB single or double. Anything else saves without Z; non-empty 3D or mixed
    # data that lose Z get a warning.
    d = emitter_ndims(smld.emitters)
    zcols = nothing
    if eltype(smld.emitters) <: Emitter3DFit
        zcols = ([e.z for e in smld.emitters], [e.σ_z for e in smld.emitters])
    elseif d == 3 && !isempty(smld.emitters) &&
           all(e -> hasproperty(e, :z) && hasproperty(e, :σ_z), smld.emitters)
        z, σz = [e.z for e in smld.emitters], [e.σ_z for e in smld.emitters]
        eltype(z) in (Float32, Float64) && eltype(σz) in (Float32, Float64) && (zcols = (z, σz))
    end
    if zcols === nothing && !isempty(smld.emitters)
        d === nothing && @warn("save_smite: mixed 2D/3D emitters: Z not saved; " *
                               "save 2D and 3D data separately")
        d == 3 && @warn("save_smite: 3D emitters without Float32/Float64 z and σ_z: Z not saved")
    end
    
    # Extract arrays from emitters. SMITE stores X, Y and their standard errors in pixels, (1, 1)
    # the center of the top-left pixel; positions here are microns from the camera's corner.
    # SMITE has one PixelSize, so the camera's pixels must be uniform and square: every edge lies
    # on the first pixel's grid to a thousandth of a pixel (edges carry floating rounding, which
    # in a Float32 camera grows with distance from the origin).
    cam = smld.camera
    ex, ey = cam.pixel_edges_x, cam.pixel_edges_y
    x0, y0 = ex[1], ey[1]
    px = ex[2] - ex[1]
    on_grid(edges, e0) = all(k -> abs(edges[k] - (e0 + (k - 1) * px)) <= px / 1000, eachindex(edges))
    px > 0 && on_grid(ex, x0) && on_grid(ey, y0) ||
        throw(ArgumentError("save_smite: SMITE stores one square pixel size; this camera's pixels " *
                            "are not uniform and square"))
    iszero(x0) && iszero(y0) || @warn("save_smite: the camera's first pixel edge is at ($x0, $y0) µm; " *
        "SMITE has no origin, so the file loads back with the camera at 0 and positions moved by that offset")
    s["X"] = _same_float([e.x for e in smld.emitters], e -> (e.x - x0) / px + 0.5, smld.emitters)
    s["Y"] = _same_float([e.y for e in smld.emitters], e -> (e.y - y0) / px + 0.5, smld.emitters)
    if zcols !== nothing
        s["Z"] = zcols[1]
    end
    
    s["Photons"] = [e.photons for e in smld.emitters]
    s["Bg"] = [e.bg for e in smld.emitters]
    
    s["X_SE"] = _same_float([e.σ_x for e in smld.emitters], e -> e.σ_x / px, smld.emitters)
    s["Y_SE"] = _same_float([e.σ_y for e in smld.emitters], e -> e.σ_y / px, smld.emitters)
    if zcols !== nothing
        s["Z_SE"] = zcols[2]
    end
    
    s["Photons_SE"] = [e.σ_photons for e in smld.emitters]
    s["Bg_SE"] = [e.σ_bg for e in smld.emitters]
    
    # Counts and indices are doubles, as SMITE's results files hold them (MATLAB integer arithmetic
    # rounds)
    s["FrameNum"] = Float64[e.frame for e in smld.emitters]
    s["DatasetNum"] = Float64[e.dataset for e in smld.emitters]
    s["ConnectID"] = Float64[e.track_id for e in smld.emitters]

    # Add metadata: what the loaders need to read the file back, from the camera
    s["NFrames"] = Float64(smld.n_frames)
    s["NDatasets"] = Float64(smld.n_datasets)
    s["XSize"] = Float64(length(ex) - 1)
    s["YSize"] = Float64(length(ey) - 1)
    s["PixelSize"] = Float64(px)
    s["NDims"] = zcols === nothing ? 2.0 : 3.0

    md = smld.metadata

    # Optional fields the file left empty, which the loader filled with defaults and that still hold
    # only the default, go back empty (SMITE treats an empty ConnectID differently from zeros). A
    # field the file held with one value per dataset is written as the default column instead: its
    # raw values are not one per localization, and a later load could take them for that.
    defaulted = get(md, "smite_defaulted_fields", String[])
    for (key, default) in _SMITE_OPTIONAL
        key in defaulted && haskey(md, key) && _smite_length(md[key]) == 0 &&
            all(v -> isequal(v, default), s[key]) && (s[key] = md[key])
    end

    # SMITE's per-localization fields the loader kept (PValue, LogLikelihood, ...) are written for
    # these emitters' file rows, the row being each emitter's id, when every emitter is still its
    # row: the same frame and photons the loader recorded for that row (smite_row_frames,
    # smite_row_photons). Otherwise (combined localizations, emitters from another file, changed
    # frames or photons, ids that are not rows) those fields are dropped with a warning, never
    # misplaced; so is a field replaced since loading that no longer has one value per row.
    perloc = filter(k -> haskey(md, k), get(md, "smite_per_localization_fields", String[]))
    rowframes = get(md, "smite_row_frames", nothing)
    rowphotons = get(md, "smite_row_photons", nothing)
    rows = Int[]
    towrite = String[]
    if !isempty(perloc) && (rowframes === nothing || rowphotons === nothing)
        @warn("save_smite: $(join(perloc, ", ")) not saved: the metadata has no record of the file's rows")
    elseif !isempty(perloc)
        nrows = length(rowframes)
        changed = filter(k -> _smite_length(md[k]) != nrows, perloc)
        isempty(changed) || @warn("save_smite: $(join(changed, ", ")) no longer holds one value per " *
                                  "file row ($nrows); not saved")
        remaining = setdiff(perloc, changed)
        isrow(e) = hasproperty(e, :id) && (r = e.id; r isa Integer && 1 <= r <= nrows &&
                   e.frame == rowframes[r] && isequal(e.photons, rowphotons[r]))
        if all(isrow, smld.emitters)
            rows = Int[e.id for e in smld.emitters]
            towrite = remaining
        elseif !isempty(remaining)
            @warn("save_smite: $(join(remaining, ", ")) not saved: these emitters are not all rows of the " *
                  "file they were loaded from (combined localizations, another file's emitters, changed " *
                  "frames or photons, or ids that are not rows)")
        end
    end

    # Add any additional fields from metadata, except the loaders' own entries
    for (key, value) in md
        (key in _SMITE_LOADER_KEYS || haskey(s, key)) && continue
        if key in perloc
            key in towrite || continue
            value = _smite_vec(value)[rows]
        end
        s[key] = value
    end
    # A 3D file loaded as 2D keeps Z in metadata, written back above for each emitter: say 3D
    zcols === nothing && !isempty(smld.emitters) && haskey(s, "Z") &&
        _smite_length(s["Z"]) == length(smld.emitters) && (s["NDims"] = 3.0)
    
    # Save to file
    matwrite(joinpath(filepath, filename), Dict("SMD" => s))
end

# The column f(e) for each emitter, in the float type of the emitters' own column `col` when that
# is Float32 or Float64 (MATLAB single or double), so Float32 data stay single.
function _same_float(col, f, emitters)
    vals = [f(e) for e in emitters]
    T = eltype(col)
    return T in (Float32, Float64) ? convert(Vector{T}, vals) : vals
end
