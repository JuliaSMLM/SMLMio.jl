"""
    save_smite(smld::AbstractSMLD, filepath::String, filename::String)

Save SMLD data to SMITE's SMD .mat format, readable by `load_smite_2d` and `load_smite_3d`.

# Arguments
- `smld::AbstractSMLD`: SMLD object to save (a `SmiteSMLD`, a `BasicSMLD` or any other SMLD
  with `emitters`, `camera`, `n_frames`, `n_datasets` and `metadata`)
- `filepath::String`: Directory path where to save the file
- `filename::String`: Name of the output .mat file

# Notes
- Saves in MATLAB v7.3 format
- Writes `X`, `Y`, `X_SE` and `Y_SE` in SMITE's pixels, (1, 1) the center of the camera's
  top-left pixel: `X = (x - x_edge) / pixel_size + 0.5`, `X_SE = σ_x / pixel_size`, where
  `x_edge` is the camera's first pixel edge. `Z` and `Z_SE` stay in microns.
- Writes `XSize`, `YSize` and `PixelSize` from the camera, and `NDims` (3 when it writes `Z`,
  else 2). SMITE has one pixel size, so a camera whose pixels are not uniform and square is an
  `ArgumentError`.
- Preserves all metadata fields whose names are not emitter columns, except the loaders' own
  `data_size`, `pixel_size` and `smite_defaulted_fields`; an emitter column wins over a metadata
  entry of the same name, so data that get `Z` and `Z_SE` here replace a metadata `Z` or `Z_SE`
  that v0.7.0 would have written
- Writes `Z` and `Z_SE` for an `Emitter3DFit` element type (as v0.7.0 did), or for non-empty 3D
  data (`emitter_ndims`) whose every emitter has properties `z` and `σ_z` and whose `z` and `σ_z`
  columns each come out as a `Vector{Float32}` or `Vector{Float64}`. Otherwise it writes v0.7.0's
  columns, with no `Z`, and warns when non-empty 3D or mixed 2D/3D data lose `Z` (a SMITE SMD
  file has one `Z` column for all localizations, so save 2D and 3D data separately to keep it).
  Apart from the camera check above, it never throws where v0.7.0 saved.
- Only `Float32` and `Float64` columns widen, because MATLAB stores those as single and double:
  MAT.jl writes `Float64`, `Float32`, `Int` and `Bool` vectors but throws on `Float16`,
  `Rational`, `BigFloat` and `missing`, and writes a mixed `Float32`/`Float64` column as a cell.
- A type whose `propertynames` lists a property that `getproperty` refuses breaks Julia's
  property interface and is not supported.
"""
function save_smite(smld::AbstractSMLD, filepath::String, filename::String)
    # Create SMD structure
    s = Dict{String,Any}()

    # Z columns for v0.7.0's rule (an Emitter3DFit element type), or for non-empty 3D data whose
    # z and σ_z columns MAT writes as MATLAB single or double. Anything else gets v0.7.0's columns
    # and never throws where v0.7.0 saved; non-empty 3D or mixed data that lose Z get a warning.
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
    
    # Counts and indices are doubles, as SMITE writes them (MATLAB integer arithmetic rounds)
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

    # Optional fields the loader filled with defaults and that still hold only the default go back
    # as the file had them (SMITE treats an empty ConnectID differently from zeros).
    defaulted = get(md, "smite_defaulted_fields", String[])
    for (key, default) in _SMITE_OPTIONAL
        key in defaulted && haskey(md, key) && all(v -> isequal(v, default), s[key]) && (s[key] = md[key])
    end

    # Fields the loader kept with one value per file row (PValue, LogLikelihood, ...) are written
    # for these emitters' rows, the file row being each emitter's id. Combined data sets, or ids
    # that are not distinct file rows, cannot be traced: those fields are dropped with a warning.
    perloc = filter(k -> haskey(md, k), get(md, "smite_per_localization_fields", String[]))
    rows = Int[]
    traceable = isempty(perloc)
    if !traceable && !haskey(md, "concatenated_from") && !haskey(md, "merged_from") &&
       all(e -> hasproperty(e, :id), smld.emitters)
        rows = Int[e.id for e in smld.emitters]
        traceable = allunique(rows) && all(k -> all(r -> 1 <= r <= _smite_length(md[k]), rows), perloc)
    end
    traceable || @warn("save_smite: $(join(perloc, ", ")) not saved: their rows cannot be traced to " *
                       "these emitters (combined data sets, or emitter ids that are not file rows)")

    # Add any additional fields from metadata, except the loaders' own entries
    for (key, value) in md
        (key in _SMITE_LOADER_KEYS || haskey(s, key)) && continue
        if key in perloc
            traceable || continue
            value = _smite_vec(value)[rows]
        end
        s[key] = value
    end
    
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
