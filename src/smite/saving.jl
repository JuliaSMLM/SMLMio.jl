"""
    save_smite(smld::AbstractSMLD, filepath::AbstractString, filename::AbstractString) -> nothing

Save an SMLD as an `SMD` variable in a SMITE `.mat` file (MATLAB v7.3), which the loaders read.

Positions convert from µm to SMITE's pixels, (1, 1) the center of the camera's top-left pixel:
`X = (x - x_edge) / pixel_size + 0.5` and `X_SE = σ_x / pixel_size` (the same for `y`), with
`x_edge` the camera's first pixel edge. `Z` and `Z_SE` are µm in both. The file holds no `SMF`.
Which fields are written, which are dropped and what a load gives back are on the
[SMITE files](@ref) page.

# Arguments
- `smld::AbstractSMLD`: an SMLD (`SmiteSMLD`, `BasicSMLD`, or your own with `emitters`, `camera`,
  `n_frames`, `n_datasets` and `metadata`) whose emitters have `x`, `y`, `photons`, `bg`, `σ_x`,
  `σ_y`, `σ_photons`, `σ_bg`, `frame`, `dataset` and `track_id`, as `Emitter2DFit` and
  `Emitter3DFit` do (plus `z` and `σ_z` for 3D), and whose camera has `pixel_edges_x` and
  `pixel_edges_y`. `Emitter2D` and `Emitter3D` cannot be saved.
- `filepath::AbstractString`: an existing folder.
- `filename::AbstractString`: name of the `.mat` file; an existing file of that name is
  rewritten.

# Returns
- `nothing`.

# Notes
- A failed save leaves an existing file as it was and no new file behind; what it refuses and
  how a file is rewritten are on the [SMITE files](@ref) page.

- `σ_xy`, and in 3D `σ_xz` and `σ_yz`, are not saved: SMITE has no covariance fields, and they
  load as 0.

# Example
```jldoctest
julia> using SMLMCore: IdealCamera, Emitter2DFit, BasicSMLD

julia> smld = BasicSMLD([Emitter2DFit{Float64}(1.0, 2.0, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0)], IdealCamera(64, 64, 0.1), 1, 1);

julia> dir = mktempdir(); save_smite(smld, dir, "localizations.mat")

julia> e = load_smite_2d(SmiteSMD(dir, "localizations.mat")).emitters[1]; (e.x, e.y)
(1.0, 2.0)
```

# See also
[`load_smite_2d`](@ref), [`load_smite_3d`](@ref), [`SmiteSMD`](@ref)
"""
function save_smite(smld::AbstractSMLD, filepath::AbstractString, filename::AbstractString)
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
    
    # Counts and indices are doubles, as most SMITE results files hold them (MATLAB integer arithmetic
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

    # Add any additional fields from metadata, except the loaders' own entries. For data loaded
    # from SMITE, a field SMITE stacks per localization that is not listed as one value per file
    # row (the file held it otherwise, often one value per dataset, or it was added after loading)
    # is dropped when the saved localizations are as many as its values, since a load would take
    # it for one value per localization; with another count it is written as it is.
    loaded = haskey(md, "smite_per_localization_fields")
    nsaved = length(smld.emitters)
    ambiguous = String[]
    for (key, value) in md
        (key in _SMITE_LOADER_KEYS || haskey(s, key)) && continue
        if key in perloc
            key in towrite || continue
            value = _smite_vec(value)[rows]
        elseif loaded && key in _SMITE_STACKABLE && nsaved > 0 && _smite_length(value) == nsaved
            push!(ambiguous, key)
            continue
        end
        s[key] = value
    end
    isempty(ambiguous) || @warn("save_smite: $(join(sort!(ambiguous), ", ")) not saved: SMITE stores these " *
        "per localization, they are not listed in smite_per_localization_fields as one value per file " *
        "row, and with as many values as saved localizations a load would take them for one value per " *
        "localization; list a field there, with one value per file row, to save it by row")
    # A 3D file loaded as 2D keeps Z in metadata, written back above for each emitter: say 3D
    zcols === nothing && !isempty(smld.emitters) && haskey(s, "Z") &&
        _smite_length(s["Z"]) == length(smld.emitters) && (s["NDims"] = 3.0)
    
    # MAT truncates the target before it checks a struct's keys. Check them all first; an
    # existing file is rewritten only once the new contents are complete.
    _check_mat_keys(s)
    _write_smd(joinpath(String(filepath), String(filename)), s)
    return nothing
end

# Writes the SMD to `target`, never removing anything that was there: a target that is a folder,
# or a name that is no file, is an ArgumentError. A new file is written directly and removed if
# the write fails. An existing file is first written in full to a temporary file elsewhere, and
# only then copied into the file in place, so a failed write leaves it as it was, and the file
# keeps its owner, group, permissions and links, as a direct write does (a symlink's target file
# is the one rewritten). Whether the file may be written is the system's answer, as for a direct
# write. Only an interruption during the final copy can leave the file partly written.
function _write_smd(target::String, s)
    dir, base = splitdir(target)
    (isempty(base) || base in (".", "..")) &&
        throw(ArgumentError("save_smite: $(repr(target)) names no file"))
    isempty(dir) && (dir = ".")
    isdir(dir) || throw(ArgumentError(ispath(dir) ? "save_smite: $(repr(dir)) is not a folder" :
                                      "save_smite: the folder $(repr(dir)) does not exist"))
    isdir(target) && throw(ArgumentError("save_smite: $(repr(target)) is a folder"))
    if !isfile(target)
        try
            matwrite(target, Dict("SMD" => s))
        catch
            # Remove the partial new file, never a symlink that was there
            islink(target) || try rm(target; force=true) catch end
            rethrow()
        end
        return nothing
    end
    try
        close(open(target, "r+"))
    catch e
        e isa SystemError || rethrow()
        throw(ArgumentError("save_smite: $(repr(target)) cannot be written: $(Libc.strerror(e.errnum))"))
    end
    tmp = tempname()
    try
        matwrite(tmp, Dict("SMD" => s))
        open(tmp, "r") do src
            open(target, "w") do dst   # truncates in place: owner, group, mode and links stay
                buf = Vector{UInt8}(undef, 1 << 20)
                while !eof(src)
                    n = readbytes!(src, buf)
                    write(dst, view(buf, 1:n))
                end
            end
        end
    finally
        rm(tmp; force=true)
    end
    return nothing
end

# MATLAB's names: an ASCII letter, then up to 62 letters, digits or underscores.
_is_mat_name(k) = k isa AbstractString && occursin(r"\A[A-Za-z][A-Za-z0-9_]{0,62}\z", k)

# Throws one ArgumentError naming every key MAT would refuse: a struct key that is not a MATLAB
# name, or a key of a Dict value that is not a string. `path` names the Dict inside the SMD.
function _check_mat_keys(s)
    bad = String[]
    _collect_bad_keys!(bad, "", s)
    isempty(bad) || throw(ArgumentError("save_smite: MATLAB cannot store these keys; a key must be " *
        "a string of an ASCII letter, then up to 62 letters, digits or underscores: " *
        join(bad, ", ")))
end

function _collect_bad_keys!(bad, path, d::AbstractDict)
    for (k, v) in d
        _is_mat_name(k) || push!(bad, isempty(path) ? "metadata key $(repr(k))" :
                                      "key $(repr(k)) of the Dict $(path)")
        name = isempty(path) ? string(k) : string(path, ".", k)
        _collect_bad_keys!(bad, name, v)
    end
end
_collect_bad_keys!(bad, path, _) = nothing

# The column f(e) for each emitter, in the float type of the emitters' own column `col` when that
# is Float32 or Float64 (MATLAB single or double), so Float32 data stay single.
function _same_float(col, f, emitters)
    vals = [f(e) for e in emitters]
    T = eltype(col)
    return T in (Float32, Float64) ? convert(Vector{T}, vals) : vals
end
