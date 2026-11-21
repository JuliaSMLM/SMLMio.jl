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
- The file is written in full to a temporary file before the target is touched, so a save that
  fails while writing leaves an existing file as it was; only an error during the final copy
  (a full disk, a lost share) can leave it partly written. What it refuses, and how a file is
  rewritten, are on the [SMITE files](@ref) page.

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

    # Z columns for an Emitter3DFit element type, or for non-empty 3D data whose every emitter has
    # z and σ_z, when every z and σ_z is a real number a double holds exactly (_mat_floats).
    # Anything else saves without Z; non-empty 3D or mixed data that lose Z get a warning.
    d = emitter_ndims(smld.emitters)
    zcols = nothing
    if eltype(smld.emitters) <: Emitter3DFit || (d == 3 && !isempty(smld.emitters) &&
           all(e -> hasproperty(e, :z) && hasproperty(e, :σ_z), smld.emitters))
        zv, σzv = [e.z for e in smld.emitters], [e.σ_z for e in smld.emitters]
        z, σz = _mat_floats(zv), _mat_floats(σzv)
        # A real z or σ_z a double would round is an error naming the emitter; z that are not
        # real numbers (missing) save without Z, with the warning below, as before
        for (name, vals, r) in (("z", zv, z), ("σ_z", σzv, σz))
            r isa Integer && _is_real(vals[r]) && _mat_column(name, vals)
        end
        z isa Integer || σz isa Integer || (zcols = (z, σz))
    end
    if zcols === nothing && !isempty(smld.emitters)
        d === nothing && @warn("save_smite: mixed 2D/3D emitters: Z not saved; " *
                               "save 2D and 3D data separately")
        d == 3 && @warn("save_smite: 3D emitters whose z or σ_z are not all real numbers: Z not saved")
    end
    
    # Extract arrays from emitters. SMITE stores X, Y and their standard errors in pixels, (1, 1)
    # the center of the top-left pixel; positions here are microns from the camera's corner.
    # SMITE has one PixelSize, so the camera's pixels must be uniform and square: every edge lies
    # on the first pixel's grid to a thousandth of a pixel (edges carry floating rounding, which
    # in a Float32 camera grows with distance from the origin).
    cam = smld.camera
    ex, ey = cam.pixel_edges_x, cam.pixel_edges_y
    nx, ny = length(ex) - 1, length(ey) - 1
    nx >= 1 && ny >= 1 || throw(ArgumentError("save_smite: the camera has $nx × $ny pixels; SMITE's " *
                                              "XSize and YSize must each be at least 1"))
    x0, y0 = ex[1], ey[1]
    px = ex[2] - ex[1]
    on_grid(edges, e0) = all(k -> abs(edges[k] - (e0 + (k - 1) * px)) <= px / 1000, eachindex(edges))
    px > 0 && on_grid(ex, x0) && on_grid(ey, y0) ||
        throw(ArgumentError("save_smite: SMITE stores one square pixel size; this camera's pixels " *
                            "are not uniform and square"))
    # The pixel size as written, single for a Float32 camera (so the file loads back with a Float32
    # camera) and otherwise double; the camera the loaders will build from it must be valid
    pxw = eltype(ex) === Float32 ? px : Float64(px)
    _smite_camera("save_smite", nx, ny, _smite_pixel_size("save_smite", "the camera's pixel size", pxw))
    # Positions in SMITE's pixels, computed in Float64 from the written pixel size
    x0f, y0f, pxf = Float64(x0), Float64(y0), Float64(pxw)
    iszero(x0) && iszero(y0) || @warn("save_smite: the camera's first pixel edge is at ($x0, $y0) µm; " *
        "SMITE has no origin, so the file loads back with the camera at 0 and positions moved by that offset")
    s["X"] = _pixel_column("x", [e.x for e in smld.emitters], v -> (v - x0f) / pxf + 0.5)
    s["Y"] = _pixel_column("y", [e.y for e in smld.emitters], v -> (v - y0f) / pxf + 0.5)
    if zcols !== nothing
        s["Z"] = zcols[1]
    end
    
    s["Photons"] = _mat_column("photons", [e.photons for e in smld.emitters])
    s["Bg"] = _mat_column("bg", [e.bg for e in smld.emitters])
    
    s["X_SE"] = _pixel_column("σ_x", [e.σ_x for e in smld.emitters], v -> v / pxf)
    s["Y_SE"] = _pixel_column("σ_y", [e.σ_y for e in smld.emitters], v -> v / pxf)
    if zcols !== nothing
        s["Z_SE"] = zcols[2]
    end
    
    s["Photons_SE"] = _mat_column("σ_photons", [e.σ_photons for e in smld.emitters])
    s["Bg_SE"] = _mat_column("σ_bg", [e.σ_bg for e in smld.emitters])
    
    # Counts and indices are doubles, as most SMITE results files hold them (MATLAB integer arithmetic
    # rounds); one a double cannot hold exactly is an error, not a rounded value
    s["FrameNum"] = _mat_counts("frame", [e.frame for e in smld.emitters])
    s["DatasetNum"] = _mat_counts("dataset", [e.dataset for e in smld.emitters])
    s["ConnectID"] = _mat_counts("track_id", [e.track_id for e in smld.emitters])

    # Add metadata: what the loaders need to read the file back, from the camera
    s["NFrames"] = _mat_count("n_frames", smld.n_frames, _SMITE_COUNT_MIN["NFrames"])
    s["NDatasets"] = _mat_count("n_datasets", smld.n_datasets, _SMITE_COUNT_MIN["NDatasets"])
    s["XSize"] = Float64(nx)
    s["YSize"] = Float64(ny)
    s["PixelSize"] = pxw
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
            # Z and Z_SE (kept by load_smite_2d from a 3D file) are coordinates the loaders read as
            # numbers: written as a number column (a cell of them, as 0.1.0 could write, included),
            # exactly or not at all. Every other field keeps its MATLAB class: SMITE's IndSMD is a
            # cell, one entry per localization, and ThreshFlag is uint32.
            if key in ("Z", "Z_SE")
                v = _mat_floats(value)
                v isa Integer && throw(ArgumentError("save_smite: $key of file row $(rows[v]) is " *
                    "$(repr(value[v])), which is not a real number a MATLAB double holds exactly; " *
                    "nothing was written"))
                value = v
            end
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
    # A 3D file loaded as 2D keeps Z and Z_SE in metadata, written back above for each emitter: say
    # 3D when both are there with one value per localization, as load_smite_3d needs
    zcols === nothing && !isempty(smld.emitters) &&
        all(k -> haskey(s, k) && _smite_length(s[k]) == length(smld.emitters), ("Z", "Z_SE")) &&
        (s["NDims"] = 3.0)
    
    # MAT truncates the target before it checks a struct's keys. Check them all first; an
    # existing file is rewritten only once the new contents are complete.
    _check_mat_keys(s)
    _write_smd(joinpath(String(filepath), String(filename)), s)
    return nothing
end

# The whole file is written to a temporary file first and copied into place only when complete,
# in place so that an existing file keeps its owner, group, permissions and links; the guide's
# SMITE files page states the rules. `copy` is a seam for tests.
function _write_smd(target::String, s; copy=_copy_into)
    dir, base = splitdir(target)
    (isempty(base) || base in (".", "..")) &&
        throw(ArgumentError("save_smite: $(repr(target)) names no file"))
    isempty(dir) && (dir = ".")
    isdir(dir) || throw(ArgumentError(ispath(dir) ? "save_smite: $(repr(dir)) is not a folder" :
                                      "save_smite: the folder $(repr(dir)) does not exist"))
    isdir(target) && throw(ArgumentError("save_smite: $(repr(target)) is a folder"))
    islink(target) && !ispath(target) &&
        throw(ArgumentError("save_smite: $(repr(target)) is a symlink to nothing"))
    ispath(target) && !isfile(target) &&
        throw(ArgumentError("save_smite: $(repr(target)) is not a regular file"))
    existed = isfile(target)
    if existed
        try
            close(open(target, "r+"))
        catch e
            e isa SystemError || rethrow()
            throw(ArgumentError("save_smite: $(repr(target)) cannot be written: " *
                                Libc.strerror(e.errnum)))
        end
    end
    tmp = joinpath(tempdir(), string("smlmio-", getpid(), "-", time_ns(), ".tmp"))
    try
        matwrite(tmp, Dict("SMD" => s))
        if existed
            _check_space(target, filesize(tmp) - filesize(target))
            # HDF5 refuses to recreate a file it holds open, before anything changes
            close(matopen(target, "w"))
        end
        try
            copy(tmp, target)
        catch e
            # A new file that could not be finished is removed
            existed || (try rm(target; force=true) catch end; rethrow())
            e isa InterruptException && rethrow()
            # An existing file was truncated: say so, and that the data can be saved again
            error("save_smite: writing $(repr(target)) failed partway, so that file is now " *
                  "incomplete and its old contents are gone. The data are unchanged in memory: " *
                  "save again once the cause is fixed. The error was: $(sprint(showerror, e))")
        end
    finally
        rm(tmp; force=true)
    end
    return nothing
end

# Copies the file `src` into `dst` in place (truncating it, or creating it).
function _copy_into(src, dst)
    open(src, "r") do from
        open(dst, "w") do to
            buf = Vector{UInt8}(undef, 1 << 20)
            while !eof(from)
                n = readbytes!(from, buf)
                write(to, view(buf, 1:n))
            end
        end
    end
end

# Refuses before an existing file is truncated when its folder lacks the extra bytes the new
# contents need (a quota can still stop the write). Skipped where the system cannot say,
# including a file system that reports no size at all.
function _check_space(target, extra)
    extra > 0 || return nothing
    st = try
        diskstat(dirname(realpath(target)))
    catch
        return nothing
    end
    st.total == 0 && return nothing
    free = st.available
    free >= extra || throw(ArgumentError("save_smite: $(repr(target)) needs $extra more bytes " *
                                         "than its folder has free ($free)"))
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
# An array of Dicts is a MATLAB struct array, whose field names come from the keys. They must be
# strings too: MAT 0.10 refuses Symbol keys there, so every allowed MAT version behaves the same.
function _collect_bad_keys!(bad, path, v::AbstractArray)
    eltype(v) <: Union{Number,AbstractString} && return nothing
    foreach(x -> _collect_bad_keys!(bad, path, x), v)
end
_collect_bad_keys!(bad, path, _) = nothing

# A numeric column as MATLAB writes it: single when every value is Float32, otherwise double, each
# value a real number (a complex one with a zero imaginary part is its real part).
# The column is converted by its values, not its element type: MAT writes a vector whose element
# type is not concrete (Float32 and Float64 emitters together, say) as a cell, which SMITE cannot
# compute with, and refuses some types (Float16, BigFloat). Returns the converted vector, or the
# index of the first value that is not a real number a double holds exactly.
function _mat_floats(col)
    isempty(col) && return eltype(col) === Float32 ? Float32[] : Float64[]
    all(x -> x isa Float32, col) && return convert(Vector{Float32}, col)
    out = Vector{Float64}(undef, length(col))
    for (i, x) in enumerate(col)
        _is_real(x) || return i
        r = real(x)
        d = Float64(r)
        r isa Union{Float32,Float64} || isequal(d, r) || return i
        out[i] = d
    end
    return out
end
# A real number: a Real, or a complex number whose imaginary part is zero
_is_real(x) = x isa Number && iszero(imag(x))

# The emitters' column `what`, or an ArgumentError naming the first emitter whose value cannot be
# written so (a value that is not real, or one a double would round, such as a BigFloat or an
# integer above 2^53), before anything is written.
function _mat_column(what, col)
    r = _mat_floats(col)
    r isa Integer && throw(ArgumentError("save_smite: $what of emitter $r is $(repr(col[r])), which " *
        "is not a real number a MATLAB double holds exactly; nothing was written"))
    return r
end

# A position column in SMITE's pixels: f of each emitter's value, computed in Float64 from the
# value as checked by _mat_column, then single when the emitters' column is all Float32, so Float32
# data stay single. A finite value whose pixel coordinate the written type cannot hold is an
# ArgumentError naming the emitter.
function _pixel_column(what, col, f)
    own = _mat_column(what, col)
    vals = Float64[f(Float64(v)) for v in own]
    out = eltype(own) === Float32 ? Float32.(vals) : vals
    bad = findfirst(i -> isfinite(own[i]) && !isfinite(out[i]), eachindex(own))
    bad === nothing || throw(ArgumentError("save_smite: $what of emitter $bad is $(repr(col[bad])), " *
        "whose value in SMITE's pixels a MATLAB $(eltype(out) === Float32 ? "single" : "double") " *
        "cannot hold; nothing was written"))
    return out
end

# A count or index as a MATLAB double, by the loaders' rule (_fits_int, _SMITE_COUNT_MIN) so that
# every file save_smite writes loads: a whole number in Int's range, at least `lo`, that a double
# holds exactly (every whole number up to 2^53, and some above). Anything else is an
# ArgumentError naming the value, before anything is written.
function _mat_count(what, v, lo=typemin(Int))
    _mat_count_ok(v, lo) || _bad_count(what, v, lo)
    return Float64(real(v))
end
function _mat_counts(what, vals)
    bad = findfirst(v -> !_mat_count_ok(v, typemin(Int)), vals)
    bad === nothing || _bad_count("$what of emitter $bad", vals[bad], typemin(Int))
    return Float64[Float64(real(v)) for v in vals]
end
_mat_count_ok(v, lo) = _fits_int(v) && real(v) >= lo && Float64(real(v)) == real(v)
function _bad_count(what, v, lo)
    bound = lo == typemin(Int) ? "" : " of at least $lo"
    throw(ArgumentError("save_smite: $what is $(repr(v)); SMITE stores it as a MATLAB double, so " *
        "it must be a whole number$bound in Int's range that a double holds exactly (every whole " *
        "number up to 2^53); nothing was written"))
end
