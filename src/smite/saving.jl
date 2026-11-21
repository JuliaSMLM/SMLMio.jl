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
  fails while writing leaves an existing file as it was; only an error once the file is being
  rewritten (a full disk, a file size limit, a lost share) can leave it partly written, and the
  error then says so. It does not coordinate with other programs using the same file. What it
  refuses, and how a file is rewritten, are on the [SMITE files](@ref) page.

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
    # z and σ_z, written by _z_columns. Anything else saves without Z; non-empty 3D or mixed data
    # that lose Z get a warning.
    d = emitter_ndims(smld.emitters)
    n = length(smld.emitters)
    zcols = nothing
    if eltype(smld.emitters) <: Emitter3DFit || (d == 3 && n > 0 &&
           all(e -> hasproperty(e, :z) && hasproperty(e, :σ_z), smld.emitters))
        zcols = _z_columns(("z", "σ_z"), ([e.z for e in smld.emitters], [e.σ_z for e in smld.emitters]),
                           "emitter", 1:n)
    end
    # The emitters' own Z was not saved: metadata Z and Z_SE, below, must not stand in for it
    lostz = zcols === nothing && n > 0 && d !== 2
    if lostz
        d === nothing && @warn("save_smite: mixed 2D/3D emitters: Z not saved; " *
                               "save 2D and 3D data separately")
        d == 3 && @warn("save_smite: 3D emitters whose z or σ_z are not all numbers: Z not saved")
    end
    
    # Extract arrays from emitters. SMITE stores X, Y and their standard errors in pixels, (1, 1)
    # the center of the top-left pixel; positions here are microns from the camera's corner.
    # SMITE has one PixelSize, so the camera must be the one the loaders will build from XSize,
    # YSize and the written pixel size (single for a Float32 camera, so the file loads back with a
    # Float32 camera; otherwise double): each edge, relative to the camera's first, within a
    # thousandth of a pixel of the rebuilt camera's. Differences of edges are taken exactly and
    # rounded once (_exact_sub: no overflow in a small integer type, no rounding of a Float16 or
    # rational edge before the subtraction, whatever BigFloat's working precision) and compared in
    # Float64.
    cam = smld.camera
    ex, ey = cam.pixel_edges_x, cam.pixel_edges_y
    nx, ny = length(ex) - 1, length(ey) - 1
    nx >= _SMITE_COUNT_MIN["XSize"] && ny >= _SMITE_COUNT_MIN["YSize"] || throw(ArgumentError(
        "save_smite: the camera has $nx × $ny pixels; SMITE's XSize and YSize must be at least " *
        "$(_SMITE_COUNT_MIN["XSize"]) and $(_SMITE_COUNT_MIN["YSize"])"))
    all(isfinite, ex) && all(isfinite, ey) ||
        throw(ArgumentError("save_smite: the camera's pixel edges are not all finite"))
    x0, y0 = ex[1], ey[1]
    px = _exact_sub(ex[2], ex[1])
    pxw = _smite_pixel_size("save_smite", "the camera's pixel size", eltype(ex) === Float32 ? Float32(px) : px)
    rebuilt = _smite_camera("save_smite", nx, ny, pxw)
    pxf = Float64(pxw)
    same(edges, e0, re) = all(k -> abs(_exact_sub(edges[k], e0) - Float64(re[k])) <= pxf / 1000,
                              eachindex(edges))
    same(ex, x0, rebuilt.pixel_edges_x) && same(ey, y0, rebuilt.pixel_edges_y) ||
        throw(ArgumentError("save_smite: SMITE stores one square pixel size; this camera's pixels are " *
                            "not uniform and square (its edges are not those of $nx × $ny pixels of " *
                            "$pxw µm)"))
    # Positions in SMITE's pixels are computed in Float64 from the written pixel size, as the
    # loaders compute them back
    iszero(x0) && iszero(y0) || @warn("save_smite: the camera's first pixel edge is at ($x0, $y0) µm; " *
        "SMITE has no origin, so the file loads back with the camera at 0 and positions moved by that offset")
    s["X"] = _pixel_column("x", [e.x for e in smld.emitters], x0, pxf)
    s["Y"] = _pixel_column("y", [e.y for e in smld.emitters], y0, pxf)
    if zcols !== nothing
        s["Z"] = zcols[1]
    end
    
    s["Photons"] = _mat_column("photons", [e.photons for e in smld.emitters])
    s["Bg"] = _mat_column("bg", [e.bg for e in smld.emitters])
    
    s["X_SE"] = _pixel_column("σ_x", [e.σ_x for e in smld.emitters], nothing, pxf)
    s["Y_SE"] = _pixel_column("σ_y", [e.σ_y for e in smld.emitters], nothing, pxf)
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
    s["NFrames"] = _mat_count("n_frames", smld.n_frames, "NFrames")
    s["NDatasets"] = _mat_count("n_datasets", smld.n_datasets, "NDatasets")
    s["XSize"] = Float64(nx)
    s["YSize"] = Float64(ny)
    s["PixelSize"] = pxw
    s["NDims"] = zcols === nothing ? 2.0 : 3.0

    md = _smite_metadata(smld.metadata)

    # Optional fields the file left empty, which the loader filled with defaults and that still hold
    # only the default, go back empty (SMITE treats an empty ConnectID differently from zeros). A
    # field the file held with one value per dataset is written as the default column instead: its
    # raw values are not one per localization, and a later load could take them for that.
    defaulted = get(md, "smite_defaulted_fields", String[])
    for (key, default) in _SMITE_OPTIONAL
        key in defaulted && haskey(md, key) && _smite_empty(md[key]) &&
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

    # Z and Z_SE in the metadata with one value per saved emitter: for data loaded from SMITE, the
    # listed per-localization fields, by row (load_smite_2d keeps them from a 3D file); for other
    # data (built by hand, or loaded by SMLMData 0.7.1, which kept a 3D file's Z that way), fields
    # with as many values as emitters, in the order given, as 0.1.0 wrote them. Such data carry no
    # record of which row each value belongs to (ids are not that: a join, a filter or frame
    # connection keeps or reuses them), so that order cannot be checked; the warning below says
    # so. On 2D emitters they are written through the same writer as the emitters' z, so a save
    # and a later save of its load agree; a value that is not a number drops them with a warning. The file says NDims 3 when both are written
    # and it has localizations, the loaded file's Z had values (SMITE writes NDims 2 for its 3D
    # fits too) or the metadata's NDims is 3. When the emitters' own Z was not saved (lostz),
    # metadata Z and Z_SE are not saved at all.
    loaded = haskey(md, "smite_per_localization_fields")
    zkeys, zsrc, zunit, zids, zhad = String[], Vector[], "metadata row", collect(1:n), n > 0
    if zcols === nothing && !lostz
        if loaded
            zkeys = filter(k -> k in towrite, ["Z", "Z_SE"])
            zsrc = [_smite_vec(md[k])[rows] for k in zkeys]
            zunit, zids, zhad = "file row", rows, any(k -> _smite_length(md[k]) > 0, zkeys)
        else
            zkeys = filter(k -> haskey(md, k) && _smite_length(md[k]) == n, ["Z", "Z_SE"])
            zsrc = [_smite_vec(md[k]) for k in zkeys]
        end
    end
    mdz = Dict{String,Any}()
    if !isempty(zkeys)
        conv = _z_columns(zkeys, zsrc, zunit, zids)
        if conv === nothing
            @warn("save_smite: $(join(zkeys, ", ")) not saved: they hold values that are not numbers")
        else
            foreach((k, c) -> mdz[k] = c, zkeys, conv)
        end
    end
    ndims3_from_md = length(mdz) == 2 && (n > 0 || zhad || _is_ndims3(get(md, "NDims", nothing)))
    # Data not loaded by SMLMio: SMITE's per-localization fields written from the metadata with one
    # value per emitter are in the order given, which save_smite cannot check. Fields written from
    # the emitters (s holds them by now) are not among them; Z and Z_SE count when written (mdz).
    unchecked = loaded || n == 0 ? String[] :
        sort!([k for (k, v) in md if k in _SMITE_STACKABLE && _smite_length(v) == n &&
               (k in ("Z", "Z_SE") ? haskey(mdz, k) : !haskey(s, k))])
    isempty(unchecked) || @warn("save_smite: $(join(unchecked, ", ")) written in the order given, one " *
        "value per emitter: these data were not loaded by SMLMio, so nothing records which " *
        "localization each value belongs to (a sort, filter, join or frame connection since the " *
        "values were made would misplace them)")

    # Add any additional fields from metadata, except the loaders' own entries. For data loaded
    # from SMITE, a field SMITE stacks per localization that is not listed as one value per file
    # row (the file held it otherwise, often one value per dataset, or it was added after loading)
    # is dropped when the saved localizations are as many as its values, since a load would take
    # it for one value per localization; with another count it is written as it is.
    ambiguous = String[]
    for (key, value) in md
        (key in _SMITE_LOADER_KEYS || haskey(s, key)) && continue
        if key in ("Z", "Z_SE") && (lostz || key in zkeys)
            haskey(mdz, key) && (s[key] = mdz[key])
            continue
        end
        if key in perloc
            key in towrite || continue
            # Every other field keeps its MATLAB class: SMITE's IndSMD is a cell, one entry per
            # localization, and ThreshFlag is uint32
            value = _smite_vec(value)[rows]
        elseif loaded && key in _SMITE_STACKABLE && n > 0 && _smite_length(value) == n
            push!(ambiguous, key)
            continue
        end
        s[key] = value
    end
    isempty(ambiguous) || @warn("save_smite: $(join(sort!(ambiguous), ", ")) not saved: SMITE stores these " *
        "per localization, they are not listed in smite_per_localization_fields as one value per file " *
        "row, and with as many values as saved localizations a load would take them for one value per " *
        "localization; list a field there, with one value per file row, to save it by row")
    ndims3_from_md && (s["NDims"] = 3.0)
    
    # MAT truncates the target before it checks a struct's keys. Check them all first; an
    # existing file is rewritten only once the new contents are complete.
    _check_mat_keys(s)
    _write_smd(joinpath(String(filepath), String(filename)), s)
    return nothing
end

# The whole file is written to a temporary file first (in a folder of its own, _staging_dir) and
# copied into place only when complete,
# in place so that an existing file keeps its owner, group, permissions and links; the guide's
# SMITE files page states the rules. `copy` and `recreate` are seams for tests.
function _write_smd(target::String, s; copy=_copy_into, recreate=_recreate)
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
    stage = _staging_dir()
    tmp = joinpath(stage, "SMD.mat")
    try
        matwrite(tmp, Dict("SMD" => s))
        # A new file is created here, exclusively, so that a failure removes only a file this save
        # made; one that appeared at the target while the save was staged is an existing file
        created = !existed && _create_new(target)
        if !created
            _check_space(target, filesize(tmp) - filesize(target))
            # HDF5 refuses to recreate a file this session holds open, before it truncates it. With
            # file locking on, it truncates a file another program holds locked and then refuses
            # it: a refusal that changed the file is reported as an incomplete file
            st0 = stat(target)
            try
                recreate(target)
            catch e
                _unchanged(st0, stat(target)) && rethrow()
                _incomplete(target, e)
            end
        end
        try
            copy(tmp, target)
        catch e
            # A new file that could not be finished is removed
            created && (try rm(target; force=true) catch end; rethrow())
            _incomplete(target, e)
        end
    finally
        rm(stage; recursive=true, force=true)
    end
    return nothing
end

_recreate(target) = close(matopen(target, "w"))
# The file's size and times are as they were
_unchanged(a, b) = a.size == b.size && a.mtime == b.mtime && a.ctime == b.ctime

# Creates `path` empty, only if nothing is there (O_EXCL): true when it did, false when a file is
# already there; any other failure is an ArgumentError naming the file.
function _create_new(path)
    F = Base.Filesystem
    f = try
        F.open(path, F.JL_O_WRONLY | F.JL_O_CREAT | F.JL_O_EXCL, 0o666)
    catch e
        e isa Base.IOError && e.code == Base.UV_EEXIST && return false
        e isa Base.IOError || rethrow()
        throw(ArgumentError("save_smite: $(repr(path)) cannot be created: $(sprint(showerror, e))"))
    end
    close(f)
    return true
end

# An existing file was truncated and not rewritten in full: say so, and that the data can be saved
# again (an interrupt in a warning, and passed on as it is)
function _incomplete(target, e)
    msg = "save_smite: writing $(repr(target)) failed partway, so that file is now incomplete and " *
          "its old contents are gone. The data are unchanged in memory: save again once the cause " *
          "is fixed."
    e isa InterruptException && (@warn("$msg The save was interrupted."); rethrow())
    error("$msg The error was: $(sprint(showerror, e))")
end

# A folder of this save's own in the system's temporary folder, created with a name no other save
# uses, of this process or another: saves on several threads, or processes that share the
# folder and a process id (containers), each get their own (mkdtemp makes a new folder, never
# an existing one).
_staging_dir() = mktempdir(tempdir(); prefix=string("smlmio-", getpid(), "-"), cleanup=false)

# Copies the file `src` into `dst` in place (truncating it, or creating it). A write can take fewer
# bytes than it was given without an error (a file size limit): each write must take them all, and
# the file must end up as long as `src` (a short write when the last bytes are flushed included).
function _copy_into(src, dst)
    open(src, "r") do from
        open(dst, "w") do to
            _copy_stream(from, to)
            flush(to)
        end
    end
    a, b = filesize(dst), filesize(src)
    a == b || error("$(repr(dst)) holds $a of the $b bytes written")
    return nothing
end
function _copy_stream(from::IO, to::IO)
    buf = Vector{UInt8}(undef, 1 << 20)
    while !eof(from)
        n = readbytes!(from, buf)
        w = write(to, view(buf, 1:n))
        w == n || error("a write took $w of $n bytes (a full disk or a file size limit)")
    end
    return nothing
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

# A numeric column as MATLAB writes it, converted by its values: single when every value is
# Float32; otherwise double when every value is a real number (a complex one with a zero imaginary
# part is its real part), or complex double when some value has a non-zero imaginary part (SMITE's
# fits can give such values, and the loaders screen them out). MAT writes a vector whose element
# type is not concrete as a cell, which SMITE cannot compute with, and refuses some types (Float16,
# BigFloat), hence the conversion. Every value must convert exactly: the result is the converted
# vector, or (i, kind) for a value that cannot be written, kind :inexact for a number a double would
# round (found wherever it is in the column) and otherwise :notnumber.
function _mat_numbers(col)
    isempty(col) && return eltype(col) === Float32 ? Float32[] : Float64[]
    all(x -> x isa Float32, col) && return convert(Vector{Float32}, col)
    bad = findfirst(x -> x isa Number && !_exact_double(x), col)
    bad === nothing || return (bad, :inexact)
    bad = findfirst(x -> !(x isa Number), col)
    bad === nothing || return (bad, :notnumber)
    all(_is_real, col) && return Float64[_float64(real(x)) for x in col]
    return ComplexF64[ComplexF64(_float64(real(x)), _float64(imag(x))) for x in col]
end
_exact_double(x) = x isa Union{Float32,Float64,ComplexF32,ComplexF64} ||
    (_exact_part(real(x)) && _exact_part(imag(x)))
_exact_part(v) = (f = _float64(v); isnan(f) ? isnan(v) : _same_value(f, v) && signbit(f) == signbit(v))
# A real number: a Real, or a complex number whose imaginary part is zero
_is_real(x) = x isa Number && iszero(imag(x))

# The emitters' column `what`, or an ArgumentError naming the first emitter whose value cannot be
# written, before anything is written.
function _mat_column(what, col)
    r = _mat_numbers(col)
    r isa Tuple && _unwritable(what, "emitter", r[1], col[r[1]])
    return r
end
_unwritable(what, unit, i, v) = throw(ArgumentError("save_smite: $what of $unit $i is $(repr(v)), " *
    "which is not a number that a MATLAB double holds exactly (a complex double for a complex " *
    "value); nothing was written"))

# Z and Z_SE as written, for the emitters (unit "emitter") or for file rows kept in metadata
# ("file row", `ids` the rows): each column through _mat_numbers. A number a double would round is
# an ArgumentError naming the field and the emitter or row; a value that is not a number
# (missing) gives nothing, and the file is saved without these columns.
function _z_columns(names, cols, unit, ids)
    conv = map(_mat_numbers, cols)
    for (name, vals, r) in zip(names, cols, conv)
        r isa Tuple && r[2] === :inexact && _unwritable(name, unit, ids[r[1]], vals[r[1]])
    end
    return any(r -> r isa Tuple, conv) ? nothing : collect(conv)
end

# A column in SMITE's pixels: a position (`origin`, the camera's first edge) X = (x - origin) / px
# + 0.5, or an uncertainty (`origin` nothing) X_SE = σ / px, each from the value as checked by
# _mat_column, in Float64 (ComplexF64 for complex values) and then single when the emitters'
# column is all Float32. The origin, of whatever number type, is subtracted from each value
# exactly and rounded once (_exact_sub). A finite value whose pixel value the loaders would read
# back as infinite, or an uncertainty that would read back as 0 when it is not, is an
# ArgumentError naming the emitter.
function _pixel_column(what, col, origin, px)
    own = _mat_column(what, col)
    T = eltype(own) <: Complex ? Complex{Float64} : Float64
    vals = if origin === nothing
        T[T(v) / px for v in own]
    else
        T[_exact_sub(v, origin) / px + 0.5 for v in own]
    end
    out = eltype(own) === Float32 ? Float32.(vals) : vals
    back(X) = origin === nothing ? X * px : (X - 0.5) * px      # as the loaders compute (px a Float64)
    ok(i) = !isfinite(own[i]) || (isfinite(back(out[i])) && (origin !== nothing || iszero(own[i]) || !iszero(out[i])))
    bad = findfirst(!ok, eachindex(own))
    bad === nothing || throw(ArgumentError("save_smite: $what of emitter $bad is $(repr(col[bad])), " *
        "which a MATLAB $(eltype(out) === Float32 ? "single" : "double") in SMITE's pixels cannot " *
        "hold so that it loads back; nothing was written"))
    return out
end

# A value as MAT can write it, changed only where it must be. MAT reads an unfilled element of a
# MATLAB cell as an empty array of element type Union{}, which its writer refuses; such an element
# is written as an empty double array, as SMITE holds an empty field. A tuple (metadata built in
# Julia) is written as the array it holds; MAT 0.10 refuses tuples. Arrays (cells, or struct
# arrays of Dicts) and plain Dicts (MATLAB structs) are searched; a changed one is a copy whose
# element type is widened only as far as the new values need, so a cell of structs stays a cell
# and an array of Dicts an array of Dicts, and one that needs no change is returned as it is.
# MAT's own types (a class object, a function handle, a struct array) are written as they are,
# keeping their MATLAB class; MAT cannot write one holding an unfilled cell element.
_mat_writable(v::AbstractArray{Union{}}) = zeros(size(v))
function _mat_writable(v::AbstractArray)
    eltype(v) <: Union{Number,AbstractString} && return v
    out = v
    for i in eachindex(v)
        isassigned(v, i) || continue
        x = v[i]
        y = _mat_writable(x)
        y === x && continue
        out === v && (out = copy(v))
        # widened only as far as the new value needs: an array of Dicts (a MATLAB struct array)
        # stays one, a cell stays a cell
        y isa eltype(out) || (out = Array{typejoin(eltype(out), typeof(y))}(out))
        out[i] = y
    end
    return out
end
_mat_writable(v::Tuple) = _mat_writable(collect(v))
function _mat_writable(v::Dict)
    out = v
    for (k, x) in v
        y = _mat_writable(x)
        y === x && continue
        out === v && (out = copy(v))
        y isa valtype(out) || (out = Dict{keytype(v),Any}(out))
        out[k] = y
    end
    return out
end
_mat_writable(v) = v

# The metadata as save_smite reads it, a copy: every value made writable (_mat_writable), and empty
# text in SMITE's own fields made an empty array, as the loaders do (_smite_empty_text!).
function _smite_metadata(md)
    out = copy(md)
    for (k, v) in md
        out[k] = _mat_writable(v)
    end
    return _smite_empty_text!(out)
end

# a - b rounded once to the nearest Float64 (a ComplexF64 for a complex a), without touching
# BigFloat's working precision or rounding mode (which Julia 1.10 shares between threads). When a
# and b both convert to Float64 exactly (every Float16, Float32 and Float64; most integers), one
# IEEE subtraction rounds the exact difference once. Otherwise finite integers, rationals and
# binary floats (BigFloat) are subtracted exactly as fractions of BigInts and rounded once
# (_round_float64). Other numbers and values that are not finite are converted, then subtracted.
function _exact_sub(a::Real, b::Real)
    fa, fb = _float64(a), _float64(b)
    _same_value(fa, a) && _same_value(fb, b) && return fa - fb
    exact = Union{Integer,Rational,Float16,Float32,Float64,BigFloat}
    (a isa exact && b isa exact && isfinite(a) && isfinite(b)) || return fa - fb
    (na, da), (nb, db) = _fraction(a), _fraction(b)
    return _round_float64(na * db - nb * da, da * db)
end
_exact_sub(a::Complex, b::Real) = complex(_exact_sub(real(a), b), _float64(imag(a)))

# x as (numerator, denominator), BigInts with the denominator positive, exactly
_fraction(x::Union{Integer,Rational}) = (BigInt(numerator(x)), BigInt(denominator(x)))
# A finite binary float through Base.decompose (internal, but exact: Rational{BigInt} of a
# BigFloat rounds at BigFloat's working precision): x == s * 2^e / d with d = ±1, the sign in d
function _fraction(x::AbstractFloat)
    s, e, d = Base.decompose(x)
    n = BigInt(s) * sign(d)
    return e >= 0 ? (n << e, big(1)) : (n, big(1) << -e)
end

# The Float64 nearest num / den (den > 0), ties to even, as IEEE rounds: subnormal results keep
# fewer bits, and a value beyond the largest double is infinite. Integer arithmetic only.
function _round_float64(num::BigInt, den::BigInt)
    iszero(num) && return 0.0
    m = abs(num)
    # A quotient of at least 56 bits, the remainder folded into its last bit (sticky), so that it
    # rounds as the exact value does: value = q * 2^e, up to the sticky bit
    k = max(0, 56 + ndigits(den; base=2) - ndigits(m; base=2))
    q, r = divrem(m << k, den)
    q = 2q + !iszero(r)
    e = -k - 1
    nq = ndigits(q; base=2)
    E = nq - 1 + e                       # 2^E <= value < 2^(E + 1)
    p = min(53, E + 1075)                # bits a double keeps at that exponent (2^-1074 the least)
    p < 0 && return copysign(0.0, num)   # below half the least subnormal
    shift = nq - p                       # at least 3
    t = q >> shift
    rest = q - (t << shift)
    half = big(1) << (shift - 1)
    (rest > half || (rest == half && isodd(t))) && (t += 1)
    v = ldexp(Float64(t), e + shift)     # exact: t has at most p bits, or is 2^p
    return num < 0 ? -v : v
end

# x as a Float64, rounded to nearest once: a Rational through _round_float64 (Julia converts a
# Rational{BigInt} through BigFloat at its working precision), anything else by Float64.
_float64(x::Rational) = isfinite(x) ? _round_float64(_fraction(x)...) : Float64(x)
_float64(x::Real) = Float64(x)

# Whether the Float64 f equals x exactly. Julia 1.10 compares a float with a Rational{BigInt}
# through BigFloat at its working precision, so a finite Rational is compared as a fraction.
_same_value(f::Float64, x::Rational) = isfinite(f) && isfinite(x) ?
    ((nf, df) = _fraction(f); nf * denominator(x) == BigInt(numerator(x)) * df) : f == x
_same_value(f::Float64, x::Real) = f == x

# A count or index as a MATLAB double, by the loaders' rule (_fits_int, _SMITE_COUNT_MIN) so that
# every file save_smite writes loads: a whole number in Int's range, at least `lo`, that a double
# holds exactly (every whole number up to 2^53, and some above). Anything else is an
# ArgumentError naming the value, before anything is written.
function _mat_count(what, v, key)
    lo = _SMITE_COUNT_MIN[key]
    _mat_count_ok(v, lo) || _bad_count(what, v, lo)
    return _float64(real(v))
end
function _mat_counts(what, vals)
    bad = findfirst(v -> !_mat_count_ok(v, typemin(Int)), vals)
    bad === nothing || _bad_count("$what of emitter $bad", vals[bad], typemin(Int))
    return Float64[_float64(real(v)) for v in vals]
end
_mat_count_ok(v, lo) = _fits_int(v) && real(v) >= lo && _same_value(_float64(real(v)), real(v))
function _bad_count(what, v, lo)
    bound = lo == typemin(Int) ? "" : " of at least $lo"
    throw(ArgumentError("save_smite: $what is $(repr(v)); SMITE stores it as a MATLAB double, so " *
        "it must be a whole number$bound in Int's range that a double holds exactly (every whole " *
        "number up to 2^53); nothing was written"))
end
