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
  `pixel_edges_y`. `Emitter2D` and `Emitter3D` cannot be saved. Emitter values must be
  `Float32` or `Float64` (metadata `Z` and `Z_SE` written as coordinates may also be complex
  ones); `frame`, `dataset` and `track_id` integers, and `n_frames` and `n_datasets` integers of
  at least 0, each of at most 2^53 in size; camera edges `Float32` or `Float64` values or such
  integers. Anything else is an `ArgumentError` naming it.
- `filepath::AbstractString`: an existing folder.
- `filename::AbstractString`: name of the `.mat` file; an existing file of that name is
  rewritten.

# Returns
- `nothing`.

# Notes
- The file is written in full to a new file in the same folder and then renamed over the
  target, so the name holds the old file or the whole new one, and a save that fails while
  writing leaves an existing file as it was. A rename that fails, or is interrupted, keeps the
  new file and names it; on the lab's network share (CIFS) a file another program has open
  cannot be replaced, and the share then deletes the old file, so the name can hold neither
  file and the kept new file holds the data. The saved file is a new file (owned by whoever
  saves, with the folder's usual mode for a new file). What it refuses, and how a file is
  replaced, are on the [SMITE files](@ref) page.

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
    # z and σ_z. Anything else saves without Z; non-empty 3D or mixed data that lose Z get a
    # warning.
    d = emitter_ndims(smld.emitters)
    n = length(smld.emitters)
    zcols = nothing
    if eltype(smld.emitters) <: Emitter3DFit || (d == 3 && n > 0 &&
           all(e -> hasproperty(e, :z) && hasproperty(e, :σ_z), smld.emitters))
        zcols = (_lab_column("z", [e.z for e in smld.emitters]),
                 _lab_column("σ_z", [e.σ_z for e in smld.emitters]))
    end
    # The emitters' own Z was not saved: metadata Z and Z_SE, below, must not stand in for it
    lostz = zcols === nothing && n > 0 && d !== 2
    if lostz
        d === nothing && @warn("save_smite: mixed 2D/3D emitters: Z not saved; " *
                               "save 2D and 3D data separately")
        d == 3 && @warn("save_smite: 3D emitters without σ_z: Z not saved")
    end

    # Extract arrays from emitters. SMITE stores X, Y and their standard errors in pixels, (1, 1)
    # the center of the top-left pixel; positions here are microns from the camera's corner.
    # SMITE has one PixelSize, so the camera must be the one the loaders will build from XSize,
    # YSize and the written pixel size (single for a Float32 camera, so the file loads back with a
    # Float32 camera; otherwise double): each edge, relative to the camera's first, within a
    # thousandth of a pixel of the rebuilt camera's, compared in Float64 (exact for the edges
    # accepted: single, double, and integers up to 2^53).
    cam = smld.camera
    ex, ey = cam.pixel_edges_x, cam.pixel_edges_y
    nx, ny = length(ex) - 1, length(ey) - 1
    nx >= _SMITE_COUNT_MIN["XSize"] && ny >= _SMITE_COUNT_MIN["YSize"] || throw(ArgumentError(
        "save_smite: the camera has $nx × $ny pixels; SMITE's XSize and YSize must be at least " *
        "$(_SMITE_COUNT_MIN["XSize"]) and $(_SMITE_COUNT_MIN["YSize"])"))
    edge_ok(v) = v isa _LAB_REAL ? isfinite(v) : _lab_whole(v)
    for (axis, edges) in (("x", ex), ("y", ey))
        bad = findfirst(!edge_ok, edges)
        bad === nothing || throw(ArgumentError("save_smite: the camera's pixel edge $bad along $axis " *
            "is $(repr(edges[bad])) ($(typeof(edges[bad]))); SMLMio saves edges that are finite " *
            "Float32 or Float64 values, or integers of at most 2^53 in size; nothing was written"))
    end
    fx, fy = Vector{Float64}(ex), Vector{Float64}(ey)
    x0, y0 = fx[1], fy[1]
    px = fx[2] - fx[1]
    pxw = _smite_pixel_size("save_smite", "the camera's pixel size", eltype(ex) === Float32 ? Float32(px) : px)
    rebuilt = _smite_camera("save_smite", nx, ny, pxw)
    pxf = Float64(pxw)
    same(f, re) = all(k -> abs((f[k] - f[1]) - Float64(re[k])) <= pxf / 1000, eachindex(f))
    same(fx, rebuilt.pixel_edges_x) && same(fy, rebuilt.pixel_edges_y) ||
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
    
    s["Photons"] = _lab_column("photons", [e.photons for e in smld.emitters])
    s["Bg"] = _lab_column("bg", [e.bg for e in smld.emitters])
    
    s["X_SE"] = _pixel_column("σ_x", [e.σ_x for e in smld.emitters], nothing, pxf)
    s["Y_SE"] = _pixel_column("σ_y", [e.σ_y for e in smld.emitters], nothing, pxf)
    if zcols !== nothing
        s["Z_SE"] = zcols[2]
    end
    
    s["Photons_SE"] = _lab_column("σ_photons", [e.σ_photons for e in smld.emitters])
    s["Bg_SE"] = _lab_column("σ_bg", [e.σ_bg for e in smld.emitters])

    # Counts and indices are doubles, as most SMITE results files hold them (MATLAB integer arithmetic
    # rounds): integers a double holds exactly
    s["FrameNum"] = _lab_counts("frame", [e.frame for e in smld.emitters])
    s["DatasetNum"] = _lab_counts("dataset", [e.dataset for e in smld.emitters])
    s["ConnectID"] = _lab_counts("track_id", [e.track_id for e in smld.emitters])

    # Add metadata: what the loaders need to read the file back, from the camera
    s["NFrames"] = _lab_count("n_frames", smld.n_frames, "NFrames")
    s["NDatasets"] = _lab_count("n_datasets", smld.n_datasets, "NDatasets")
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
    # so. On 2D emitters they are written through the same rule as the emitters' z, so a save and
    # a later save of its load agree, except that complex values (SMITE's fits give them; a 2D load
    # does not screen Z) are written as complex double. The file says NDims 3 when both are written
    # and it has localizations, the loaded file's Z had values (SMITE writes NDims 2 for its 3D
    # fits too) or the metadata's NDims is 3. When the emitters' own Z was not saved (lostz),
    # metadata Z and Z_SE are not saved at all.
    loaded = haskey(md, "smite_per_localization_fields")
    zkeys, zsrc, zunit, zids, zhad = String[], Vector[], "emitter", collect(1:n), n > 0
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
    for (k, src) in zip(zkeys, zsrc)
        mdz[k] = _lab_column("metadata $k", src; complex=true, unit=zunit, ids=zids,
                             hint=" (remove it from the metadata to save without it)")
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
    
    # Every key is checked first, so one error names them all
    _check_mat_keys(s)
    _write_smd(joinpath(String(filepath), String(filename)), s)
    return nothing
end

# The file is written in full to a new file beside the target and then renamed over it, so the
# name holds the old file or the whole new one, never a partial one, and a program that has the
# old file open goes on reading the old contents (local disks and NFS; the CIFS share refuses the
# rename, below). The new file is created as any new file in
# that folder is (owned by whoever saves, with the mode the folder's defaults give); a symbolic
# link keeps pointing at the saved file. The guide's SMITE files page states the rules.
# `rename` is a seam for tests.
function _write_smd(target::String, s; rename=_rename)
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
    path = islink(target) ? realpath(target) : target
    # An existing file you may not write is refused, as a write would be: the rename needs only
    # the folder's permission, and would replace a file someone made read-only
    if isfile(path)
        try
            close(open(path, "r+"))
        catch e
            e isa SystemError || rethrow()
            throw(ArgumentError("save_smite: $(repr(target)) cannot be written: " *
                                Libc.strerror(e.errnum)))
        end
    end
    tmp, mode = _new_beside(path)
    try
        matwrite(tmp, Dict("SMD" => s))
        mode === nothing || chmod(tmp, mode)
    catch
        try
            rm(tmp; force=true)
        catch   # a lost share: the write's own error is the one to see
        end
        rethrow()
    end
    # A rename that fails keeps the new file, complete, and names it: on the lab's network share
    # (CIFS), a refused replace (another program has the file open, for one) makes the share delete
    # the old file, at once or when that program closes it, so the new file may then be the only
    # copy of the data. So does an interrupt. A new file no longer beside the target means the
    # rename may have completed (NFS can report an error for a rename it carried out); one that
    # cannot be looked for (the folder cannot be read) leaves that unknown. Either way the target is
    # named to check, and an interrupt is passed on after a warning that says the same.
    try
        rename(tmp, path)
    catch e
        kept = try
            ispath(tmp)
        catch
            nothing
        end
        msg = if kept === true
            "the new file could not be renamed over $(repr(target)): $(sprint(showerror, e)). The " *
            "data are saved in full in $(repr(tmp)); rename it to $(repr(basename(path))) once the " *
            "cause is fixed. On the lab's network share (CIFS) a file another program has open " *
            "cannot be replaced, and the share then deletes the old file."
        else
            "renaming the new file $(repr(tmp)) over $(repr(target)) reported " *
            "$(sprint(showerror, e)), and " * (kept === false ?
                "the new file is no longer beside it, so the save may have completed" :
                "the folder cannot be read, so whether the save completed is not known") *
            "; check $(repr(target))."
        end
        e isa InterruptException && (@warn("save_smite: interrupted; $msg"); rethrow())
        error("save_smite: $msg")
    end
    return nothing
end

# A rename as the one system call. Base.Filesystem.rename on Julia 1.10 and 1.11 copies the file
# and removes the source when the call fails, which could put a partial file at the name.
function _rename(src, dst)
    err = ccall(:jl_fs_rename, Int32, (Cstring, Cstring), src, dst)
    err < 0 && Base.uv_error("rename($(repr(src)), $(repr(dst)))", err)
    return nothing
end

# A new, empty file in `path`'s folder under a name no other file has, created exclusively (so
# no other save can take the same name) with the mode any new file there gets: the system applies
# the user's umask and the folder's default permissions. It is closed. Returns its name, and the
# mode to give it back once written when that mode does not let its owner read and write it (a
# folder whose default permissions make new files read-only): it is made so meanwhile.
function _new_beside(path)
    dir = dirname(path)
    isempty(dir) && (dir = ".")
    F = Base.Filesystem
    for _ in 1:100
        tmp = tempname(dir; cleanup=false)
        f = try
            F.open(tmp, F.JL_O_WRONLY | F.JL_O_CREAT | F.JL_O_EXCL, 0o666)
        catch e
            e isa Base.IOError && e.code == Base.UV_EEXIST && continue
            e isa Base.IOError || rethrow()
            throw(ArgumentError("save_smite: cannot create a file in $(repr(dir)) to save into: " *
                                sprint(showerror, e)))
        end
        try
            close(f)
            m = filemode(tmp) & 0o7777
            m & 0o600 == 0o600 && return tmp, nothing
            chmod(tmp, m | 0o600)
            return tmp, m
        catch
            try
                rm(tmp; force=true)
            catch   # a lost share: the first error is the one to see
            end
            rethrow()
        end
    end
    error("save_smite: found no unused file name in $(repr(dir))")
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

# The metadata as save_smite reads it, a copy that takes any value (metadata of a narrower value
# type, such as Dict{String,String}, included): every value made writable (_mat_writable), and
# empty text in SMITE's own fields made an empty array, as the loaders do (_smite_empty_text!).
function _smite_metadata(md)
    out = Dict{keytype(md),Any}(md)
    for (k, v) in md
        out[k] = _mat_writable(v)
    end
    return _smite_empty_text!(out)
end

# What save_smite writes: the values SMITE files and the lab's packages hold. A measured value (a
# position, an uncertainty, photons, background, z) is a single or double, Float32 or Float64; Z
# and Z_SE kept in metadata from a file may also be complex, as SMITE's fits can give them. A
# count or index (frame, dataset, track id, n_frames, n_datasets), or a camera edge given as an
# integer, is an integer a double holds exactly: at most 2^53 in size. Anything else is one
# ArgumentError naming the field, before anything is written.
const _LAB_REAL = Union{Float32,Float64}
const _LAB_NUMBER = Union{Float32,Float64,ComplexF32,ComplexF64}

# The column `what` as MATLAB writes it: single when every value is Float32, otherwise double
# (complex double when a value is complex and `complex` allows it). MAT writes a vector whose
# element type is not concrete as a cell, which SMITE cannot compute with, hence the conversion,
# exact for these types. `unit` and `ids` say where a refused value is: an emitter and its index,
# or a file row and its row.
function _lab_column(what, col; complex=false, unit="emitter", ids=eachindex(col), hint="")
    ok = complex ? _LAB_NUMBER : _LAB_REAL
    bad = findfirst(v -> !(v isa ok), col)
    bad === nothing || throw(ArgumentError("save_smite: $what of $unit $(ids[bad]) is " *
        "$(repr(col[bad])) ($(typeof(col[bad]))); SMLMio saves " *
        (complex ? "single or double values, real or complex" : "Float32 or Float64 values") *
        "$hint; nothing was written"))
    isempty(col) && return eltype(col) === Float32 ? Float32[] : Float64[]
    all(v -> v isa Float32, col) && return Vector{Float32}(col)
    all(v -> v isa _LAB_REAL, col) && return Vector{Float64}(col)
    return Vector{ComplexF64}(col)
end

# A column in SMITE's pixels from the emitters' values: a position (`origin`, the camera's first
# edge) X = (x - origin) / px + 0.5, or an uncertainty (`origin` nothing) X_SE = σ / px, computed
# in Float64, then single when the emitters' column is all Float32. A finite value that the
# loaders would read back as infinite, or a non-zero uncertainty they would read back as 0 (values
# far outside any microscope's range), is an ArgumentError naming the emitter.
function _pixel_column(what, col, origin, px)
    own = _lab_column(what, col)
    x = Vector{Float64}(own)
    vals = origin === nothing ? x ./ px : (x .- origin) ./ px .+ 0.5
    out = eltype(own) === Float32 ? Vector{Float32}(vals) : vals
    back(X) = origin === nothing ? Float64(X) * px : (Float64(X) - 0.5) * px   # as the loaders compute
    ok(i) = !isfinite(own[i]) ||
        (isfinite(back(out[i])) && (origin !== nothing || iszero(own[i]) || !iszero(back(out[i]))))
    bad = findfirst(!ok, eachindex(own))
    bad === nothing || throw(ArgumentError("save_smite: $what of emitter $bad is $(own[bad]), which " *
        "a MATLAB $(eltype(out) === Float32 ? "single" : "double") in SMITE's pixels cannot hold so " *
        "that it loads back; nothing was written"))
    return out
end

# An integer a MATLAB double holds exactly: at most 2^53 in size, and at least `lo`.
_lab_whole(v, lo=-2^53) = v isa Integer && max(lo, -2^53) <= v <= 2^53

# A count as a MATLAB double, at least the loaders' bound for `key` (_SMITE_COUNT_MIN), so that
# every file save_smite writes loads; and the emitters' indices `what`, one per emitter.
function _lab_count(what, v, key)
    lo = _SMITE_COUNT_MIN[key]
    _lab_whole(v, lo) || _not_whole(what, v, lo)
    return Float64(v)
end
function _lab_counts(what, vals)
    bad = findfirst(v -> !_lab_whole(v), vals)
    bad === nothing || _not_whole("$what of emitter $bad", vals[bad], nothing)
    return Float64[Float64(v) for v in vals]
end
_not_whole(what, v, lo) = throw(ArgumentError("save_smite: $what is $(repr(v)) ($(typeof(v))); " *
    "SMITE stores it as a MATLAB double, so SMLMio saves an integer" *
    (lo === nothing ? "" : " of at least $lo") * " of at most 2^53 in size; nothing was written"))
