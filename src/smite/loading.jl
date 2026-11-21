# For each of `fields` in `cols` (vectors with one value per localization) that holds a non-zero
# imaginary part, the rows where it does.
function _complex_rows(cols, fields)
    rows = Dict{String,Vector{Int}}()
    for f in fields
        v = cols[f]
        eltype(v) <: Real && continue      # a cell's values may be complex in a Vector{Number}
        r = findall(x -> !iszero(imag(x)), v)
        isempty(r) || (rows[f] = r)
    end
    return rows
end

# MAT reads a 1x1 field as a scalar and an empty field as a 0x0 matrix.
_smite_length(v::AbstractArray) = length(v)
_smite_length(::Number) = 1
_smite_length(_) = -1
_smite_vec(v::AbstractArray) = vec(v)
_smite_vec(v::Number) = [v]

# A per-localization field as a vector of numbers, or nothing when it is not numeric. MAT reads a
# MATLAB cell of numbers (single and double values together, say) as an array of Any; its values
# are kept as stored, in a Vector{Number}, and each is checked and converted on its own where it
# is used: promoting them to one type could round them (an Int64 with a Float32 gives Float32).
function _smite_numeric(v)
    v isa Number && return [v]
    v isa AbstractArray || return nothing
    v isa AbstractArray{<:Number} && return vec(v)
    all(x -> x isa Number, v) || return nothing
    return Number[x for x in vec(v)]
end

# A field SMITE left empty: a 0 x 0 matrix (or an empty cell or text).
_smite_empty(v) = (v isa AbstractArray || v isa AbstractString) && isempty(v)

# A scalar field (XSize, PixelSize, ...) as the file stores it, or nothing when it is absent or
# empty; the caller checks it before converting it. Several values, or a value that is not a
# number, is an ArgumentError naming the field.
function _smite_scalar(where, s, key)
    v = get(s, key, nothing)
    v === nothing && return nothing
    k = _smite_length(v)
    k == 0 && return nothing
    k == 1 && first(v) isa Number && return first(v)
    throw(ArgumentError(k > 1 ? "$where: $key holds $k values; SMLMio loads files with one $key " *
                                "(SMITE's catSMD stores a list of a field that differs between the " *
                                "datasets it joins, such as NFrames when they differ in length)" :
                                "$where: $key is not a number"))
end

# The pixel size as MATLAB's single or double: a Float32 (a MATLAB single PixelSize) stays
# Float32 and gives a Float32 camera, as in SMLMData 0.7.1; every other real number becomes a
# Float64. It must be real, finite and positive.
function _smite_pixel_size(where, what, v)
    v isa Number || throw(ArgumentError("$where: $what is not a number"))
    iszero(imag(v)) || throw(ArgumentError("$where: $what $v is not a real number"))
    r = real(v)
    px = r isa Float32 ? r : Float64(r)
    isfinite(px) || throw(ArgumentError("$where: $what $px is not finite"))
    px > 0 || throw(ArgumentError("$where: pixel size $px is not positive"))
    return px
end

# Per-localization fields every SMITE fit writes; a file without one of them cannot be loaded.
const _SMITE_REQUIRED = ["X", "Y", "Photons", "Bg", "X_SE", "Y_SE", "FrameNum"]
const _SMITE_REQUIRED_3D = ["Z", "Z_SE"]
# Per-localization fields SMITE may leave empty (ConnectID before frame connection) or that some
# older SMITE versions wrote with one value per dataset (Photons_SE, Bg_SE). Such a field loads as
# its default here, and its raw value stays in the metadata.
const _SMITE_OPTIONAL = ("Photons_SE" => NaN, "Bg_SE" => NaN, "DatasetNum" => 1, "ConnectID" => 0)
# Metadata entries the loaders make; save_smite never writes them, and a file field of one of
# these names is not loaded.
const _SMITE_LOADER_KEYS = ("original_file", "data_size", "pixel_size", "smite_defaulted_fields",
    "smite_per_localization_fields", "smite_row_frames", "smite_row_photons", "complex_fields_removed",
    "complex_fields", "original_emitter_count", "removed_emitter_count")
# The SMD fields SMITE stacks with one value per localization (catSMD.m's StackableFields).
# Those the loaders keep in metadata, when they hold one value per localization, follow the
# emitters through save_smite.
const _SMITE_STACKABLE = ("XBoxCorner", "YBoxCorner", "X", "Y", "Z", "X_SE", "Y_SE", "Z_SE",
    "Photons", "Photons_SE", "Bg", "Bg_SE", "PSFSigma", "PSFSigmaX", "PSFSigmaY", "PSFSigma_SE",
    "PSFSigmaX_SE", "PSFSigmaY_SE", "PValue", "LogLikelihood", "ThreshFlag", "NCombined", "IndSMD",
    "DatasetNum", "FrameNum", "ConnectID", "PreGCConnectID", "RegError")

# Shared body of load_smite_2d and load_smite_3d.
function _load_smite(smd::SmiteSMD, three_d::Bool, pixel_size)
    p = joinpath(smd.filepath, smd.filename)
    # The do form closes the file on any error
    s = matopen(p) do file
        smd.varname in keys(file) ||
            throw(ArgumentError("$(smd.filename) has no variable $(smd.varname) (it has $(join(sort(collect(keys(file))), ", ")))"))
        read(file, smd.varname)
    end
    s isa AbstractDict || throw(ArgumentError("$(smd.filename): $(smd.varname) is not a structure"))
    where = "$(smd.filename) ($(smd.varname))"

    haskey(s, "X") || throw(ArgumentError("$where has no X field"))
    n = _smite_length(s["X"])

    # Per-localization columns
    cols = Dict{String,Vector}()
    for key in (three_d ? vcat(_SMITE_REQUIRED, _SMITE_REQUIRED_3D) : _SMITE_REQUIRED)
        v = get(s, key, nothing)
        v === nothing && throw(ArgumentError("$where has no $key field"))
        c = _smite_numeric(v)
        c === nothing && throw(ArgumentError("$where: $key is not numeric"))
        length(c) == n || throw(ArgumentError("$where: $key holds $(length(c)) values for $n localizations"))
        cols[key] = c
    end
    defaulted = String[]
    mismatched = String[]
    for (key, default) in _SMITE_OPTIONAL
        v = get(s, key, nothing)
        c = v === nothing ? nothing : _smite_numeric(v)
        # Present and not empty, an optional field must be numeric
        v !== nothing && c === nothing && !_smite_empty(v) &&
            throw(ArgumentError("$where: $key is not numeric"))
        if c !== nothing && length(c) == n
            cols[key] = c
        else
            cols[key] = fill(default, n)
            push!(defaulted, key)
            c !== nothing && !isempty(c) && push!(mismatched, key)
        end
    end
    !three_d && _is_ndims3(get(s, "NDims", nothing)) && n > 0 && _smite_length(get(s, "Z", nothing)) == n &&
        @warn("$where has NDims 3; load_smite_2d leaves out Z (kept in metadata), load_smite_3d loads it")
    isempty(mismatched) || @warn("$where: $(join(mismatched, ", ")) do not hold one value per " *
        "localization ($n); loaded as their defaults (NaN, or 1 for DatasetNum, 0 for ConnectID), " *
        "raw values kept in metadata")

    # Pixel size (microns) and image size
    raw_px = _smite_scalar(where, s, "PixelSize")
    file_px = raw_px === nothing ? nothing : _smite_pixel_size(where, "PixelSize", raw_px)
    arg_px = pixel_size === nothing ? nothing : _smite_pixel_size(where, "pixel_size", pixel_size)
    px = if file_px === nothing
        arg_px === nothing && throw(ArgumentError(
            "$where has no PixelSize; pass pixel_size (microns per pixel) to load it"))
        arg_px
    else
        # rtol 1e-6: a PixelSize stored as MATLAB single matches the decimal it was written from
        arg_px === nothing || isapprox(arg_px, file_px; rtol=1e-6) || throw(ArgumentError(
            "$where has PixelSize $file_px µm, but pixel_size = $pixel_size was passed"))
        file_px
    end
    xsize = _smite_scalar(where, s, "XSize")
    ysize = _smite_scalar(where, s, "YSize")
    (xsize === nothing || ysize === nothing) && throw(ArgumentError(
        "$where has no image size (XSize, YSize absent or empty; $n localizations), so its camera " *
        "cannot be built"))
    nframes = _smite_scalar(where, s, "NFrames")
    ndatasets = _smite_scalar(where, s, "NDatasets")
    (nframes === nothing || ndatasets === nothing) && throw(ArgumentError(
        "$where has no NFrames or NDatasets"))
    xsize, ysize = _smite_count(where, "XSize", xsize), _smite_count(where, "YSize", ysize)
    nframes = _smite_count(where, "NFrames", nframes, _SMITE_COUNT_MIN["NFrames"])
    ndatasets = _smite_count(where, "NDatasets", ndatasets, _SMITE_COUNT_MIN["NDatasets"])
    camera = IdealCamera(1:xsize, 1:ysize, px)
    # In Float32, a large image or pixel size can give edges that are not finite or that round
    # onto each other
    increasing(e) = all(isfinite, e) && all(k -> e[k+1] > e[k], 1:length(e)-1)
    increasing(camera.pixel_edges_x) && increasing(camera.pixel_edges_y) || throw(ArgumentError(
        "$where: $xsize × $ysize pixels of $px µm do not give finite, increasing pixel edges " *
        "in $(typeof(px))"))

    # Screen complex values in the per-localization columns read from the file
    fields_to_check = filter(k -> !(k in defaulted),
        ["X", "Y", "Photons", "Bg", "X_SE", "Y_SE", "Photons_SE", "Bg_SE"])
    three_d && append!(fields_to_check, ["Z", "Z_SE"])
    complex_indices = _complex_rows(cols, fields_to_check)
    has_complex = !isempty(complex_indices)
    if has_complex
        valid_indices = setdiff(1:n, reduce(vcat, values(complex_indices)))
        removed_count = n - length(valid_indices)
        @warn "Found $(removed_count) emitters with non-zero imaginary components. These will be excluded from the result." fields_with_complex=collect(keys(complex_indices))
    else
        valid_indices = 1:n
    end

    # Frames, datasets and track ids become Int: each value must be a whole number. FrameNum is
    # kept for every file row (the row record save_smite checks), the others for the loaded rows.
    for (key, rows) in (("FrameNum", 1:n), ("DatasetNum", valid_indices), ("ConnectID", valid_indices))
        bad = findfirst(i -> !_fits_int(cols[key][i]), rows)
        bad === nothing || throw(ArgumentError(
            "$where: $key must hold whole numbers an Int holds, but row $(rows[bad]) is $(cols[key][rows[bad]])"))
    end

    # SMITE stores X, Y and their standard errors in pixels, with (1, 1) the center of the top-left
    # pixel; SMLMCore's positions are microns from the camera's corner. Z and Z_SE are microns in both.
    f(key, i) = Float64(real(cols[key][i]))
    E = three_d ? Emitter3DFit{Float64} : Emitter2DFit{Float64}
    emitters = Vector{E}(undef, length(valid_indices))
    for (new_idx, i) in enumerate(valid_indices)
        x = (f("X", i) - 0.5) * px
        y = (f("Y", i) - 0.5) * px
        σ_x = f("X_SE", i) * px
        σ_y = f("Y_SE", i) * px
        kw = (frame=Int(cols["FrameNum"][i]), dataset=Int(cols["DatasetNum"][i]),
              track_id=Int(cols["ConnectID"][i]), id=i)
        emitters[new_idx] = if three_d
            E(x, y, f("Z", i), f("Photons", i), f("Bg", i), σ_x, σ_y, f("Z_SE", i),
              f("Photons_SE", i), f("Bg_SE", i); kw...)
        else
            E(x, y, f("Photons", i), f("Bg", i), σ_x, σ_y, f("Photons_SE", i), f("Bg_SE", i); kw...)
        end
    end

    # Metadata
    data_size = [ysize, xsize]
    # SMITE's SMD has no ZSize. load_smite_3d appends one number to data_size (a whole number of
    # at least 0); a ZSize of several values or text is only kept in the metadata, as is every
    # ZSize, so that save_smite writes it back.
    zsize = get(s, "ZSize", nothing)
    three_d && _smite_length(zsize) == 1 && first(zsize) isa Number &&
        push!(data_size, _smite_count(where, "ZSize", first(zsize), 0))
    # Every other SMITE field, and the raw value of an optional field loaded as its default
    used = Set(vcat(collect(keys(cols)), ["XSize", "YSize", "NFrames", "NDatasets", "PixelSize"]))
    setdiff!(used, defaulted)
    metadata = Dict{String,Any}(k => v for (k, v) in s if !(k in used || k in _SMITE_LOADER_KEYS))
    # SMITE's per-localization fields among them, with one value per file row: save_smite writes
    # them for each emitter's row (its id) when every emitter is still that row, which the row
    # record (each row's frame and photons) lets it check.
    perloc = sort!([k for (k, v) in metadata if k in _SMITE_STACKABLE && _smite_length(v) == n])
    # The loaders' own entries, after the file's fields so that none of these is a file's.
    metadata["original_file"] = smd.filename
    metadata["data_size"] = data_size
    metadata["pixel_size"] = px
    metadata["smite_per_localization_fields"] = perloc
    metadata["smite_row_frames"] = Int[Int(v) for v in cols["FrameNum"]]
    metadata["smite_row_photons"] = Float64[Float64(real(v)) for v in cols["Photons"]]
    isempty(defaulted) || (metadata["smite_defaulted_fields"] = defaulted)
    if has_complex
        metadata["complex_fields_removed"] = true
        metadata["complex_fields"] = collect(keys(complex_indices))
        metadata["original_emitter_count"] = n
        metadata["removed_emitter_count"] = n - length(valid_indices)
    end

    SmiteSMLD{Float64,E}(emitters, camera, nframes, ndatasets, metadata)
end

# SMITE stores counts and indices as doubles; one becomes an Int only when it is a whole number
# an Int holds, compared exactly (no rounding through Float64). A complex value counts when its
# imaginary part is zero.
_fits_int(x) = x isa Number && iszero(imag(x)) && isinteger(real(x)) &&
               typemin(Int) <= real(x) <= typemax(Int)

# Lower bounds of the counts, which save_smite applies too, so that every file it writes loads.
# An image needs a pixel; n_frames 0 means "unknown" in the lab's packages (frame connection,
# trajectories), and n_datasets 0 an SMLD without datasets.
const _SMITE_COUNT_MIN = Dict("XSize" => 1, "YSize" => 1, "NFrames" => 0, "NDatasets" => 0)

# An image size or count as an Int: a whole number of at least `lo`.
function _smite_count(where, key, v, lo=1)
    _fits_int(v) && real(v) >= lo ||
        throw(ArgumentError("$where: $key must be a whole number of at least $lo, but is $v"))
    return Int(real(v))
end

_is_ndims3(v) = _smite_length(v) == 1 && first(v) isa Number && first(v) == 3

"""
    load_smite_2d(smd::SmiteSMD; pixel_size=nothing) -> SmiteSMLD

Load a 2D SMITE SMD into a [`SmiteSMLD`](@ref) of `Emitter2DFit{Float64}`.

SMITE stores `X`, `Y`, `X_SE` and `Y_SE` in pixels, with (1, 1) the center of the top-left pixel.
They load in µm from the camera's corner: `x = (X - 0.5) * PixelSize`,
`σ_x = X_SE * PixelSize`. The camera is an `IdealCamera` of `XSize` by `YSize` pixels:
`IdealCamera{Float32}` for a MATLAB single `PixelSize`, otherwise `IdealCamera{Float64}`. How
the loaders treat empty and per-dataset fields, and what is kept in the metadata, is on the
[SMITE files](@ref) page.

# Arguments
- `smd::SmiteSMD`: the file and variable to load.

# Keyword Arguments
- `pixel_size=nothing`: µm per pixel, needed only when the file has no `PixelSize` (absent or
  empty); when the file has one, a different `pixel_size` is an error.

# Returns
- [`SmiteSMLD`](@ref): the emitters in µm, with `n_frames` frames per dataset and the SMD's
  other fields in `metadata`.

# Notes
- Throws an `ArgumentError` when the file has no pixel size and none is passed, no image size
  (`XSize`, `YSize`; SMITE leaves them empty in a result without localizations), no `NFrames` or
  `NDatasets` or several of them, or a required field (`X`, `Y`, `Photons`, `Bg`, `X_SE`,
  `Y_SE`, `FrameNum`) without one value per localization; when `XSize` or `YSize` is not one
  whole number of at least 1, `NFrames` or `NDatasets` not one of at least 0, `FrameNum`, `DatasetNum` or
  `ConnectID` holds a value that is not a whole number in `Int`'s range (`DatasetNum` and
  `ConnectID` in the loaded rows); and when the pixel size (the file's or `pixel_size`) is
  not a positive, finite real number.
- Emitters with non-zero imaginary components are excluded with a warning.

# Example
```jldoctest
julia> using SMLMCore: IdealCamera, Emitter2DFit, BasicSMLD

julia> smld = BasicSMLD([Emitter2DFit{Float64}(1.0, 2.0, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0)], IdealCamera(64, 64, 0.1), 1, 1);

julia> dir = mktempdir(); save_smite(smld, dir, "localizations.mat")

julia> loaded = load_smite_2d(SmiteSMD(dir, "localizations.mat")); loaded.metadata["pixel_size"]
0.1

julia> e = loaded.emitters[1]; (e.x, e.y)
(1.0, 2.0)
```

# See also
[`load_smite_3d`](@ref), [`save_smite`](@ref), [`SmiteSMD`](@ref)
"""
load_smite_2d(smd::SmiteSMD; pixel_size=nothing) = _load_smite(smd, false, pixel_size)

"""
    load_smite_3d(smd::SmiteSMD; pixel_size=nothing) -> SmiteSMLD

Load a 3D SMITE SMD into a [`SmiteSMLD`](@ref) of `Emitter3DFit{Float64}`.

Positions load as in [`load_smite_2d`](@ref); `Z` and `Z_SE` are µm in SMITE and load
unchanged. The file also needs `Z` and `Z_SE` with one value per localization. SMITE's SMD has
no `ZSize`; when a file has one, it is appended to `metadata["data_size"]`, and a `ZSize` that
is not a whole number is an `ArgumentError`. Everything else (optional fields, metadata,
`pixel_size`, the errors it throws) is as in [`load_smite_2d`](@ref).

# Arguments
- `smd::SmiteSMD`: the file and variable to load.

# Keyword Arguments
- `pixel_size=nothing`: µm per pixel, needed only when the file has no `PixelSize`; when the
  file has one, a different `pixel_size` is an error.

# Returns
- [`SmiteSMLD`](@ref): the emitters in µm.
"""
load_smite_3d(smd::SmiteSMD; pixel_size=nothing) = _load_smite(smd, true, pixel_size)
