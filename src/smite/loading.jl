# For each of `fields` in `cols` (vectors with one value per localization) that holds a non-zero
# imaginary part, the rows where it does.
function _complex_rows(cols, fields)
    rows = Dict{String,Vector{Int}}()
    for f in fields
        v = cols[f]
        eltype(v) <: Complex || continue
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

# A scalar field (XSize, PixelSize, ...) as a Float64, or nothing when it is absent or empty.
function _smite_scalar(s, key)
    v = get(s, key, nothing)
    _smite_length(v) == 1 || return nothing
    return Float64(real(first(v)))
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
    "smite_per_localization_fields", "complex_fields_removed", "complex_fields",
    "original_emitter_count", "removed_emitter_count")
# SMITE's SMD fields that hold one value for the whole file.
const _SMITE_SCALARS = ("NDims", "NFrames", "NDatasets", "FrameRate", "PixelSize", "XSize", "YSize",
    "ZSize", "ZOffset", "IsTransformed", "RegError")

# A field with one value per localization: a vector (or n x 1 or 1 x n) of n values for n >= 2;
# for a single localization, which MAT reads as scalars, a number that is not one of SMITE's
# whole-file fields. With n equal to NDatasets, a per-dataset field cannot be told apart.
function _smite_per_localization(key, v, n)
    n >= 2 && return v isa AbstractArray && ndims(v) <= 2 && length(v) == n && (v isa AbstractVector || 1 in size(v))
    return n == 1 && v isa Number && !(key in _SMITE_SCALARS)
end

# Shared body of load_smite_2d and load_smite_3d.
function _load_smite(smd::SmiteSMD, three_d::Bool, pixel_size)
    p = joinpath(smd.filepath, smd.filename)
    file = matopen(p)
    s = read(file, smd.varname)
    close(file)
    where = "$(smd.filename) ($(smd.varname))"

    haskey(s, "X") || throw(ArgumentError("$where has no X field"))
    n = _smite_length(s["X"])

    # Per-localization columns
    cols = Dict{String,Vector}()
    for key in (three_d ? vcat(_SMITE_REQUIRED, _SMITE_REQUIRED_3D) : _SMITE_REQUIRED)
        v = get(s, key, nothing)
        v === nothing && throw(ArgumentError("$where has no $key field"))
        _smite_length(v) == n || throw(ArgumentError(
            "$where: $key holds $(_smite_length(v)) values for $n localizations"))
        cols[key] = _smite_vec(v)
    end
    defaulted = String[]
    mismatched = String[]
    for (key, default) in _SMITE_OPTIONAL
        v = get(s, key, nothing)
        if v !== nothing && _smite_length(v) == n
            cols[key] = _smite_vec(v)
        else
            cols[key] = fill(default, n)
            push!(defaulted, key)
            v !== nothing && _smite_length(v) > 0 && push!(mismatched, key)
        end
    end
    !three_d && _smite_scalar(s, "NDims") == 3 && n > 0 && _smite_length(get(s, "Z", nothing)) == n &&
        @warn("$where has NDims 3; load_smite_2d leaves out Z (kept in metadata), load_smite_3d loads it")
    isempty(mismatched) || @warn("$where: $(join(mismatched, ", ")) do not hold one value per " *
        "localization ($n); loaded as their defaults (NaN, or 1 for DatasetNum, 0 for ConnectID), " *
        "raw values kept in metadata")

    # Pixel size (microns) and image size
    file_px = _smite_scalar(s, "PixelSize")
    px = if file_px === nothing
        pixel_size === nothing && throw(ArgumentError(
            "$where has no PixelSize; pass pixel_size (microns per pixel) to load it"))
        Float64(pixel_size)
    else
        # rtol 1e-6: a PixelSize stored as MATLAB single matches the decimal it was written from
        pixel_size === nothing || isapprox(pixel_size, file_px; rtol=1e-6) || throw(ArgumentError(
            "$where has PixelSize $file_px µm, but pixel_size = $pixel_size was passed"))
        file_px
    end
    px > 0 || throw(ArgumentError("$where: pixel size $px is not positive"))
    xsize = _smite_scalar(s, "XSize")
    ysize = _smite_scalar(s, "YSize")
    (xsize === nothing || ysize === nothing) && throw(ArgumentError(
        "$where has no image size (XSize, YSize absent or empty; $n localizations), so its camera " *
        "cannot be built"))
    nframes = _smite_scalar(s, "NFrames")
    ndatasets = _smite_scalar(s, "NDatasets")
    (nframes === nothing || ndatasets === nothing) && throw(ArgumentError(
        "$where has no NFrames or NDatasets"))
    camera = IdealCamera(1:Int(xsize), 1:Int(ysize), px)

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
    data_size = [Int(ysize), Int(xsize)]
    zsize = _smite_scalar(s, "ZSize")   # SMITE's SMD has no ZSize; kept when a file has one
    three_d && zsize !== nothing && push!(data_size, Int(zsize))
    # Every other SMITE field, and the raw value of an optional field loaded as its default
    used = Set(vcat(collect(keys(cols)), ["XSize", "YSize", "NFrames", "NDatasets", "PixelSize"],
                    three_d ? ["ZSize"] : String[]))
    setdiff!(used, defaulted)
    metadata = Dict{String,Any}(k => v for (k, v) in s if !(k in used || k in _SMITE_LOADER_KEYS))
    # The loaders' own entries, after the file's fields so that none of these is a file's.
    # smite_per_localization_fields lists the fields kept above that hold one value per file row,
    # which save_smite writes back by each emitter's row (its id).
    metadata["original_file"] = smd.filename
    metadata["data_size"] = data_size
    metadata["pixel_size"] = px
    metadata["smite_per_localization_fields"] =
        sort!([k for (k, v) in metadata if _smite_per_localization(k, v, n)])
    isempty(defaulted) || (metadata["smite_defaulted_fields"] = defaulted)
    if has_complex
        metadata["complex_fields_removed"] = true
        metadata["complex_fields"] = collect(keys(complex_indices))
        metadata["original_emitter_count"] = n
        metadata["removed_emitter_count"] = n - length(valid_indices)
    end

    SmiteSMLD{Float64,E}(emitters, camera, Int(nframes), Int(ndatasets), metadata)
end

"""
    load_smite_2d(smd::SmiteSMD; pixel_size=nothing)

Load a 2D SMITE SMD .mat file into a `SmiteSMLD` of `Emitter2DFit{Float64}`.

SMITE stores `X`, `Y`, `X_SE` and `Y_SE` in pixels, with (1, 1) the center of the top-left pixel.
They load in microns from the camera's corner: `x = (X - 0.5) * PixelSize`,
`σ_x = X_SE * PixelSize`. The camera is an `IdealCamera` of `XSize` by `YSize` pixels.

# Arguments
- `smd::SmiteSMD`: SmiteSMD object specifying the file to load
- `pixel_size`: microns per pixel, needed only when the file has no `PixelSize` (absent or
  empty); when the file has one, a different `pixel_size` is an error

# Returns
SmiteSMLD containing 2D localizations

# Notes
- `Photons_SE` and `Bg_SE` load as `NaN`, `DatasetNum` as 1 and `ConnectID` as 0 when the file
  leaves them empty (SMITE leaves `ConnectID` empty before frame connection) or holds them with
  one value per dataset (some older SMITE versions); their raw values stay in the metadata, and
  `metadata["smite_defaulted_fields"]` lists them. A non-empty field of the wrong length warns.
- Every other SMD field is kept in the metadata as SMITE wrote it (for example `DriftX` in pixels).
- Emitters with non-zero imaginary components are excluded with a warning.
- Throws an `ArgumentError` when the file has no pixel size and none is passed, no image size
  (`XSize`, `YSize`; SMITE leaves them empty in a result without localizations), no `NFrames`
  or `NDatasets`, or a required field (`X`, `Y`, `Photons`, `Bg`, `X_SE`, `Y_SE`, `FrameNum`)
  without one value per localization.
"""
load_smite_2d(smd::SmiteSMD; pixel_size=nothing) = _load_smite(smd, false, pixel_size)

"""
    load_smite_3d(smd::SmiteSMD; pixel_size=nothing)

Load a 3D SMITE SMD .mat file into a `SmiteSMLD` of `Emitter3DFit{Float64}`.

Positions load as in [`load_smite_2d`](@ref); `Z` and `Z_SE` are microns in SMITE and load
unchanged. The file also needs `Z` and `Z_SE` with one value per localization. SMITE's SMD has
no `ZSize`; when a file has one, it is appended to `metadata["data_size"]`.

# Arguments
- `smd::SmiteSMD`: SmiteSMD object specifying the file to load
- `pixel_size`: microns per pixel, needed only when the file has no `PixelSize`

# Returns
SmiteSMLD containing 3D localizations
"""
load_smite_3d(smd::SmiteSMD; pixel_size=nothing) = _load_smite(smd, true, pixel_size)
