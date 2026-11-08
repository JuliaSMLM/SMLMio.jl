"""
    has_nonzero_imag(value)

Check if a value has a non-zero imaginary component.
Works for both scalar values and arrays.
"""
function has_nonzero_imag(value)
    if isa(value, Complex) || eltype(value) <: Complex
        # For arrays, check if any element has non-zero imaginary part
        if isa(value, AbstractArray)
            return any(imag.(value) .!= 0)
        else
            # For scalar complex values
            return imag(value) != 0
        end
    end
    return false  # Not complex
end

"""
    check_complex_fields(s, fields)

Check if any of the given fields in s are complex and have non-zero imaginary components.
Returns a tuple with:
1. Boolean indicating if any fields are complex with non-zero imaginary parts
2. Dict mapping field names to the linear indices with non-zero imaginary parts
"""
function check_complex_fields(s, fields)
    has_complex = false
    complex_indices = Dict{String, Vector{Int}}()

    for field in fields
        if haskey(s, field) && has_nonzero_imag(s[field])
            has_complex = true

            # Find indices with non-zero imaginary parts
            if isa(s[field], AbstractArray)
                complex_indices[field] = findall(!iszero, imag.(vec(s[field])))
            end
        end
    end

    return has_complex, complex_indices
end

"""
    get_valid_indices(s, complex_indices)

Get indices of elements that don't have complex values with non-zero imaginary parts
in any field.
"""
function get_valid_indices(s, complex_indices)
    n = _smite_length(s["X"])

    # Start with all indices
    valid_indices = collect(1:n)

    # Remove indices that have complex values in any field
    all_complex_indices = Set{Int}()
    for indices in values(complex_indices)
        union!(all_complex_indices, indices)
    end

    setdiff!(valid_indices, all_complex_indices)
    return valid_indices
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
        pixel_size === nothing || isapprox(pixel_size, file_px) || throw(ArgumentError(
            "$where has PixelSize $file_px µm, but pixel_size = $pixel_size was passed"))
        file_px
    end
    px > 0 || throw(ArgumentError("$where: pixel size $px is not positive"))
    xsize, ysize = s["XSize"], s["YSize"]
    nframes, ndatasets = s["NFrames"], s["NDatasets"]
    camera = IdealCamera(1:Int(xsize), 1:Int(ysize), px)

    # Screen complex values in the per-localization columns read from the file
    fields_to_check = filter(k -> !(k in defaulted),
        ["X", "Y", "Photons", "Bg", "X_SE", "Y_SE", "Photons_SE", "Bg_SE"])
    three_d && append!(fields_to_check, ["Z", "Z_SE"])
    has_complex, complex_indices = check_complex_fields(cols, fields_to_check)
    if has_complex
        valid_indices = get_valid_indices(cols, complex_indices)
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
    three_d && push!(data_size, Int(s["ZSize"]))
    metadata = Dict{String,Any}(
        "original_file" => smd.filename,
        "data_size" => data_size,
        "pixel_size" => px
    )
    isempty(defaulted) || (metadata["smite_defaulted_fields"] = defaulted)
    if has_complex
        metadata["complex_fields_removed"] = true
        metadata["complex_fields"] = collect(keys(complex_indices))
        metadata["original_emitter_count"] = n
        metadata["removed_emitter_count"] = n - length(valid_indices)
    end
    # Every other SMITE field, and the raw value of an optional field loaded as its default
    used = Set(vcat(collect(keys(cols)), ["XSize", "YSize", "NFrames", "NDatasets", "PixelSize"],
                    three_d ? ["ZSize"] : String[]))
    setdiff!(used, defaulted)
    for key in keys(s)
        key in used || (metadata[key] = s[key])
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
- Throws an `ArgumentError` when the file has no pixel size and none is passed, or when a
  required field (`X`, `Y`, `Photons`, `Bg`, `X_SE`, `Y_SE`, `FrameNum`) does not hold one value
  per localization.
"""
load_smite_2d(smd::SmiteSMD; pixel_size=nothing) = _load_smite(smd, false, pixel_size)

"""
    load_smite_3d(smd::SmiteSMD; pixel_size=nothing)

Load a 3D SMITE SMD .mat file into a `SmiteSMLD` of `Emitter3DFit{Float64}`.

Positions load as in [`load_smite_2d`](@ref); `Z` and `Z_SE` are microns in SMITE and load
unchanged. The file also needs `Z` and `Z_SE` with one value per localization.

# Arguments
- `smd::SmiteSMD`: SmiteSMD object specifying the file to load
- `pixel_size`: microns per pixel, needed only when the file has no `PixelSize`

# Returns
SmiteSMLD containing 3D localizations
"""
load_smite_3d(smd::SmiteSMD; pixel_size=nothing) = _load_smite(smd, true, pixel_size)
