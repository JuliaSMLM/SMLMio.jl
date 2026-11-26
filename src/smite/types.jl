"""
    SmiteSMD(filepath, filename, varname="SMD")

Name of an SMD variable in a SMITE `.mat` file: folder, file name and variable name.

# Fields
- `filepath::String`: folder holding the file.
- `filename::String`: name of the `.mat` file.
- `varname::String`: name of the SMD variable in the file.

# Constructors
- `SmiteSMD(filepath, filename)`: the variable is `"SMD"`, as SMITE names it.
- `SmiteSMD(filepath, filename, varname)`: a variable of another name.

Each argument can be any `AbstractString`; the fields are `String`s.

# Example
```jldoctest
julia> SmiteSMD("path/to/data", "localizations.mat")
SmiteSMD:
  Filepath: "path/to/data"
  Filename: "localizations.mat"
  Variable name: "SMD"
```

# See also
[`load_smite_2d`](@ref), [`load_smite_3d`](@ref)
"""
mutable struct SmiteSMD 
    filepath::String
    filename::String
    varname::String
end

# Any string type for the paths, and the default variable name
SmiteSMD(filepath::AbstractString, filename::AbstractString, varname::AbstractString) =
    SmiteSMD(String(filepath), String(filename), String(varname))
SmiteSMD(filepath::AbstractString, filename::AbstractString) = SmiteSMD(filepath, filename, "SMD")

# Show methods for SmiteSMD: the constructor call that rebuilds it
function Base.show(io::IO, smd::SmiteSMD)
    print(io, "SmiteSMD(", repr(smd.filepath), ", ", repr(smd.filename), ", ", repr(smd.varname), ")")
end

function Base.show(io::IO, ::MIME"text/plain", smd::SmiteSMD)
    println(io, "SmiteSMD:")
    println(io, "  Filepath: \"$(smd.filepath)\"")
    println(io, "  Filename: \"$(smd.filename)\"")
    print(io, "  Variable name: \"$(smd.varname)\"")
end

"""
    SmiteSMLD{T,E<:AbstractEmitter} <: AbstractSMLD

SMLD of the fit emitters of a SMITE SMD; the SMD's other fields are kept in `metadata`.

# Fields
- `emitters::Vector{E}`: the localizations; positions and uncertainties in µm.
- `camera::AbstractCamera`: an `IdealCamera` of `XSize` × `YSize` pixels of `PixelSize` µm, built
  from the file (an sCMOS calibration is not in the file); `Float32` for a MATLAB single
  `PixelSize` (or a `Float32` `pixel_size` when the file has none), otherwise `Float64`.
- `n_frames::Int`: frames per dataset (SMITE's `NFrames`); in SMITE's convention `frame` runs
  from 1 to `n_frames` in each dataset, which the loaders do not check. 0 means "unknown".
- `n_datasets::Int`: number of datasets (SMITE's `NDatasets`).
- `metadata::Dict{String,Any}`: the SMD's other fields and the loaders' entries; see
  [SMITE files](@ref).

# Constructors
Returned by [`load_smite_2d`](@ref) and [`load_smite_3d`](@ref), with `T` = `Float64`, the
emitters' number type (the camera may be `Float32`, as above). To build
one, give both type parameters: `SmiteSMLD{Float64,E}(emitters, camera, n_frames, n_datasets,
metadata)`. [`save_smite`](@ref) takes any SMLD of fit emitters, so building one is rarely needed.
"""
struct SmiteSMLD{T,E<:AbstractEmitter} <: AbstractSMLD
    emitters::Vector{E}
    camera::AbstractCamera
    n_frames::Int
    n_datasets::Int
    metadata::Dict{String,Any}
end

# Display helpers, copies of SMLMCore's private ones so SMLMio depends only on SMLMCore's public names.
# Dimension word for display; mixed or empty non-concrete data are labelled as such rather than as 3D.
function _dim_label(emitters)
    d = emitter_ndims(emitters)
    d === nothing && return isempty(emitters) ? "unknown-dimension" : "mixed 2D/3D"
    return string(d, "D")
end

# An integer with thousands separators, for display.
function format_with_commas(n::Integer)
    return replace(string(n), r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
end

# Show method for SmiteSMLD
function Base.show(io::IO, ::MIME"text/plain", smld::SmiteSMLD{T,E}) where {T,E}
    n_emitters = length(smld.emitters)
    
    # Camera info
    cam = smld.camera
    n_pixels_x = length(cam.pixel_edges_x) - 1
    n_pixels_y = length(cam.pixel_edges_y) - 1
    
    # Find mean photons per emitter
    if n_emitters > 0
        mean_photons = sum(e.photons for e in smld.emitters) / n_emitters
        mean_photons = round(mean_photons, digits=1)
    else
        mean_photons = "N/A"
    end
    
    # Get localization dimensions (2D or 3D)
    dim_str = _dim_label(smld.emitters)
    
    # Extract original file if available
    orig_file = get(smld.metadata, "original_file", "unknown")
    
    println(io, "SmiteSMLD{$T,$E}:")
    println(io, "  Emitters: $(format_with_commas(n_emitters)) $dim_str localizations")
    println(io, "  Camera: $(n_pixels_x)×$(n_pixels_y) pixels")
    println(io, "  Frames: $(format_with_commas(smld.n_frames))")
    println(io, "  Datasets: $(smld.n_datasets)")
    println(io, "  Mean photons: $mean_photons")
    println(io, "  Original file: $orig_file")
    
    # Add complex field info if present
    if get(smld.metadata, "complex_fields_removed", false)
        removed = get(smld.metadata, "removed_emitter_count", 0)
        fields = get(smld.metadata, "complex_fields", [])
        println(io, "  Complex fields removed: $removed emitters")
        print(io, "  Fields with complex values: $(join(fields, ", "))")
    end
end
