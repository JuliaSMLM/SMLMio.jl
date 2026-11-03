# SMITE files

SMITE is a MATLAB-based format commonly used in single molecule localization microscopy. SMLMio provides types and functions for interfacing with this format.

## Types

```julia
# Helper type for loading SMITE .mat files
struct SmiteSMD
    filepath::String    # Path to the directory containing the .mat file
    filename::String    # Name of the .mat file
    varname::String     # Variable name in the .mat file (default: "SMD")
end

# Constructor with default variable name
SmiteSMD(filepath::String, filename::String) = SmiteSMD(filepath, filename, "SMD")

# SMLD type compatible with SMITE format
struct SmiteSMLD{T,E<:AbstractEmitter} <: AbstractSMLD
    emitters::Vector{E}
    camera::AbstractCamera
    n_frames::Int
    n_datasets::Int
    metadata::Dict{String,Any}
end
```

## Loading SMITE Data

SMLMio provides separate functions for loading 2D and 3D data from SMITE format:

```julia
using SMLMio

# Load 2D data
smd = SmiteSMD("path/to/data", "localizations.mat")
smld_2d = load_smite_2d(smd)

# Load 3D data
smld_3d = load_smite_3d(smd)
```

Both functions handle:

- Conversion of coordinates to microns
- Creation of appropriate camera model
- Preservation of metadata

Example:

```julia
using SMLMio

# Load SMITE data
smd = SmiteSMD("/data/microscopy", "cell1_localizations.mat")
smld = load_smite_2d(smd)

# Inspect metadata
println("Loaded $(length(smld.emitters)) emitters")
println("Image size: $(smld.metadata["data_size"])")
println("Pixel size: $(smld.metadata["pixel_size"]) μm")
```

## Saving SMITE Data

Save a SmiteSMLD object back to SMITE format:

```julia
# Save to SMITE format
save_smite(smld, "output/directory", "processed_results.mat")
```

The `save_smite` function:

- Converts SmiteSMLD data to SMITE's structure
- Preserves metadata fields whose names are not emitter columns
- Saves in MATLAB v7.3 format

Example workflow:

```julia
using SMLMio
using SMLMCore: @filter, filter_roi
using Dates

# Load data
smd = SmiteSMD("raw_data", "experiment1.mat")
smld = load_smite_2d(smd)

# Process data
bright = @filter(smld, photons > 1000)
roi = filter_roi(bright, 10.0:20.0, 10.0:20.0)

# Add analysis info to metadata
roi.metadata["analysis_date"] = Dates.now()
roi.metadata["selection_criteria"] = "photons > 1000, ROI (10-20,10-20)"

# Save processed data
save_smite(roi, "processed_data", "experiment1_processed.mat")
```
