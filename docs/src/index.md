# SMLMio.jl

SMLMio reads and writes the lab's file formats for single-molecule localization microscopy (SMLM)
data. It provides SMITE `.mat` file I/O. It was split from SMLMData 0.7.1.

!!! info "Units"
    Positions and their uncertainties are in microns, as in SMLMCore.

## Where it fits

SMLMio takes and returns SMLMCore's types. SMLMData loads and re-exports SMLMCore and SMLMio,
so existing `using SMLMData` code is unchanged by the split (SMLMData 0.7.2 and later).

## Installation

SMLMio and SMLMCore are not registered yet. Install both from GitHub, SMLMCore first:

```julia
using Pkg
Pkg.add(url="https://github.com/JuliaSMLM/SMLMCore.jl")
Pkg.add(url="https://github.com/JuliaSMLM/SMLMio.jl")
```

## Quick start

Save a small dataset to a temporary folder. Positions and their uncertainties are in microns.

```@example quickstart
using SMLMio
using SMLMCore: IdealCamera, Emitter2DFit

camera = IdealCamera(1:64, 1:64, 0.1) # pixel size in microns
emitters = [Emitter2DFit{Float64}(1.0, 2.0, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0)]
smld = SmiteSMLD{Float64,Emitter2DFit{Float64}}(
    emitters, camera, 1, 1, Dict{String,Any}(),
)
mktempdir() do dir
    save_smite(smld, dir, "localizations.mat")
    println(isfile(joinpath(dir, "localizations.mat")))
end
```

See [SMITE files](smite.md) for file-based workflows and the [API](api.md) for all exported names.
