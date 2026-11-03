# SMLMio

SMLMio reads and writes the lab's file formats for single-molecule localization microscopy (SMLM)
data. It provides SMITE `.mat` file I/O using SMLMCore's emitters and cameras.

[![CI](https://github.com/JuliaSMLM/SMLMio.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaSMLM/SMLMio.jl/actions/workflows/CI.yml)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSMLM.github.io/SMLMio.jl/dev/)

Status: SMLMio 0.1 holds SMLMData 0.7.1's SMITE code unchanged. A November 2026 review found that the SMITE loaders fail on files written by current SMITE; those fixes are the next release's work.

## Installation

SMLMio and SMLMCore are not registered yet. Install both from GitHub, SMLMCore first:

```julia
using Pkg
Pkg.add(url="https://github.com/JuliaSMLM/SMLMCore.jl")
Pkg.add(url="https://github.com/JuliaSMLM/SMLMio.jl")
```

## Quick start

Save a small dataset to a temporary folder. Positions and their uncertainties are in microns.

```julia
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

## Where it fits

SMLMio takes and returns SMLMCore's types. SMLMData loads and re-exports SMLMCore and SMLMio,
so existing `using SMLMData` code is unchanged by the split (SMLMData 0.7.2 and later).

## Documentation

See the [documentation](https://JuliaSMLM.github.io/SMLMio.jl/dev/) for
[SMITE files](https://JuliaSMLM.github.io/SMLMio.jl/dev/smite/) and the
[API](https://JuliaSMLM.github.io/SMLMio.jl/dev/api/).

## License

[MIT](LICENSE).
