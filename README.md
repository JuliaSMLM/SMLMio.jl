# SMLMio

SMLMio reads and writes the lab's file formats for single-molecule localization microscopy (SMLM)
data. It provides SMITE `.mat` file I/O using SMLMCore's emitters and cameras.

[![CI](https://github.com/JuliaSMLM/SMLMio.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaSMLM/SMLMio.jl/actions/workflows/CI.yml)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSMLM.github.io/SMLMio.jl/dev/)

SMLMio 0.1 holds SMLMData 0.7.1's SMITE I/O with its bugs fixed: positions convert between
SMITE's pixels and microns, the loaders read files whose optional fields are empty or one per
dataset, and files written by `save_smite` load back. See the [CHANGELOG](CHANGELOG.md).

## Installation

SMLMio and SMLMCore are not registered yet. Install both from GitHub, SMLMCore first:

```julia
using Pkg
Pkg.add(url="https://github.com/JuliaSMLM/SMLMCore.jl")
Pkg.add(url="https://github.com/JuliaSMLM/SMLMio.jl")
```

## Quick start

Load a SMITE results file. Positions and their uncertainties come out in microns; SMITE stores
them in pixels.

```julia
using SMLMio

smld = load_smite_2d(SmiteSMD("path/to/Results/ResultsStructs", "Cell_01_Label_01_Results.mat"))
smld.emitters[1].x            # microns from the camera's corner
```

Save a small dataset as a SMITE file and load it back:

```julia
using SMLMio
using SMLMCore: IdealCamera, Emitter2DFit, BasicSMLD

camera = IdealCamera(1:64, 1:64, 0.1) # pixel size in microns
emitters = [Emitter2DFit{Float64}(1.0, 2.0, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0)]
smld = BasicSMLD(emitters, camera, 1, 1, Dict{String,Any}())
mktempdir() do dir
    save_smite(smld, dir, "localizations.mat")
    loaded = load_smite_2d(SmiteSMD(dir, "localizations.mat"))
    e = loaded.emitters[1]
    println((e.x, e.y))   # (1.0, 2.0): microns, as saved
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
