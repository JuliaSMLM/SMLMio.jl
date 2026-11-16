# SMLMio

SMLMio reads and writes the lab's file formats for single-molecule localization microscopy (SMLM)
data. It provides SMITE `.mat` file I/O using SMLMCore's emitters and cameras.

[![CI](https://github.com/JuliaSMLM/SMLMio.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaSMLM/SMLMio.jl/actions/workflows/CI.yml)
[![Docs: stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSMLM.github.io/SMLMio.jl/stable)

## Installation

SMLMio is not registered yet. Install it from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/JuliaSMLM/SMLMio.jl")
Pkg.add("SMLMCore")   # the examples build data with SMLMCore's types
```

An environment that holds SMLMData 0.7.1 or earlier, directly or through a package that depends
on it (GaussMLE, SMLMSim, SMLMBoxer, SMLMAnalysis and others), cannot add SMLMio. Those versions
define the same types and their own SMITE loaders, and Pkg reports unsatisfiable requirements for
SMLMData. Until SMLMData 0.7.2 is released, use SMLMio in an environment without them.

## Quick start

Save a small dataset as a SMITE file and load it back. Positions and their uncertainties are in
microns; the file holds them in SMITE's pixels.

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
    println((e.x, e.y))   # microns, as saved
end
```

## Where it fits

SMLMio takes and returns SMLMCore's types ([SMLMCore](https://JuliaSMLM.github.io/SMLMCore.jl/stable)).
From SMLMData 0.7.2 (not yet released), SMLMData loads and re-exports SMLMCore and SMLMio, so
`using SMLMData` code keeps working. Its SMITE results change as the
[CHANGELOG](https://github.com/JuliaSMLM/SMLMio.jl/blob/main/CHANGELOG.md) lists: positions are
in microns, and a file without `PixelSize` needs `pixel_size`.

## Documentation

See the [documentation](https://JuliaSMLM.github.io/SMLMio.jl/stable/) for the
[tutorial](https://JuliaSMLM.github.io/SMLMio.jl/stable/tutorial/),
[SMITE files](https://JuliaSMLM.github.io/SMLMio.jl/stable/guide/smite/) and the
[API](https://JuliaSMLM.github.io/SMLMio.jl/stable/api/).

## License

[MIT](LICENSE).
