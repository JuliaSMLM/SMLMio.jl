# SMLMio.jl

SMLMio reads and writes the lab's file formats for single-molecule localization microscopy (SMLM)
data. It provides SMITE `.mat` file I/O.

!!! info "Units"
    Positions and their uncertainties are in microns, as in SMLMCore.

## Where it fits

SMLMio takes and returns [SMLMCore](https://JuliaSMLM.github.io/SMLMCore.jl/stable)'s types. From
SMLMData 0.7.2 (not yet released), SMLMData re-exports both, so `using SMLMData` code keeps working,
with SMITE results changed as the [CHANGELOG](https://github.com/JuliaSMLM/SMLMio.jl/blob/main/CHANGELOG.md) lists.

## Installation

Install SMLMio from Julia's General registry:

```julia
using Pkg
Pkg.add("SMLMio")
Pkg.add("SMLMCore")   # the examples build data with SMLMCore's types
```

An environment that holds SMLMData 0.7.1 or earlier, directly or through a package that depends
on it (GaussMLE, SMLMSim, SMLMBoxer, SMLMAnalysis and others), cannot add SMLMio. Those versions
define the same types and their own SMITE loaders, and Pkg reports unsatisfiable requirements for
SMLMData. Until SMLMData 0.7.2 is released, use SMLMio in an environment without them.

## Quick start

Save a small dataset as a SMITE file and load it back. Positions and their uncertainties are in
microns; the file holds them in SMITE's pixels.

```@example quickstart
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

See the [Tutorial](tutorial.md) for a worked example, [SMITE files](guide/smite.md) for file-based
workflows and the [API](api.md) for all exported names.
