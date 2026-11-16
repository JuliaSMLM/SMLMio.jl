# Tutorial

This tutorial builds a small dataset, saves it as a SMITE file, loads it back, filters it and
saves the result. It uses constructed emitters and a temporary folder, so it needs no data.
Positions are in µm; the file holds them in SMITE's pixels, and the loader and saver convert
([SMITE files](@ref)).

## Build a dataset

Three fitted emitters on a 64 × 64 pixel camera with 0.1 µm pixels, one per frame:

```@example tutorial
using SMLMio
using SMLMCore: IdealCamera, Emitter2DFit, BasicSMLD, @filter

camera = IdealCamera(1:64, 1:64, 0.1)   # pixel size in µm
emitters = [
    Emitter2DFit{Float64}(1.0, 2.0, 900.0, 10.0, 0.01, 0.01, 50.0, 2.0; frame=1),
    Emitter2DFit{Float64}(3.0, 4.0, 1200.0, 12.0, 0.02, 0.02, 60.0, 2.0; frame=2),
    Emitter2DFit{Float64}(5.0, 6.0, 1500.0, 12.0, 0.02, 0.02, 60.0, 2.0; frame=3),
]
smld = BasicSMLD(emitters, camera, 3, 1, Dict{String,Any}())   # 3 frames, 1 dataset
```

## Save and load

[`save_smite`](@ref) writes the SMLD into a folder, and [`SmiteSMD`](@ref) names the file for
[`load_smite_2d`](@ref):

```@example tutorial
dir = mktempdir()
save_smite(smld, dir, "localizations.mat")

loaded = load_smite_2d(SmiteSMD(dir, "localizations.mat"))
```

The loaders return a [`SmiteSMLD`](@ref). Its emitters carry the saved values, in µm:

```@example tutorial
[(e.x, e.y, e.photons) for e in loaded.emitters]
```

## Filter and save again

SMLMCore's `@filter` keeps the emitters that meet a condition, here more than 1000 photons.
The result is still a `SmiteSMLD`, and saving it writes only the kept localizations:

```@example tutorial
bright = @filter(loaded, photons > 1000)
save_smite(bright, dir, "bright.mat")

again = load_smite_2d(SmiteSMD(dir, "bright.mat"))
[e.x for e in again.emitters]
```

## Next

- [SMITE files](@ref): what the loaders read and what `save_smite` writes.
- [API](api.md): every exported name.
