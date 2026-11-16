"""
    SMLMio

File reading and writing for single-molecule localization microscopy (SMLM) data.

The emitters and cameras it reads and writes are SMLMCore's. The loaders return SMLMio's
[`SmiteSMLD`](@ref), an SMLMCore `AbstractSMLD`, and [`save_smite`](@ref) writes any SMLD of fit
emitters. Positions are in µm, as in SMLMCore.

# Formats
- SMITE, the lab's MATLAB package: [`SmiteSMD`](@ref) names a `.mat` file, [`load_smite_2d`](@ref)
  and [`load_smite_3d`](@ref) read it into a [`SmiteSMLD`](@ref), and [`save_smite`](@ref) writes
  one back. See [SMITE files](@ref).

# Example
```julia
using SMLMio

smd = SmiteSMD("path/to/data", "localizations.mat")
smld = load_smite_2d(smd)
save_smite(smld, "path/to/output", "copy.mat")
```
"""
module SMLMio

using MAT: matopen, matwrite
using SMLMCore: AbstractCamera, AbstractEmitter, AbstractSMLD, Emitter2DFit, Emitter3DFit,
    IdealCamera, emitter_ndims

export SmiteSMD, SmiteSMLD, load_smite_2d, load_smite_3d, save_smite

include("smite/types.jl")
include("smite/loading.jl")
include("smite/saving.jl")

end
