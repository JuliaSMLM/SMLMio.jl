"""
    SMLMio

File reading and writing for single-molecule localization microscopy (SMLM) data.

SMLMio holds the file formats; the data types it reads into and writes from (emitters, cameras,
`BasicSMLD`) live in SMLMCore. Positions are in microns, as in SMLMCore.

# Formats
- SMITE, the lab's MATLAB package: `SmiteSMD` names a `.mat` file, `load_smite_2d` and
  `load_smite_3d` read it into a `SmiteSMLD`, and `save_smite` writes one back.

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
# Internal helpers SMLMCore uses for its own show methods; SmiteSMLD's show reuses them.
using SMLMCore: _dim_label, format_with_commas

export SmiteSMD, SmiteSMLD, load_smite_2d, load_smite_3d, save_smite

include("smite/types.jl")
include("smite/loading.jl")
include("smite/saving.jl")

end
