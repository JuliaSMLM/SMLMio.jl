# SMLMio

SMLMio provides the lab's SMLM file formats, currently SMITE `.mat` I/O.
It holds SMLMData 0.7.1's SMITE code with its bugs fixed (see CHANGELOG.md); shared data types
live in SMLMCore. SMLMData loads and re-exports both packages.

SMITE stores X, Y, X_SE, Y_SE in pixels, (1, 1) the center of the top-left pixel; SMLMio's
positions are microns from the camera's corner.

## Tests

Run from the package root:

```bash
GROUP=Core julia --project -e 'using Pkg; Pkg.test()'
GROUP=QA julia --project -e 'using Pkg; Pkg.test()'
```

Record all declared test groups and the strict docs build from a clean, pushed checkout:

```bash
julia /home/kalidke/julia_shared_dev/kalidke/admiral/conventions/templates/package/record_tests.jl
```

julia = 1.11 is temporary: [sources] (Julia 1.11+) resolves the unregistered SMLMCore. When
SMLMCore is registered, drop [sources] and lower the floor to 1.9 (SMLMData's floor) with a floor record, before
SMLMio is registered.
