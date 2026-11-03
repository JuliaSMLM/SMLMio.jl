# SMLMio

SMLMio provides the lab's SMLM file formats, currently SMITE `.mat` I/O.
It holds SMLMData 0.7.1's SMITE code unchanged; shared data types live in SMLMCore.
SMLMData loads and re-exports both packages.

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
