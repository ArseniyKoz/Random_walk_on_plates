# YAML config examples

Examples live in this folder:

- `box_wop.yaml`
- `box_wos.yaml`

Each file contains a top-level `command` field with the intended CLI invocation.
The current parser treats `command` as metadata, so it is safe to keep it in the config.
Boundary data is now selected through built-in C++ functions such as
`constant`, `x`, `y`, `z`, and `coulomb`.
The legacy `--example box` path is still the fastest specialized path.

Run from `external_wop_cpp/` on Windows with the Visual Studio generator:

```powershell
.\build\Release\wop_cli.exe --config examples/box_wop.yaml --json
.\build\Release\wop_cli.exe --config examples/box_wos.yaml --json
```

If you build another configuration, replace `Release` with `Debug` or the configuration you actually built.
