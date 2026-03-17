# YAML config examples

Examples live in this folder:

- `box_wop.yaml`
- `box_wop_legacy_equivalent.yaml`
- `box_wos.yaml`

Each file contains a top-level `command` field with the intended CLI invocation.
The current parser treats `command` as metadata, so it is safe to keep it in the config.
Boundary data for config-mode runs is defined in `app/problem_functions.cpp`.
Edit that file, rebuild `wop_cli`, then run one of the YAML configs below.

Run from `external_wop_cpp/` on Windows with the Visual Studio generator:

```powershell
.\build\Release\wop_cli.exe --config examples/box_wop.yaml --json
.\build\Release\wop_cli.exe --config examples/box_wop_legacy_equivalent.yaml --json
.\build\Release\wop_cli.exe --config examples/box_wos.yaml --json
```

If you build another configuration, replace `Release` with `Debug` or the configuration you actually built.

`box_wop_legacy_equivalent.yaml` matches the legacy `--example box` WOP defaults for
`x0`, `n`, `seed`, `max_steps`, `r_max`, `r_max_mode`, and `r_max_factor`.
For a full 1-to-1 comparison, keep `app/problem_functions.cpp` on the same boundary
function as the legacy path uses.

The active boundary and reference functions come from `app/problem_functions.cpp`.
Update that file to the function you want to benchmark before rebuilding `wop_cli`.
