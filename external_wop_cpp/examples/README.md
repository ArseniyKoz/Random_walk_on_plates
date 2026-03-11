# YAML config examples

Examples live in this folder:

- `box_wop.yaml`
- `box_wos.yaml`

Each file contains a top-level `command` field with the intended CLI invocation.
The current parser treats `command` as metadata, so it is safe to keep it in the config.

Run from `external_wop_cpp/` after building:

```bash
./build/wop_cli --config examples/box_wop.yaml --json
./build/wop_cli --config examples/box_wos.yaml --json
```

On Windows with the Visual Studio generator:

```powershell
.\build_config_work\Debug\wop_cli.exe --config examples/box_wop.yaml --json
.\build_config_work\Debug\wop_cli.exe --config examples/box_wos.yaml --json
```
