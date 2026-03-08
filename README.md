# WOP C++: `r_max` Escape vs Project

Этот файл описывает, как в C++ `wop_cli` работают режимы обработки дальнего поля:

- `escape`: траектория обрывается при достижении `r_max`;
- `project`: используется weighted far-sphere correction (без геометрической проекции точки на сферу).

## Полная документация WoP/WoS

Детальное описание общей C++-реализации (геометрия, плотности/сэмплинг, WoP, WoS, тесты) теперь вынесено в:

- `docs/index.md`

## Как работает `project` (детально)

### 1) Что задаёт CLI

В `external_wop_cpp/app/main.cpp`:

- `--r-max-mode` принимает `escape` или `project`;
- `--r-max` парсится так: если значение `<= 0`, то это `std::nullopt`;
- `--r-max-factor` должен быть `> 1.0`.

В режиме `project`:

- `r_max` трактуется как базовый радиус `rho`;
- `r_max_factor` трактуется как множитель для дальнего порога `rho1 = r_max_factor * rho`.

При `--r-max 0` (то есть `nullopt`) `rho` выбирается автоматически из геометрии.

### 2) Что происходит внутри solver

В `trace_wop_trajectory(...)` (`external_wop_cpp/src/solver/wop_solver.cpp`):

1. Для `escape` используется legacy-проверка `||x||^2 >= r_max^2 -> escaped`.
2. Для `project` строится конфиг (`center`, `rho`, `rho1`).
3. На шаге траектории:
   - если `r = ||x-center|| > rho1`, выполняется far-sphere jump:
     `x <- sample_far_sphere_step(x, center, rho, rng)`;
   - обновляется вес `eta <- eta * (rho / r)`;
   - обычный WOP-plane шаг продолжается уже из новой точки.
4. При попадании на границу возвращается взвешенное значение:

```text
u = u_inf + eta * (f(y) - u_inf)
```

Именно эта формула убирает систематический bias, который был у старой геометрической проекции.

### 3) Как выбираются `rho` и `rho1`

В `resolve_r_max_projection(...)` (`external_wop_cpp/src/solver/wop_solver_internal.cpp`):

- если `r_max` задан явно: `rho = r_max`, `rho1 = r_max_factor * rho`;
- если `r_max` не задан:
  - `rho` стартует от `r_max_factor * poly_radius`,
  - затем поднимается минимум до `1.01 * dist(x0, center)` для устойчивого старта,
  - `rho1 = r_max_factor * rho`.

## Почему это корректнее старого `project`

Старая схема «жёстко проецировать точку на сферу и идти дальше» меняла задачу и давала устойчивый сдвиг оценки.
Новая схема использует вероятностную коррекцию с весом (как в far-field WoS-подходах), поэтому `J` остаётся согласованным с `exact` при росте `n`.

## Сборка `wop_cli`

```powershell
cmake -S external_wop_cpp -B external_wop_cpp/build_compare_modes -DCMAKE_BUILD_TYPE=Release -DWOP_BUILD_TESTS=OFF
cmake --build external_wop_cpp/build_compare_modes --config Release
```

Бинарник (Windows): `external_wop_cpp/build_compare_modes/Release/wop_cli.exe`

## Примеры запуска

### 1) Стандартный режим `escape` (обрыв)

```powershell
external_wop_cpp/build_compare_modes/Release/wop_cli.exe `
  --method wop --example box --x0 "3 0 0" --n 50000 --seed 12345 `
  --r-max-mode escape --r-max 1000000 --json
```

### 2) `project` с фиксированным `r_max` (ручной радиус)

```powershell
external_wop_cpp/build_compare_modes/Release/wop_cli.exe `
  --method wop --example box --x0 "3 0 0" --n 50000 --seed 12345 `
  --r-max-mode project --r-max 2.5 --r-max-factor 2.0 --json
```

Здесь `r_max` это `rho`, а `r_max-factor` задаёт `rho1/rho`.

### 3) `project` с auto `r_max`

```powershell
external_wop_cpp/build_compare_modes/Release/wop_cli.exe `
  --method wop --example box --x0 "3 0 0" --n 50000 --seed 12345 `
  --r-max-mode project --r-max 0 --r-max-factor 2.0 --json

### 4) WoS режим для сопоставления с WOP

```powershell
external_wop_cpp/build_compare_modes/Release/wop_cli.exe `
  --method wos --example box --x0 "3 0 0" --n 50000 --seed 12345 `
  --delta 1e-3 --rho-scale 1.0 --rho1-scale 2.0 --json
```
```

Ключевой момент: `--r-max 0` => `nullopt` => `rho` вычисляется автоматически.

### 4) Быстрое сравнение `escape` vs `project` в PowerShell

```powershell
$cli = "external_wop_cpp/build_compare_modes/Release/wop_cli.exe"

$escape = & $cli --method wop --example box --x0 "3 0 0" --n 100000 --seed 314159 `
  --r-max-mode escape --r-max 1000000 --json | ConvertFrom-Json

$project = & $cli --method wop --example box --x0 "3 0 0" --n 100000 --seed 314159 `
  --r-max-mode project --r-max 0 --r-max-factor 2.0 --json | ConvertFrom-Json

[pscustomobject]@{
  J_escape = $escape.J
  J_project = $project.J
  abs_error_escape = $escape.abs_error
  abs_error_project = $project.abs_error
  n_truncated_escape = $escape.n_truncated
  n_truncated_project = $project.n_truncated
  runtime_note = "runtime измеряй внешним таймером (Measure-Command)"
}
```

## Готовый notebook для сравнения

Для систематического сравнения запусти:

- `compare_wop_cpp_rmax_modes.ipynb`

Там уже есть:

- sanity checks (tests-first),
- матрица запусков по `n` и `seed`,
- таблицы и графики для `escape` vs `project`.
