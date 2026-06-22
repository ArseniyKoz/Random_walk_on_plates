# Walk on Planes / Walk on Spheres

Проект про численное решение внешней задачи Дирихле для гармонической функции в трехмерной области вне выпуклого многогранника. Основная идея: сравнить и реализовать Monte Carlo подходы `Walk on Planes` и `Walk on Spheres`, довести их до C++/Python кода, CLI, тестов, notebooks и proof-oriented документации.

## English Summary

This repository implements and documents Monte Carlo solvers for an exterior Dirichlet problem in 3D: Walk on Planes and Walk on Spheres. The public surface includes a C++ core, a Python reference implementation, CLI experiments, validation tests, notebooks, and proof-oriented notes for the far-field handling modes.

## Что внутри

- `external_wop_cpp/` - основная C++ реализация: геометрия, sampling, WoP/WoS solvers, CLI и C++ tests.
- `external_wop/` - Python reference implementation и pytest tests для проверки геометрии, sampling и parity с C++.
- `docs/` - разбор постановки задачи, геометрии, распределений, алгоритмов, валидации и теорем.
- `compare_wop_cpp_rmax_modes.ipynb`, `compare_wop_cpp_python.ipynb`, `wop-vs-wos-cpp-unit-cube.ipynb` - notebooks для сравнения режимов и реализаций.

## Evidence

- **What this proves:** исследовательский алгоритм не оставлен в виде одного notebook; у него есть C++/Python реализации, тесты статистических инвариантов, CLI и документация по математике.
- **Where to verify:** [docs/index.md](docs/index.md), [docs/validation-and-tests.md](docs/validation-and-tests.md), [docs/wop-proof-escape-project.md](docs/wop-proof-escape-project.md), [external_wop_cpp/tests](external_wop_cpp/tests), [external_wop/tests](external_wop/tests), [compare_wop_cpp_rmax_modes.ipynb](compare_wop_cpp_rmax_modes.ipynb).
- **Limits:** это исследовательская ВКР-работа, а не production library API; часть экспериментов живет в notebooks, а численные результаты зависят от `n`, `seed` и выбранного режима дальнего поля.

## Документация

Основная карта документации:

- [Постановка и геометрия](docs/problem-and-geometry.md)
- [Сэмплинг и плотности](docs/sampling-and-densities.md)
- [Алгоритм WoP](docs/wop-method.md)
- [Алгоритм WoS](docs/wos-method.md)
- [Валидация и тесты](docs/validation-and-tests.md)
- [Строгие теоремы для WoP: `escape` и `project`](docs/wop-proof-escape-project.md)
- [WoP: постановка, алгоритм и теорема о конечном поглощении](docs/wop-problem.md)

## Режимы дальнего поля в WoP

C++ CLI поддерживает два режима обработки дальнего поля:

- `escape` - траектория обрывается при достижении `r_max`;
- `project` - используется weighted far-sphere correction без геометрической проекции точки на сферу.

В режиме `project`:

- `r_max` трактуется как базовый радиус `rho`;
- `r_max_factor` задает дальний порог `rho1 = r_max_factor * rho`;
- при `--r-max 0` радиус `rho` выбирается автоматически из геометрии.

Ключевая формула возврата при попадании на границу:

```text
u = u_inf + eta * (f(y) - u_inf)
```

Именно weighted-схема убирает систематический bias, который возникал у старой геометрической проекции.

## Сборка C++ CLI

```bash
cmake -S external_wop_cpp -B external_wop_cpp/build_compare_modes -DCMAKE_BUILD_TYPE=Release -DWOP_BUILD_TESTS=OFF
cmake --build external_wop_cpp/build_compare_modes --config Release
```

Типичные пути к бинарнику:

- Linux/macOS: `external_wop_cpp/build_compare_modes/wop_cli`
- Windows: `external_wop_cpp/build_compare_modes/Release/wop_cli.exe`

## Примеры запуска

### WoP: `escape`

```bash
external_wop_cpp/build_compare_modes/wop_cli \
  --method wop --example box --x0 "3 0 0" --n 50000 --seed 12345 \
  --r-max-mode escape --r-max 1000000 --json
```

### WoP: `project` с ручным радиусом

```bash
external_wop_cpp/build_compare_modes/wop_cli \
  --method wop --example box --x0 "3 0 0" --n 50000 --seed 12345 \
  --r-max-mode project --r-max 2.5 --r-max-factor 2.0 --json
```

### WoP: `project` с auto `r_max`

```bash
external_wop_cpp/build_compare_modes/wop_cli \
  --method wop --example box --x0 "3 0 0" --n 50000 --seed 12345 \
  --r-max-mode project --r-max 0 --r-max-factor 2.0 --json
```

### WoS для сопоставления

```bash
external_wop_cpp/build_compare_modes/wop_cli \
  --method wos --example box --x0 "3 0 0" --n 50000 --seed 12345 \
  --delta 1e-3 --rho-scale 1.0 --rho1-scale 2.0 --json
```

Для систематического сравнения режимов используй [compare_wop_cpp_rmax_modes.ipynb](compare_wop_cpp_rmax_modes.ipynb). В notebook собраны sanity checks, матрица запусков по `n` и `seed`, таблицы и графики для `escape` vs `project`.
