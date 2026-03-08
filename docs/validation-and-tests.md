# Валидация И Тесты

Этот раздел фиксирует, как корректность WoP/WoS подтверждается в C++.

Основная test-база:

- `external_wop_cpp/tests/`

## 1. Базовые smoke-инварианты

Файл: `external_wop_cpp/tests/test_core_smoke.cpp`

Проверяет:

- базовую геометрию box (`inside/outside`, closest boundary),
- корректный ранний `escaped` в WoP при `r_max`-границе,
- корректность агрегатора оценок (`J`, `S2`, `eps`, `n_truncated`),
- базовый smoke для WoS.

## 2. Статистические тесты sampling слоя

Файл: `external_wop_cpp/tests/test_sampling_stats.cpp`

Проверяет:

1. `sample_tangent_direction`:
   - `||w|| = 1`,
   - `dot(w, nu) = 0`.
2. Второй момент касательного сэмплера:
   - эмпирическая матрица близка к теоретической `0.5 * (I - nu nu^T)`.
3. Квантили радиального распределения для попаданий в плоскость:
   - сравнение с аналитическим `F_R(r)`.

Это прямая проверка того, что «плотность моделируется правильно».

## 3. Поэтапные тесты WoP

Файл: `external_wop_cpp/tests/test_wop_pipeline_stages.cpp`

Проверяет:

- флаги и индексы в `scan_distances`,
- валидность/область определения `sample_plane_radius`,
- численную согласованность `estimate_wop` с точным гармоническим решением в box-задаче.

## 4. Режимы дальнего поля WoP (`escape/project`)

Файл: `external_wop_cpp/tests/test_wop_rmax_modes.cpp`

Проверяет:

- что `escape` корректно обрывает траектории,
- что `project` работает через far-sphere схему без систематического немедленного `escaped`,
- корректность конфигурации `rho/rho1` при auto/manual настройках.

## 5. Внешняя валидация WoS

Файл: `external_wop_cpp/tests/test_wos_lecture_external.cpp`

Проверяет:

- остановку траектории при попадании в `delta`-окрестность границы,
- near-boundary сценарий (возврат граничного значения),
- согласованность с гармоническим референсом `u(x)=1/||x-a||`,
- согласованность общего `estimate_wos(poly, ...)` и box-wrapper `estimate_wos_box(...)`.

## 6. Что означают метрики отчёта

См. `external_wop_cpp/include/wop/estimation/estimation.hpp`:

- `J` — итоговая оценка `u(x0)`,
- `S2` — эмпирическая дисперсия значений траекторий,
- `eps = 3 * sqrt(S2/n)` — оценка ошибки,
- `n_truncated` — число `timeout + escaped`,
- `mean_steps` — среднее число шагов.

## 7. Практический сценарий проверки через CLI

CLI: `external_wop_cpp/app/main.cpp`.

Пример WoP:

```powershell
external_wop_cpp/build/Release/wop_cli.exe --method wop --example box --x0 "3 0 0" --n 100000 --seed 314159 --json
```

Пример WoS:

```powershell
external_wop_cpp/build/Release/wop_cli.exe --method wos --example box --x0 "3 0 0" --n 100000 --seed 314159 --delta 1e-3 --rho-scale 1.0 --rho1-scale 2.0 --json
```

Для сравнения качества смотреть одновременно:

- `abs_error`,
- `eps`,
- `n_truncated`,
- `mean_steps`,
- wall-clock время.
