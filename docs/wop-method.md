# Алгоритм WoP (Walk on Planes)

## 1. API и входные параметры

Основной API:

- `trace_wop_trajectory(...)`
- `estimate_wop(...)`

Источник: `external_wop_cpp/include/wop/solver/wop_solver.hpp`.

Ключевые параметры:

- `x0` — стартовая точка во внешней области,
- `boundary_f(y, face_idx)` — граничное значение,
- `eps_in`, `eps_plane` — допуски для принадлежности грани/области,
- `max_steps` — лимит шагов,
- `u_inf` — значение на бесконечности,
- `r_max`, `r_max_mode`, `r_max_factor` — управление дальним полем.

## 2. Инициализация траектории

В начале:

1. валидируются параметры (`max_steps > 0`, допуски неотрицательные и т.д.),
2. вычисляется начальная активная внешняя грань через `closest_outside_face_index`.

Если стартовая точка не снаружи, но внутри/на границе, траектория завершается сразу как `hit_face` с `steps = 0`.

Источник: `external_wop_cpp/src/solver/wop_solver.cpp`.

## 3. Одна итерация WoP: геометрический смысл

Пусть текущая точка `x` и активная грань `i`.

1. `d_i = dot(nu_i, x) - b_i`.
2. Проекция на плоскость грани: `p = x - d_i * nu_i`.
3. Сэмплинг касательного единичного `w`.
4. Сэмплинг радиуса `rho` по inverse-CDF.
5. Кандидат на плоскости: `y = p + rho * w`.
6. Стабилизация на плоскости:
   `y = y - (dot(nu_i, y) - b_i) * nu_i`.
7. Проверка попадания на грань через `scan_distances(...)`.

Если попадание подтверждено, возвращается `boundary_f(y, i)` с весом `eta`.
Если нет, выбирается следующая внешняя грань и цикл продолжается из `x = y`.

## 4. Как проверяется попадание на грань

`scan_distances(values, eps_in, target_idx, eps_plane)` вычисляет:

- есть ли внешние ограничения (`any_outside`),
- индекс ближайшей нарушенной грани (`argmin_outside`),
- индекс по минимальному `|d_i|` (`argmin_abs`),
- флаг `hit_target_face`.

`hit_target_face` истинно, когда:

- `|d_target| <= eps_plane`,
- одновременно выполнено `all_inside`.

Это эквивалентно допусковой версии `y in F_target`.

Источник: `external_wop_cpp/src/solver/wop_solver_internal.cpp`.

## 5. Обработка дальнего поля (`r_max`)

### 5.1 Режим `Escape`

Если задан `r_max`, то при `||x||^2 >= r_max^2` траектория завершается как `escaped` с возвратом `u_inf`.

### 5.2 Режим `Project` (weighted far-sphere)

Вычисляется конфиг `RMaxProjectionConfig`:

- `center` и базовый `rho`,
- дальний порог `rho1 = r_max_factor * rho`.

Если `r = ||x-center|| > rho1`, выполняется far-sphere шаг и вес:

- `x <- sample_far_sphere_step(x, center, rho, rng)`,
- `eta <- eta * (rho / r)`.

Итоговое значение на попадании:

`u = u_inf + eta * (f(y) - u_inf)`.

Источники:

- `external_wop_cpp/src/solver/wop_solver.cpp`
- `external_wop_cpp/src/solver/wop_solver_internal.cpp`.

## 6. Статусы завершения траектории

`TrajectoryResult.status` принимает:

- `hit_face` — достижение границы,
- `escaped` — выход по правилам дальнего поля/численной деградации,
- `timeout` — превышен `max_steps` или возникла неустойчивая ветка сэмплинга.

## 7. Псевдокод WoP

```text
trace_wop(poly, x0):
    init x <- x0, eta <- 1
    init face_idx <- closest_outside_face_index(x)

    for step = 1..max_steps:
        if not finite(x, eta): return (u_inf, escaped)

        if r_max_mode == project and r(x, center) > rho1:
            x <- sample_far_sphere_step(x, center, rho)
            eta <- eta * (rho / r_old)
            update face_idx
            continue

        if r_max_mode == escape and ||x|| >= r_max:
            return (u_inf, escaped)

        d <- signed_distance_to_active_face(x, face_idx)
        if d <= eps_in:
            rescan distances
            if no outside: return (weighted f(x), hit_face)
            else face_idx <- closest outside; continue

        w <- sample_tangent_direction(nu_face)
        rho <- sample_plane_radius(d)
        y <- project_to_plane(p + rho * w)

        scan <- scan_distances(y)
        if scan.hit_target_face or no outside:
            return (weighted f(y), hit_face)

        x <- y
        face_idx <- scan.argmin_outside

    return (u_inf, timeout)
```

## 8. Что агрегирует `estimate_wop`

`estimate_wop(...)` многократно вызывает `trace_wop_trajectory(...)` и считает:

- `J = mean(value)`,
- `S2 = mean(value^2) - J^2`,
- `eps = 3 * sqrt(S2 / n_total)`,
- `n_truncated = n_timeout + n_escaped`,
- `mean_steps`.

Источник: `external_wop_cpp/src/estimation/estimation.cpp`.
