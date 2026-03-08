# Алгоритм WoS (Walk on Spheres) для Общего Многогранника

## 1. API и архитектура

Основной API:

- `trace_wos_trajectory(...)`
- `estimate_wos(...)`

Источник: `external_wop_cpp/include/wop/solver/wos_solver.hpp`.

`wos_box_solver` сейчас только обёртка, которая:

1. строит box как `Polyhedron`,
2. ориентирует нормали,
3. вызывает общий `estimate_wos(...)`.

Источник: `external_wop_cpp/src/solver/wos_box_solver.cpp`.

## 2. Предвычисляемый WoS-контекст

Перед трассировкой собирается `WosProjectionContext`:

- ссылка на `poly`,
- `center`, `rho`, `rho1` (сфера дальнего поля),
- список `active_sets` для проекции на границу.

Источник: `resolve_wos_projection_context(...)` в `external_wop_cpp/src/solver/wos_solver.cpp`.

### 2.1 Bounding sphere

1. Строятся вершины многогранника как пересечения троек плоскостей.
2. Отбрасываются вырожденные и невалидные точки.
3. Берётся AABB вершин, центр — середина AABB.
4. Радиус — максимум расстояний от центра до вершин.

Далее:

- `rho = rho_scale * base_radius`,
- `rho1 = rho1_scale * rho`, где `rho1_scale > 1`.

### 2.2 Active-set данные

Готовятся все наборы активных ограничений размером 1, 2, 3.
Для каждого набора предварительно инвертируется малая Gram-матрица нормалей.

Это ускоряет многократные проекции точек на `∂Q`.

## 3. Проекция точки на границу многогранника

`project_to_polyhedron_boundary(ctx, x)` перебирает active sets и решает KKT-подобную схему:

1. `rhs_i = dot(nu_i, x) - b_i`,
2. `lambda = inv(Gram) * rhs`,
3. проверка `lambda_i >= -lambda_tol`,
4. `y = x - sum_i lambda_i * nu_i`,
5. проверка допустимости `dot(nu_j, y)-b_j <= feas_tol` для всех `j`.

Среди допустимых `y` выбирается минимальное `||x-y||`.

Если ни один кандидат не прошёл проверки, выбрасывается исключение.

## 4. Один шаг WoS

Внутри `trace_wos_with_context(...)`:

1. Проверяется, что `x0` во внешней области.
2. Поддерживается вес `eta` (инициализация `1.0`).
3. На каждом шаге:
   - `proj = project_to_polyhedron_boundary(ctx, x)`,
   - если `proj.distance <= delta`, траектория завершена:
     `value = eta * boundary_f(proj.point, nullopt)`.

Если до границы ещё далеко:

- Near region (`r = ||x-center|| <= rho1`):
  - `omega ~ Unif(S^2)`,
  - `x <- x + proj.distance * omega`.

- Far region (`r > rho1`):
  - `x <- sample_far_sphere_step(x, center, rho, rng)`,
  - `eta <- eta * (rho / r)`.

Если достигнут `max_steps`, статус `timeout`.

## 5. Где здесь «попадание на окружность/сферу»

В WoS есть два разных события:

1. Попадание в `delta`-окрестность границы многогранника:
   - проверка через расстояние до ближайшей граничной точки.
2. Переход на сферу дальнего поля:
   - точка генерируется конструктивно на сфере радиуса `rho` (far-sphere formula),
   - отдельный тест равенства радиуса не нужен.

То есть в 3D главным критерием остаётся расстояние до границы (`delta`), а не проверка «равно радиусу до машинного нуля».

## 6. Псевдокод WoS

```text
prepare_context(poly, rho_scale, rho1_scale):
    center, base_rho <- bounding_sphere(poly)
    rho <- rho_scale * base_rho
    rho1 <- rho1_scale * rho
    active_sets <- build_active_sets(poly)
    return context

trace_wos(context, x0):
    x <- x0
    eta <- 1
    assert x0 in exterior

    for step = 1..max_steps:
        proj <- project_to_polyhedron_boundary(context, x)
        if proj.distance <= delta:
            return (eta * f(proj.point), hit_face)

        r <- ||x - center||
        if r <= rho1:
            omega <- uniform_unit_sphere()
            x <- x + proj.distance * omega
        else:
            x <- sample_far_sphere_step(x, center, rho)
            eta <- eta * (rho / r)

    return (u_inf, timeout)
```

## 7. Отличие WoS от WoP в этой реализации

- WoP: шаги по случайным плоскостям граней и строгая face-проверка через signed distances.
- WoS: шаги по сферам вокруг текущей точки с приближением через проекцию на `∂Q` и `delta`.
- Оба используют far-sphere weighted correction для устойчивого дальнего поля.

## 8. Оценка `estimate_wos`

Агрегация статистики полностью общая с WoP через `estimate_from_trajectories(...)`:

- среднее `J`,
- дисперсионная оценка `S2`,
- доверительная ширина `eps`,
- счётчики усечённых траекторий.
