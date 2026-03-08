# Постановка И Геометрия

## 1. Математическая постановка задачи

Рассматривается внешняя задача Дирихле для гармонической функции:

- `Δu = 0` в `D = R^3 \ overline(Q)`,
- `u = f` на `Γ = ∂Q`,
- `u(x) -> u_inf` при `||x|| -> infinity`.

Вычисляется приближение `u(x0)` для стартовой точки `x0` во внешней области.

Практически это реализовано в:

- `external_wop_cpp/src/solver/wop_solver.cpp`
- `external_wop_cpp/src/solver/wos_solver.cpp`

## 2. Как задаётся пространство (многогранник)

### 2.1 Представление через плоскости

Каждая грань задаётся плоскостью:

`dot(nu_i, x) = b_i`, где `||nu_i|| = 1`.

Многогранник:

`Q = intersection_i { x : dot(nu_i, x) - b_i <= 0 }`.

В коде это хранится как:

- `std::vector<Vec3> nu_`
- `std::vector<double> b_`

Источник: `external_wop_cpp/include/wop/geometry/polyhedron.hpp`.

### 2.2 Построение из `Plane{p, nu}`

При построении из плоскостей:

1. нормаль нормируется,
2. вычисляется `b_i = dot(nu_i, p_i)`.

Источник: `build_polyhedron_from_planes(...)` в `external_wop_cpp/src/geometry/polyhedron.cpp`.

## 3. Signed distance и принадлежность области

Для точки `x`:

`d_i(x) = dot(nu_i, x) - b_i`.

Интерпретация:

- `d_i(x) <= 0` — точка в полупространстве грани,
- `d_i(x) > 0` — точка снаружи относительно i-й грани.

Точка считается внутри/на границе многогранника, если:

`d_i(x) <= eps_in` для всех `i`.

Это реализовано в:

- `Polyhedron::signed_distances_inplace(...)`
- `Polyhedron::is_inside_or_on(...)`.

## 4. Ориентация нормалей

Чтобы все неравенства имели единый смысл, нормали разворачиваются относительно внутренней точки `interior_point`:

- если `dot(nu_i, interior_point) - b_i > 0`, меняются знаки у `nu_i` и `b_i`.

Псевдокод:

```text
for each face i:
    if dot(nu[i], interior_point) - b[i] > 0:
        nu[i] = -nu[i]
        b[i] = -b[i]
```

Источник: `orient_normals(...)` в `external_wop_cpp/src/geometry/polyhedron.cpp`.

## 5. Выбор активной внешней грани

Для внешней точки `x` WoP выбирает грань с минимальным положительным расстоянием:

`argmin_{i: d_i(x) > eps_in} d_i(x)`.

В коде:

- `Polyhedron::closest_outside_face_index(...)`.

Если такой грани нет, выбрасывается исключение и далее solver проверяет, не находится ли точка уже на границе/внутри (это обработка near-boundary случая).

## 6. Что значит «попали на границу»

Для точки `y` и целевой грани `i` в WoP используется допусковая версия условия `y in F_i`:

- `|d_i(y)| <= eps_plane`,
- `d_j(y) <= eps_in` для всех `j`.

Это проверяется через `scan_distances(...)` (см. `external_wop_cpp/src/solver/wop_solver_internal.cpp`).

В WoS проверка другая: сначала считается ближайшая граничная точка `proj.point` и расстояние `proj.distance`; остановка, если:

- `proj.distance <= delta`.

## 7. Специальный случай Box и общий случай

`make_axis_aligned_box(min_corner, max_corner)` строит 6 плоскостей, но после этого box трактуется как обычный `Polyhedron`.

Источники:

- `external_wop_cpp/src/geometry/box.cpp`
- `external_wop_cpp/src/solver/wos_box_solver.cpp`

Сейчас `wos_box_solver` является thin-wrapper вокруг общего `estimate_wos(...)`:

1. строит box как polyhedron,
2. ориентирует нормали,
3. вызывает общий WoS solver.

## 8. Характеристический масштаб и допуски

`Polyhedron::characteristic_length()` берёт `max(1, max_i |b_i|)` и используется для дефолтов допусков в WoP:

- `eps_in = 1e-12 * L`
- `eps_plane = 1e-12 * L`.

Это делает допуски масштабирующимися с размером геометрии.
