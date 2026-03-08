# Сэмплинг И Плотности

## 1. Генератор случайных чисел

В C++ используется класс `wop::rng::Rng`:

- движок: `std::mt19937_64`,
- `uniform01()` -> `U(0,1)`,
- `normal01()` -> `N(0,1)`,
- `normal3()` -> 3 независимых `N(0,1)`.

Источник: `external_wop_cpp/include/wop/rng/rng.hpp`, `external_wop_cpp/src/rng/rng.cpp`.

## 2. Равномерное направление на сфере

Функция `sample_unit_sphere(rng)`:

1. берёт `w ~ N(0, I3)`,
2. возвращает `w / ||w||`.

Это даёт равномерное направление на `S^2` за счёт сферической симметрии гауссова распределения.

Источник: `external_wop_cpp/src/sampling/sampling.cpp`.

## 3. Тангенциальное направление к плоскости

Для нормали `nu`:

1. `g ~ N(0, I3)`,
2. `w_raw = g - dot(g, nu_unit) * nu_unit`,
3. `w = w_raw / ||w_raw||`.

Если `||w_raw|| <= min_norm` или значение нечисловое, идёт ресэмплинг (`max_resample` попыток).

Факт: `w` лежит в касательной плоскости (`dot(w, nu)=0`) и имеет единичную длину.

Источник:

- `sample_tangent_direction(...)` в `external_wop_cpp/src/sampling/sampling.cpp`,
- статистические проверки в `external_wop_cpp/tests/test_sampling_stats.cpp`.

## 4. Распределение радиуса для WoP шага по плоскости

В WoP после выбора активной плоскости используется:

- `alpha ~ U(0,1)`,
- `rho = d * sqrt(1/(1-alpha)^2 - 1)`.

Это inverse-CDF для радиуса `rho` в плоскости (Poisson-kernel форма):

- `F_R(r) = 1 - d / sqrt(d^2 + r^2)`,
- `p_R(r) = d * r / (d^2 + r^2)^(3/2)`.

Реализация:

- `sample_plane_radius(d, rng, min_one_minus_alpha)` в `external_wop_cpp/src/solver/wop_solver_internal.cpp`.

Численная защита:

- вместо `1-alpha` используется `max(1-alpha, min_one_minus_alpha)`.

## 5. Сэмплинг попадания луча в плоскость (вспомогательная функция)

`sample_hit_on_plane_from_point(...)`:

1. `omega ~ Unif(S^2)`,
2. принудительно выбирается полупространство к плоскости (`dot(nu, omega) <= 0`),
3. `t = -d / dot(nu, omega)`,
4. `y = x + t * omega`.

Если `|dot(nu, omega)| < min_abs_denom`, `t<=0`, или `y` не конечна, идёт ресэмплинг.

Источник: `external_wop_cpp/src/sampling/sampling.cpp`.

## 6. Как в коде проверяется «попадание на окружность / сферу»

В текущей 3D реализации проверка делается не через точное равенство, а через геометрию + допуски:

1. Для граней многогранника: signed distances (`eps_in`, `eps_plane`).
2. Для сферы дальнего поля:
   - условие выхода в far-ветку: `r = ||x-center|| > rho1`,
   - новая точка генерируется сразу на сфере радиуса `rho` (`||x_new-center|| = rho` конструктивно).

То есть проверка «попали на сферу» встроена в формулу построения точки, а не в отдельный `==`.

## 7. Far-sphere сэмплинг

И WoP `project`, и WoS используют один тип перехода:

1. `e = (x-center)/||x-center||`,
2. `alpha ~ U(0,1)`,
3. вычисляются `q`, `z1`,
4. берётся случайный `w`, ортогональный `e`,
5. направление `dir = z1*e + sqrt(1-z1^2)*w`,
6. точка `center + rho*dir`.

Источник:

- WoP: `sample_far_sphere_step(...)` в `external_wop_cpp/src/solver/wop_solver_internal.cpp`
- WoS: `sample_far_sphere_step(...)` в `external_wop_cpp/src/solver/wos_solver.cpp`.

## 8. Мини-псевдокод sampling слоя

```text
sample_unit_sphere():
    repeat:
        g <- N(0, I3)
        if ||g|| > 0: return g / ||g||

sample_tangent_direction(nu):
    nu <- normalize(nu)
    repeat up to max_resample:
        g <- N(0, I3)
        w <- g - dot(g, nu) * nu
        if ||w|| > min_norm: return w / ||w||
    fail

sample_plane_radius(d):
    alpha <- U(0,1)
    t <- max(1-alpha, min_one_minus_alpha)
    return d * sqrt(1/t^2 - 1)
```
