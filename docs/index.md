# WoP/WoS в C++: Детальная Документация

Этот раздел описывает, как в текущей C++-реализации устроены методы:

- `WoP` (Walk on Planes),
- `WoS` (Walk on Spheres),
- общая геометрия, сэмплинг, статистическая оценка и валидация.

Основной источник истины: `external_wop_cpp/`.

Ключевые точки входа:

- CLI: `external_wop_cpp/app/main.cpp`
- WoP solver: `external_wop_cpp/src/solver/wop_solver.cpp`
- WoS solver: `external_wop_cpp/src/solver/wos_solver.cpp`
- Геометрия: `external_wop_cpp/src/geometry/polyhedron.cpp`
- Сэмплинг: `external_wop_cpp/src/sampling/sampling.cpp`
- Оценка Монте-Карло: `external_wop_cpp/src/estimation/estimation.cpp`

## Что именно покрыто

В документации разобраны:

- постановка внешней задачи Дирихле,
- задание пространства через полупространства,
- signed distance и критерии попадания на грань/границу,
- как моделируются распределения (плотности) для переходов,
- как работают `r_max`-режимы в WoP,
- как устроена проекция на границу и far-sphere шаг в WoS,
- какие инварианты проверяются тестами.

## Карта документов

1. [Постановка И Геометрия](./problem-and-geometry.md)
2. [Сэмплинг И Плотности](./sampling-and-densities.md)
3. [Алгоритм WoP](./wop-method.md)
4. [Алгоритм WoS](./wos-method.md)
5. [Валидация И Тесты](./validation-and-tests.md)
6. [Строгие Теоремы Для WoP: `escape` И `project`](./wop-proof-escape-project.md)
7. [WoP: Постановка, Алгоритм И Теорема О Конечном Поглощении](wop-problem.md)

## Термины и обозначения

- `Q` — выпуклый многогранник.
- `Γ = ∂Q` — его граница.
- `D = R^3 \ overline(Q)` — внешняя область, где трассируются траектории.
- `nu_i`, `b_i` — нормаль и смещение плоскости `dot(nu_i, x) = b_i`.
- `d_i(x) = dot(nu_i, x) - b_i` — signed distance к i-й плоскости.
- `eps_in`, `eps_plane`, `delta` — численные допуски.

## Короткая выдержка из локальной теории

Из `external_wop/docs/walk_on_planes_dirichlet_R3.md`:

> `Q = intersect_j { x : d_j(x) <= 0 }` и  
> `F_i = { x : d_i(x)=0 and d_j(x)<=0 for all j }`.

Это напрямую объясняет, почему в коде проверка попадания на грань делается через условия на signed distances.
