#include "wop/solver/wos_box_solver.hpp"

#include "wop/geometry/box.hpp"
#include "wop/solver/wos_solver.hpp"

namespace wop::solver {

estimation::TrajectoryResult trace_wos_box_trajectory(
    const math::Vec3& x0,
    const math::Vec3& box_min,
    const math::Vec3& box_max,
    const WosBoundaryFunc& boundary_f,
    rng::Rng& rng,
    double delta,
    double rho_scale,
    double rho1_scale,
    int max_steps,
    double u_inf) {
    auto poly = geometry::make_axis_aligned_box(box_min, box_max);
    const math::Vec3 center = 0.5 * (box_min + box_max);
    poly = geometry::orient_normals(poly, center);
    return trace_wos_trajectory(poly, x0, boundary_f, rng, delta, rho_scale, rho1_scale, max_steps, u_inf);
}

estimation::EstimateResult estimate_wos_box(
    const math::Vec3& x0,
    const math::Vec3& box_min,
    const math::Vec3& box_max,
    const WosBoundaryFunc& boundary_f,
    int n_paths,
    rng::Rng& rng,
    double delta,
    double rho_scale,
    double rho1_scale,
    int max_steps,
    double u_inf) {
    auto poly = geometry::make_axis_aligned_box(box_min, box_max);
    const math::Vec3 center = 0.5 * (box_min + box_max);
    poly = geometry::orient_normals(poly, center);
    return estimate_wos(poly, x0, boundary_f, n_paths, rng, delta, rho_scale, rho1_scale, max_steps, u_inf);
}

}  // namespace wop::solver

