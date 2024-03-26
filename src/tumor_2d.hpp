#ifndef TUMOR_2D_HPP
#define TUMOR_2D_HPP

#include <galois/Timer.h>

#include "DiffusionMap.hpp"
#include "InitialCondition.hpp"
#include "Treatment.hpp"
#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"

namespace ads::problems
{

    class tumor_2d : public simulation_2d
    {
    private:
        using Base = simulation_2d;
        vector_type u, u_prev;

        output_manager<2> output;
        galois_executor executor{1};
        galois::StatTimer integration_timer{"integration"};
        const InitialCondition &initial_condition;
        const DiffusionMap &diffusion_map;
        const Treatment &treatment;
        lin::tensor<double, 4> D;
        int lastDay = -1;

    public:
        explicit tumor_2d(const config_2d &config,
                          const InitialCondition &initial_condition,
                          const DiffusionMap &diffusion_map,
                          const Treatment &treatment)
            : Base{config}, u{shape()}, u_prev{shape()},
              output{x.B, y.B, 200}, initial_condition{initial_condition},
              diffusion_map{diffusion_map}, treatment{treatment},
              D{{x.basis.elements, y.basis.elements, x.basis.quad_order + 1,
                 y.basis.quad_order + 1}} {}

        double init_state(double x, double y)
        {
            return initial_condition(x, y);
        };

    private:
        void solve(vector_type &v)
        {
            // lin::vector buf{{y.dofs()}};
            // compute_projection(buf, y.basis,
            //                    [](double y) { return std::sin(y * M_PI); });
            // for (int i = 0; i < y.dofs(); ++i) {
            //   v(0, i) = buf(i);
            // }
            Base::solve(v);
        }

        void prepare_matrices()
        {
            x.fix_left();
            Base::prepare_matrices();
        }

        void before() override
        {
            fill_diffusion_map();
            prepare_matrices();

            auto init = [this](double x, double y)
            {
                return init_state(x, y);
            };
            projection(u, init);
            solve(u);

            output.to_file(u, "init.data");
        }

        void fill_diffusion_map()
        {
            for (auto e : elements())
            {
                for (auto q : quad_points())
                {
                    auto x = point(e, q);
                    D(e[0], e[1], q[0], q[1]) = diffusion_map(x[0], x[1]);
                }
            }
        }

        void before_step(int /*iter*/, double /*t*/) override
        {
            using std::swap;
            swap(u, u_prev);
        }

        void step(int /*iter*/, double t) override
        {
            compute_rhs(t);
            solve(u);
        }

        void after_step(int /*iter*/, double t) override
        {
            if (t > lastDay)
            {
                lastDay++;
                output.to_file(u, "out/step_%d.data", lastDay);
            }
        }

        void compute_rhs(double t)
        {
            integration_timer.start();
            auto &rhs = u;

            zero(rhs);

            executor.for_each(elements(), [&](index_type e)
                              {
      auto U = element_rhs();

      double J = jacobian(e);
      for (auto q : quad_points()) {
        double w = weight(q);
        auto x = point(e, q);

        value_type u = eval_fun(u_prev, e, q);

        constexpr double rho = 0.025;
        double R = treatment(x[0], x[1], t);
        double D_val = D(e[0], e[1], q[0], q[1]);

        for (auto a : dofs_on_element(e)) {
          auto aa = dof_global_to_local(e, a);
          value_type v = eval_basis(e, q, a);

          double diffusive_term = -D_val * grad_dot(u, v);
          double growth_term = rho * u.val * (1 - u.val) * v.val;
          double inhibition_term = -R * u.val * v.val;

          double all_terms = diffusive_term + growth_term + inhibition_term;
          double val = u.val * v.val + steps.dt * all_terms;
          U(aa[0], aa[1]) += val * w * J;
        }
      }

      executor.synchronized([&]() { update_global_rhs(rhs, U, e); }); });
            integration_timer.stop();
        }

        void after() override
        {
            std::cout << "integration: "
                      << static_cast<double>(integration_timer.get())
                      << std::endl;
        }
    };

} // namespace ads::problems

#endif // TUMOR_2D_HPP
