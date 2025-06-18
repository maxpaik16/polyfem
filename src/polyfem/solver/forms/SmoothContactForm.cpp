#include "SmoothContactForm.hpp"
#include <polyfem/utils/Logger.hpp>
#include <polyfem/utils/Types.hpp>
#include <polyfem/utils/Timer.hpp>
#include <polyfem/utils/MatrixUtils.hpp>
#include <polyfem/utils/MaybeParallelFor.hpp>
#include <polyfem/io/OBJWriter.hpp>

#include <ipc/utils/eigen_ext.hpp>
#include <ipc/barrier/adaptive_stiffness.hpp>
#include <ipc/utils/world_bbox_diagonal_length.hpp>

namespace polyfem::solver
{
	template <int _dim>
	SmoothContactForm<_dim>::SmoothContactForm(const ipc::CollisionMesh &collision_mesh,
											   const double dhat,
											   const double avg_mass,
											   const double alpha_t,
											   const double alpha_n,
											   const bool use_adaptive_dhat,
											   const double min_distance_ratio,
											   const bool use_adaptive_barrier_stiffness,
											   const bool is_time_dependent,
											   const bool enable_shape_derivatives,
											   const ipc::BroadPhaseMethod broad_phase_method,
											   const double ccd_tolerance,
											   const int ccd_max_iterations) : ContactForm(collision_mesh, dhat, avg_mass, use_adaptive_barrier_stiffness, is_time_dependent, enable_shape_derivatives, broad_phase_method, ccd_tolerance, ccd_max_iterations), params(dhat_, alpha_t, 0, alpha_n, 0, _dim == 3 ? 2 : 1), use_adaptive_dhat(use_adaptive_dhat)
	{
		collision_set_ = std::make_shared<ipc::SmoothCollisions<_dim>>();
		collision_set_->set_are_shape_derivatives_enabled(enable_shape_derivatives);

		contact_potential_ = std::make_shared<ipc::SmoothContactPotential<ipc::SmoothCollisions<_dim>>>(params);
		params.set_adaptive_dhat_ratio(min_distance_ratio);
		if (use_adaptive_dhat)
		{
			get_smooth_collision_set().compute_adaptive_dhat(collision_mesh, collision_mesh.rest_positions(), params, broad_phase_);
			if (use_adaptive_barrier_stiffness)
				logger().error("Adaptive dhat is not compatible with adaptive barrier stiffness");
		}
	}

	template <int _dim>
	void SmoothContactForm<_dim>::update_barrier_stiffness(const Eigen::VectorXd &x, const Eigen::MatrixXd &grad_energy)
	{
		if (!use_adaptive_barrier_stiffness())
			return;

		log_and_throw_error("Adaptive barrier stiffness not implemented for SmoothContactForm!");
	}

	template <int _dim>
	void SmoothContactForm<_dim>::force_shape_derivative(ipc::CollisionsBase *collision_set, const Eigen::MatrixXd &solution, const Eigen::VectorXd &adjoint_sol, Eigen::VectorXd &term)
	{
		if (!collision_set)
		{
			term.setZero(solution.size());
			return;
		}
		StiffnessMatrix hessian = contact_potential_->hessian(*dynamic_cast<ipc::SmoothCollisions<_dim> *>(collision_set), collision_mesh_, compute_displaced_surface(solution), ipc::PSDProjectionMethod::NONE);
		term = barrier_stiffness() * collision_mesh_.to_full_dof(hessian) * adjoint_sol;
	}

	template <int _dim>
	void SmoothContactForm<_dim>::update_collision_set(const Eigen::MatrixXd &displaced_surface)
	{
		// Store the previous value used to compute the constraint set to avoid duplicate computation.
		static Eigen::MatrixXd cached_displaced_surface;
		if (cached_displaced_surface.size() == displaced_surface.size() && cached_displaced_surface == displaced_surface)
			return;

		if (use_cached_candidates_)
			get_smooth_collision_set().build(
				candidates_, collision_mesh_, displaced_surface, params, use_adaptive_dhat);
		else
			get_smooth_collision_set().build(
				collision_mesh_, displaced_surface, params, use_adaptive_dhat, broad_phase_);
		cached_displaced_surface = displaced_surface;
	}

	template <int _dim>
	double SmoothContactForm<_dim>::value_unweighted(const Eigen::VectorXd &x) const
	{
		return (*contact_potential_)(get_smooth_collision_set(), collision_mesh_, compute_displaced_surface(x));
	}

	template <int _dim>
	Eigen::VectorXd SmoothContactForm<_dim>::value_per_element_unweighted(const Eigen::VectorXd &x) const
	{
		log_and_throw_error("value_per_element_unweighted not implemented!");
	}

	template <int _dim>
	void SmoothContactForm<_dim>::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const
	{
		gradv = contact_potential_->gradient(get_smooth_collision_set(), collision_mesh_, compute_displaced_surface(x));
		gradv = collision_mesh_.to_full_dof(gradv);
	}

	template <int _dim>
	void SmoothContactForm<_dim>::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const
	{
		// {
		// 	io::OBJWriter::write(
		// 		"collision_mesh.obj",
		// 		compute_displaced_surface(x),
		// 		collision_mesh_.edges(), collision_mesh_.faces());
		// }
		POLYFEM_SCOPED_TIMER("barrier hessian");
		hessian = contact_potential_->hessian(get_smooth_collision_set(), collision_mesh_, compute_displaced_surface(x), project_to_psd_ ? ipc::PSDProjectionMethod::CLAMP : ipc::PSDProjectionMethod::NONE);
		hessian = collision_mesh_.to_full_dof(hessian);
	}

	template <int _dim>
	void SmoothContactForm<_dim>::post_step(const polysolve::nonlinear::PostStepData &data)
	{
		const Eigen::MatrixXd displaced_surface = compute_displaced_surface(data.x);

		const double curr_distance = get_smooth_collision_set().compute_minimum_distance(collision_mesh_, displaced_surface);
		const double curr_active_distance = get_smooth_collision_set().compute_active_minimum_distance(collision_mesh_, displaced_surface);
		if (!std::isinf(curr_distance))
		{
			const double ratio = sqrt(curr_distance) / dhat();
			const auto log_level = (ratio < 1e-6) ? spdlog::level::err : ((ratio < 1e-4) ? spdlog::level::warn : spdlog::level::debug);
			polyfem::logger().log(log_level, "Minimum distance during solve: {}, active distance: {}, dhat: {}", sqrt(curr_distance), sqrt(curr_active_distance), dhat());
		}

		if (data.iter_num == 0)
			return;

		if (use_adaptive_barrier_stiffness_)
		{
			if (is_time_dependent_)
			{
				const double prev_barrier_stiffness = barrier_stiffness();

				barrier_stiffness_ = ipc::update_barrier_stiffness(
					prev_distance_, curr_distance, max_barrier_stiffness_,
					barrier_stiffness(), ipc::world_bbox_diagonal_length(displaced_surface), 1e-7);

				if (barrier_stiffness() != prev_barrier_stiffness)
				{
					polyfem::logger().debug(
						"updated barrier stiffness from {:g} to {:g}",
						prev_barrier_stiffness, barrier_stiffness());
				}
			}
			else
			{
				// TODO: missing feature
				// update_barrier_stiffness(data.x);
			}
		}

		prev_distance_ = curr_distance;
	}

	template class SmoothContactForm<2>;
	template class SmoothContactForm<3>;
} // namespace polyfem::solver