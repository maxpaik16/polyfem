#include "NLProblem.hpp"

#include <polyfem/io/OBJWriter.hpp>
#include <fstream>
#include <polyfem/solver/forms/ContactForm.hpp>

/*
m \frac{\partial^2 u}{\partial t^2} = \psi = \text{div}(\sigma[u])\newline
u^{t+1} = u(t+\Delta t)\approx u(t) + \Delta t \dot u + \frac{\Delta t^2} 2 \ddot u \newline
= u(t) + \Delta t \dot u + \frac{\Delta t^2}{2} \psi\newline
M u^{t+1}_h \approx M u^t_h + \Delta t M v^t_h + \frac{\Delta t^2} {2} A u^{t+1}_h \newline
%
M (u^{t+1}_h - (u^t_h + \Delta t v^t_h)) - \frac{\Delta t^2} {2} A u^{t+1}_h
*/
// mü = ψ = div(σ[u])
// uᵗ⁺¹ = u(t + Δt) ≈ u(t) + Δtu̇ + ½Δt²ü = u(t) + Δtu̇ + ½Δt²ψ
// Muₕᵗ⁺¹ ≈ Muₕᵗ + ΔtMvₕᵗ ½Δt²Auₕᵗ⁺¹
// Root-finding form:
// M(uₕᵗ⁺¹ - (uₕᵗ + Δtvₕᵗ)) - ½Δt²Auₕᵗ⁺¹ = 0

namespace polyfem::solver
{
	NLProblem::NLProblem(
		const int full_size,
		const std::vector<int> &boundary_nodes,
		const std::vector<std::shared_ptr<Form>> &forms)
		: FullNLProblem(forms),
		  boundary_nodes_(boundary_nodes),
		  full_size_(full_size),
		  reduced_size_(full_size_ - boundary_nodes.size()),
		  t_(0),
		  rhs_assembler_(nullptr),
		  local_boundary_(nullptr),
		  n_boundary_samples_(0)
	{
		use_reduced_size();
	}

	NLProblem::NLProblem(
		const int full_size,
		const std::vector<int> &boundary_nodes,
		const std::vector<mesh::LocalBoundary> &local_boundary,
		const int n_boundary_samples,
		const assembler::RhsAssembler &rhs_assembler,
		const std::shared_ptr<utils::PeriodicBoundary> &periodic_bc,
		const double t,
		const std::vector<std::shared_ptr<Form>> &forms,
		const std::vector<std::set<int>> &neighbors,
		const bool use_neighbors_for_precond)
		: FullNLProblem(forms),
		  full_boundary_nodes_(boundary_nodes),
		  boundary_nodes_(periodic_bc ? periodic_bc->full_to_periodic(boundary_nodes) : boundary_nodes),
		  full_size_(full_size),
		  reduced_size_((periodic_bc ? periodic_bc->n_periodic_dof() : full_size) - boundary_nodes_.size()),
		  periodic_bc_(periodic_bc),
		  t_(t),
		  rhs_assembler_(&rhs_assembler),
		  local_boundary_(&local_boundary),
		  n_boundary_samples_(n_boundary_samples),
		  neighbors_(neighbors),
		  use_neighbors_for_precond_(use_neighbors_for_precond)
	{
		assert(std::is_sorted(boundary_nodes.begin(), boundary_nodes.end()));
		assert(boundary_nodes.size() == 0 || (boundary_nodes.front() >= 0 && boundary_nodes.back() < full_size_));
		use_reduced_size();
	}

	NLProblem::NLProblem(
		const int full_size,
		const std::vector<int> &boundary_nodes,
		const std::vector<mesh::LocalBoundary> &local_boundary,
		const int n_boundary_samples,
		const assembler::RhsAssembler &rhs_assembler,
		const std::shared_ptr<utils::PeriodicBoundary> &periodic_bc,
		const double t,
		const std::vector<std::shared_ptr<Form>> &forms)
		: FullNLProblem(forms),
		  full_boundary_nodes_(boundary_nodes),
		  boundary_nodes_(periodic_bc ? periodic_bc->full_to_periodic(boundary_nodes) : boundary_nodes),
		  full_size_(full_size),
		  reduced_size_((periodic_bc ? periodic_bc->n_periodic_dof() : full_size) - boundary_nodes_.size()),
		  periodic_bc_(periodic_bc),
		  t_(t),
		  rhs_assembler_(&rhs_assembler),
		  local_boundary_(&local_boundary),
		  n_boundary_samples_(n_boundary_samples)
	{
		assert(std::is_sorted(boundary_nodes.begin(), boundary_nodes.end()));
		assert(boundary_nodes.size() == 0 || (boundary_nodes.front() >= 0 && boundary_nodes.back() < full_size_));
		use_reduced_size();
	}

	void NLProblem::init_lagging(const TVector &x)
	{
		FullNLProblem::init_lagging(reduced_to_full(x));
	}

	void NLProblem::update_lagging(const TVector &x, const int iter_num)
	{
		FullNLProblem::update_lagging(reduced_to_full(x), iter_num);
	}

	void NLProblem::update_quantities(const double t, const TVector &x)
	{
		t_ = t;
		const TVector full = reduced_to_full(x);
		for (auto &f : forms_)
			f->update_quantities(t, full);
	}

	void NLProblem::line_search_begin(const TVector &x0, const TVector &x1)
	{
		FullNLProblem::line_search_begin(reduced_to_full(x0), reduced_to_full(x1));
	}

	double NLProblem::max_step_size(const TVector &x0, const TVector &x1)
	{
		return FullNLProblem::max_step_size(reduced_to_full(x0), reduced_to_full(x1));
	}

	bool NLProblem::is_step_valid(const TVector &x0, const TVector &x1)
	{
		return FullNLProblem::is_step_valid(reduced_to_full(x0), reduced_to_full(x1));
	}

	bool NLProblem::is_step_collision_free(const TVector &x0, const TVector &x1)
	{
		return FullNLProblem::is_step_collision_free(reduced_to_full(x0), reduced_to_full(x1));
	}

	double NLProblem::value(const TVector &x)
	{
		// TODO: removed fearure const bool only_elastic
		return FullNLProblem::value(reduced_to_full(x));
	}

	void NLProblem::gradient(const TVector &x, TVector &grad)
	{
		TVector full_grad;
		FullNLProblem::gradient(reduced_to_full(x), full_grad);
		grad = full_to_reduced_grad(full_grad);
	}

	void NLProblem::hessian(const TVector &x, THessian &hessian)
	{
		THessian full_hessian;
		FullNLProblem::hessian(reduced_to_full(x), full_hessian);

		full_hessian_to_reduced_hessian(full_hessian, hessian);
	}

	void NLProblem::get_problematic_indices(std::vector<std::set<int>> &bad_indices)
	{
		for (auto &f : forms_)
		{
			std::vector<std::set<int>> full_bad_indices;
			if (f->get_problematic_indices(full_bad_indices))
			{
				const int dim = 3;
				dof_to_func_mapping.clear();
				std::set<int> neighbors_to_add;
				if (use_neighbors_for_precond_)
				{
					//std::cout << "SIZE: " << neighbors_.size() << std::endl; 
					//std::cout << "FSIZE: " << full_size_ << std::endl; 
					for (int i = 0; i < full_size_; ++i)
					{
						if (full_bad_indices[0].count(i) > 0)
						{
							neighbors_to_add.insert(neighbors_[i].begin(), neighbors_[i].end());
						}
					}
					full_bad_indices[0].insert(neighbors_to_add.begin(), neighbors_to_add.end());
				}

				bad_indices.clear();
				bad_indices.resize(1);
				long j = 0;
				size_t k = 0;
				for (int i = 0; i < full_size_; ++i)
				{
					if (k < boundary_nodes_.size() && boundary_nodes_[k] == i)
					{
						++k;
						continue;
					}
					if (full_bad_indices[0].count(i) > 0)
					{
						bad_indices[0].insert(j);
						dof_to_func_mapping.push_back(i % dim);
					}
					++j; 
				}

				return;
			}
		}
	}

	void NLProblem::solution_changed(const TVector &newX)
	{
		FullNLProblem::solution_changed(reduced_to_full(newX));
	}

	void NLProblem::post_step(const polysolve::nonlinear::PostStepData &data)
	{
		FullNLProblem::post_step(polysolve::nonlinear::PostStepData(data.iter_num, data.solver_info, reduced_to_full(data.x), reduced_to_full(data.grad)));

		if (state_->args["output"]["advanced"]["save_selected_nodes"])
		{
			std::ofstream grad_file(state_->resolve_output_path(fmt::format("grad_{:05d}.txt", total_step)));
			grad_file << reduced_to_full(data.grad);
			grad_file.close();


			for (auto &f : forms_)
			{
				ContactForm* contact_form = dynamic_cast<ContactForm*>(f.get());
				if (contact_form)
				{
					Eigen::VectorXd contact_grad = contact_form->last_grad;
					if (contact_grad.size() > 0)
					{
						std::ofstream contact_grad_file(state_->resolve_output_path(fmt::format("contact_grad_{:05d}.txt", total_step)));
						contact_grad_file << reduced_to_full(contact_grad);
						contact_grad_file.close();
					}
				}
			}

		}

		// TODO: add me back
		if (state_->args["output"]["advanced"]["save_nl_solve_sequence"])
		{
		 	const Eigen::MatrixXd displacements = utils::unflatten(reduced_to_full(data.x), state_->mesh->dimension());
		 	
			Eigen::MatrixXi edges = state_->collision_mesh.edges();
			Eigen::MatrixXi faces = state_->collision_mesh.faces();

			for (int i = 0; i < edges.rows(); ++i)
			{
				edges(i, 0) = state_->collision_mesh.to_full_vertex_id(edges(i, 0));
				edges(i, 1) = state_->collision_mesh.to_full_vertex_id(edges(i, 1));
			}

			for (int i = 0; i < faces.rows(); ++i)
			{
				faces(i, 0) = state_->collision_mesh.to_full_vertex_id(faces(i, 0));
				faces(i, 1) = state_->collision_mesh.to_full_vertex_id(faces(i, 1));
				faces(i, 2) = state_->collision_mesh.to_full_vertex_id(faces(i, 2));
			}

			Eigen::MatrixXd vertices = state_->collision_mesh.displace_vertices(displacements);
			Eigen::MatrixXd full_vertices(vertices.rows(), vertices.cols());
			for (int i = 0; i < vertices.rows(); ++i)
			{
				full_vertices.row(state_->collision_mesh.to_full_vertex_id(i)) = vertices.row(i);
			}
			
			io::OBJWriter::write(
		 		state_->resolve_output_path(fmt::format("nonlinear_solve_iter{:05d}.obj", total_step)),
		 		full_vertices,
		 		edges, faces);
		}
		++total_step;
	}

	void NLProblem::set_apply_DBC(const TVector &x, const bool val)
	{
		TVector full = reduced_to_full(x);
		for (auto &form : forms_)
			form->set_apply_DBC(full, val);
	}

	NLProblem::TVector NLProblem::full_to_reduced(const TVector &full) const
	{
		TVector reduced;
		full_to_reduced_aux(boundary_nodes_, full_size(), current_size(), full, reduced);
		return reduced;
	}

	NLProblem::TVector NLProblem::full_to_reduced_grad(const TVector &full) const
	{
		TVector reduced;
		full_to_reduced_aux_grad(boundary_nodes_, full_size(), current_size(), full, reduced);
		return reduced;
	}

	NLProblem::TVector NLProblem::reduced_to_full(const TVector &reduced) const
	{
		TVector full;
		reduced_to_full_aux(boundary_nodes_, full_size(), current_size(), reduced, boundary_values(), full);
		return full;
	}

	Eigen::MatrixXd NLProblem::boundary_values() const
	{
		Eigen::MatrixXd result = Eigen::MatrixXd::Zero(full_size(), 1);
		// rhs_assembler->set_bc(*local_boundary_, boundary_nodes_, n_boundary_samples_, local_neumann_boundary_, result, t_);
		rhs_assembler_->set_bc(*local_boundary_, full_boundary_nodes_, n_boundary_samples_, std::vector<mesh::LocalBoundary>(), result, Eigen::MatrixXd(), t_);
		return result;
	}

	template <class FullMat, class ReducedMat>
	void NLProblem::full_to_reduced_aux(const std::vector<int> &boundary_nodes, const int full_size, const int reduced_size, const FullMat &full, ReducedMat &reduced) const
	{
		using namespace polyfem;

		// Reduced is already at the full size
		if (full_size == reduced_size || full.size() == reduced_size)
		{
			reduced = full;
			return;
		}

		assert(full.size() == full_size);
		assert(full.cols() == 1);
		reduced.resize(reduced_size, 1);

		Eigen::MatrixXd mid;
		if (periodic_bc_)
			mid = periodic_bc_->full_to_periodic(full, false);
		else
			mid = full;
		
		assert(std::is_sorted(boundary_nodes.begin(), boundary_nodes.end()));

		long j = 0;
		size_t k = 0;
		for (int i = 0; i < mid.size(); ++i)
		{
			if (k < boundary_nodes.size() && boundary_nodes[k] == i)
			{
				++k;
				continue;
			}

			assert(j < reduced.size());
			reduced(j++) = mid(i);
		}
	}

	template <class ReducedMat, class FullMat>
	void NLProblem::reduced_to_full_aux(const std::vector<int> &boundary_nodes, const int full_size, const int reduced_size, const ReducedMat &reduced, const Eigen::MatrixXd &rhs, FullMat &full) const
	{
		using namespace polyfem;

		// Full is already at the reduced size
		if (full_size == reduced_size || full_size == reduced.size())
		{
			full = reduced;
			return;
		}

		assert(reduced.size() == reduced_size);
		assert(reduced.cols() == 1);
		full.resize(full_size, 1);

		assert(std::is_sorted(boundary_nodes.begin(), boundary_nodes.end()));

		long j = 0;
		size_t k = 0;
		Eigen::MatrixXd mid(reduced_size + boundary_nodes.size(), 1);
		for (int i = 0; i < mid.size(); ++i)
		{
			if (k < boundary_nodes.size() && boundary_nodes[k] == i)
			{
				++k;
				mid(i) = rhs(i);
				continue;
			}

			mid(i) = reduced(j++);
		}
		
		full = periodic_bc_ ? periodic_bc_->periodic_to_full(full_size, mid) : mid;
	}

	template <class FullMat, class ReducedMat>
	void NLProblem::full_to_reduced_aux_grad(const std::vector<int> &boundary_nodes, const int full_size, const int reduced_size, const FullMat &full, ReducedMat &reduced) const
	{
		using namespace polyfem;

		// Reduced is already at the full size
		if (full_size == reduced_size || full.size() == reduced_size)
		{
			reduced = full;
			return;
		}

		assert(full.size() == full_size);
		assert(full.cols() == 1);
		reduced.resize(reduced_size, 1);

		Eigen::MatrixXd mid;
		if (periodic_bc_)
			mid = periodic_bc_->full_to_periodic(full, true);
		else
			mid = full;

		long j = 0;
		size_t k = 0;
		for (int i = 0; i < mid.size(); ++i)
		{
			if (k < boundary_nodes.size() && boundary_nodes[k] == i)
			{
				++k;
				continue;
			}

			reduced(j++) = mid(i);
		}
	}

	void NLProblem::full_hessian_to_reduced_hessian(const THessian &full, THessian &reduced) const
	{
		// POLYFEM_SCOPED_TIMER("\tfull hessian to reduced hessian");
		THessian mid = full;
		
		if (periodic_bc_)
			periodic_bc_->full_to_periodic(mid);

		if (current_size() < full_size())
			utils::full_to_reduced_matrix(mid.rows(), mid.rows() - boundary_nodes_.size(), boundary_nodes_, mid, reduced);
		else
			reduced = mid;
	}
} // namespace polyfem::solver
