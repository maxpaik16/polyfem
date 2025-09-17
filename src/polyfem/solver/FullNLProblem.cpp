#include "FullNLProblem.hpp"
#include <polyfem/utils/Logger.hpp>

#include <polyfem/solver/forms/ElasticForm.hpp>
#include <polyfem/solver/forms/ContactForm.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace polyfem::solver
{
	FullNLProblem::FullNLProblem(const std::vector<std::shared_ptr<Form>> &forms)
		: forms_(forms)
	{
	}

	double FullNLProblem::normalize_forms()
	{
		double total_weight = 0;
		for (const auto &f : forms_)
			total_weight += f->weight();

		logger().debug("Normalizing forms with scale: {}", total_weight);

		for (auto &f : forms_)
			f->set_scale(total_weight);

		return total_weight;
	}

	void FullNLProblem::init(const TVector &x)
	{
		for (auto &f : forms_)
			f->init(x);
	}

	void FullNLProblem::set_project_to_psd(bool project_to_psd)
	{
		project_to_psd_ = project_to_psd;
		for (auto &f : forms_)
			f->set_project_to_psd(project_to_psd);
	}

	void FullNLProblem::init_lagging(const TVector &x)
	{
		for (auto &f : forms_)
			f->init_lagging(x);
	}

	void FullNLProblem::update_lagging(const TVector &x, const int iter_num)
	{
		for (auto &f : forms_)
			f->update_lagging(x, iter_num);
	}

	int FullNLProblem::max_lagging_iterations() const
	{
		int max_lagging_iterations = 1;
		for (auto &f : forms_)
			max_lagging_iterations = std::max(max_lagging_iterations, f->max_lagging_iterations());
		return max_lagging_iterations;
	}

	bool FullNLProblem::uses_lagging() const
	{
		for (auto &f : forms_)
			if (f->uses_lagging())
				return true;
		return false;
	}

	void FullNLProblem::line_search_begin(const TVector &x0, const TVector &x1)
	{
		for (auto &f : forms_)
			f->line_search_begin(x0, x1);
	}

	void FullNLProblem::line_search_end()
	{
		for (auto &f : forms_)
			f->line_search_end();
	}

	double FullNLProblem::max_step_size(const TVector &x0, const TVector &x1)
	{
		double step = 1;
		for (auto &f : forms_)
			if (f->enabled())
				step = std::min(step, f->max_step_size(x0, x1));
		return step;
	}

	bool FullNLProblem::is_step_valid(const TVector &x0, const TVector &x1)
	{
		for (auto &f : forms_)
			if (f->enabled() && !f->is_step_valid(x0, x1))
				return false;
		return true;
	}

	bool FullNLProblem::is_step_collision_free(const TVector &x0, const TVector &x1)
	{
		for (auto &f : forms_)
			if (f->enabled() && !f->is_step_collision_free(x0, x1))
				return false;
		return true;
	}

	double FullNLProblem::value(const TVector &x)
	{
		double val = 0;
		for (auto &f : forms_)
			if (f->enabled())
				val += f->value(x);
		return val;
	}

	void FullNLProblem::gradient(const TVector &x, TVector &grad)
	{
		grad = TVector::Zero(x.size());
		for (auto &f : forms_)
		{
			if (!f->enabled())
				continue;
			TVector tmp;
			f->first_derivative(x, tmp);

			std::cout << f->name() << ": " << tmp.norm() << std::endl;

			grad += tmp;
		}
	}

	void FullNLProblem::hessian(const TVector &x, THessian &hessian)
	{
		std::vector<std::set<int>> global_element_indices;
		hessian.resize(x.size(), x.size());
		THessian contact_hessian;
		contact_hessian.resize(x.size(), x.size());
		for (auto &f : forms_)
		{
			if (!f->enabled())
				continue;

			THessian tmp;

			switch (projection_setting) {
				case 0:
					f->set_project_to_psd(project_to_psd_);
					f->second_derivative(x, tmp);
					hessian += tmp;
					break;
				case 1:
					if (f->name() == "contact")
					{
						f->set_project_to_psd(true);
						f->second_derivative(x, tmp);
						contact_hessian += tmp;
					}
					else
					{
						f->set_project_to_psd(false);
						f->second_derivative(x, tmp);

						global_element_indices.insert(global_element_indices.end(), f->global_element_indices.begin(), f->global_element_indices.end());
						
						hessian += tmp;
					}
					break;
				case 2:
					f->set_project_to_psd(false);
					f->second_derivative(x, tmp);

					global_element_indices.insert(global_element_indices.end(), f->global_element_indices.begin(), f->global_element_indices.end());
					
					hessian += tmp;
					break;
				default:
					assert(false);
			}	
		}

		if (project_to_psd_ && projection_setting > 0 && global_element_indices.size() > 0)
		{
			std::set<int> elements_to_project;
			Eigen::MatrixXd dense_hessian = hessian;
			for (int e = 0; e < global_element_indices.size(); ++e)
			{
				std::set<int> indices = global_element_indices[e];
				Eigen::MatrixXd local_hessian(indices.size(), indices.size());
				int i = 0;
				for (auto gi : indices)
				{
					int j = 0;
					for (auto gj : indices)
					{
						local_hessian(i, j) = dense_hessian(gi, gj);
						++j;
					}
					++i;
				}
				auto projected_local_hessian = ipc::project_to_psd(local_hessian);
				if (!projected_local_hessian.isApprox(local_hessian))
				{
					elements_to_project.insert(e);
				}
			}

			for (int e : elements_to_project)
			{
				std::set<int> indices = global_element_indices[e];
				Eigen::MatrixXd local_hessian(indices.size(), indices.size());
				int i = 0;
				for (auto gi : indices)
				{
					int j = 0;
					for (auto gj : indices)
					{
						local_hessian(i, j) = dense_hessian(gi, gj);
						++j;
					}
					++i;
				}
				auto projected_local_hessian = ipc::project_to_psd(local_hessian);
				i = 0;
				for (auto gi : indices)
				{
					int j = 0;
					for (auto gj : indices)
					{
						dense_hessian(gi, gj) = projected_local_hessian(i, j);
						++j;
					}
					++i;
				}
			}
			hessian = dense_hessian.sparseView();
			global_element_indices.clear();
		}

		if (projection_setting == 1)
		{
			hessian += contact_hessian;
		}
	}

	void FullNLProblem::solution_changed(const TVector &x)
	{
		for (auto &f : forms_)
			f->solution_changed(x);
	}

	void FullNLProblem::post_step(const polysolve::nonlinear::PostStepData &data)
	{
		for (auto &f : forms_)
			f->post_step(data);
	}
} // namespace polyfem::solver
