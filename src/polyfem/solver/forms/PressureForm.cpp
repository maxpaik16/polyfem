#include "PressureForm.hpp"

#include <polyfem/io/Evaluator.hpp>
#include <polyfem/utils/MaybeParallelFor.hpp>
#include <polyfem/utils/BoundarySampler.hpp>

#include <polyfem/basis/ElementBases.hpp>
#include <polyfem/assembler/AssemblerUtils.hpp>
#include <polyfem/assembler/AssemblyValsCache.hpp>

#include <polyfem/autogen/auto_p_bases.hpp>

#include <polyfem/solver/AdjointTools.hpp>

#include <polyfem/utils/Logger.hpp>

namespace polyfem::solver
{
	namespace
	{
		class LocalThreadVecStorage
		{
		public:
			Eigen::MatrixXd vec;
			assembler::ElementAssemblyValues vals, gvals;
			QuadratureVector da;

			LocalThreadVecStorage(const int size)
			{
				vec.resize(size, 1);
				vec.setZero();
			}
		};
	} // namespace
	PressureForm::PressureForm(const int ndof,
							   const std::vector<mesh::LocalBoundary> &local_pressure_boundary,
							   const std::vector<int> dirichlet_nodes,
							   const int n_boundary_samples,
							   const assembler::PressureAssembler &pressure_assembler,
							   const bool is_time_dependent)
		: ndof_(ndof),
		  local_pressure_boundary_(local_pressure_boundary),
		  dirichlet_nodes_(dirichlet_nodes),
		  n_boundary_samples_(n_boundary_samples),
		  pressure_assembler_(pressure_assembler)
	{
		t_ = 0;
	}

	double PressureForm::value_unweighted(const Eigen::VectorXd &x) const
	{
		return -1 * pressure_assembler_.compute_energy(x, local_pressure_boundary_, n_boundary_samples_, t_);
	}

	void PressureForm::first_derivative_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const
	{
		// REMEMBER -!!!!!
		pressure_assembler_.compute_energy_grad(x, local_pressure_boundary_, dirichlet_nodes_, n_boundary_samples_, t_, gradv);
		gradv *= -1;
	}

	void PressureForm::second_derivative_unweighted(const Eigen::VectorXd &x, StiffnessMatrix &hessian) const
	{
		pressure_assembler_.compute_energy_hess(x, local_pressure_boundary_, dirichlet_nodes_, n_boundary_samples_, t_, project_to_psd_, hessian);
		hessian *= -1;
	}

	void PressureForm::update_quantities(const double t, const Eigen::VectorXd &x)
	{
		this->t_ = t;
	}

} // namespace polyfem::solver
