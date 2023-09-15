#pragma once

#include <polyfem/assembler/Assembler.hpp>
#include <polyfem/mesh/Obstacle.hpp>

#include <polyfem/assembler/Problem.hpp>
#include <polyfem/assembler/MatParams.hpp>
#include <polyfem/mesh/LocalBoundary.hpp>

namespace polyfem
{
	namespace assembler
	{
		// computes the rhs of a problem by \int \phi rho rhs
		class PressureAssembler
		{
		public:
			// initialization with assembler factory mesh
			// size of the problem, bases
			// and solver used internally
			PressureAssembler(const Assembler &assembler, const mesh::Mesh &mesh, const mesh::Obstacle &obstacle,
							  const std::vector<mesh::LocalBoundary> &local_pressure_boundary,
							  const int n_basis, const int size,
							  const std::vector<basis::ElementBases> &bases, const std::vector<basis::ElementBases> &gbases, const Problem &problem);

			double compute_energy(
				const Eigen::MatrixXd &displacement,
				const std::vector<mesh::LocalBoundary> &local_pressure_boundary,
				const int resolution,
				const double t) const;
			void compute_energy_grad(
				const Eigen::MatrixXd &displacement,
				const std::vector<mesh::LocalBoundary> &local_pressure_boundary,
				const std::vector<int> dirichlet_nodes,
				const int resolution,
				const double t,
				Eigen::VectorXd &grad) const;
			void compute_energy_hess(
				const Eigen::MatrixXd &displacement,
				const std::vector<mesh::LocalBoundary> &local_pressure_boundary,
				const std::vector<int> dirichlet_nodes,
				const int resolution,
				const double t,
				const bool project_to_psd,
				StiffnessMatrix &hess) const;

			inline const Problem &problem() const { return problem_; }
			inline const mesh::Mesh &mesh() const { return mesh_; }
			inline const std::vector<basis::ElementBases> &bases() const { return bases_; }
			inline const std::vector<basis::ElementBases> &gbases() const { return gbases_; }
			inline const Assembler &assembler() const { return assembler_; }

		private:
			const Assembler &assembler_;
			const mesh::Mesh &mesh_;
			const mesh::Obstacle &obstacle_;
			const int n_basis_;
			const int size_;
			const std::vector<basis::ElementBases> &bases_;
			const std::vector<basis::ElementBases> &gbases_;
			const Problem &problem_;
		};
	} // namespace assembler
} // namespace polyfem
