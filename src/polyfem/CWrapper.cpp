
#pragma once

#include <polyfem/State.hpp>
#include <polyfem/solver/NLProblem.hpp>
#include <polyfem/time_integrator/ImplicitTimeIntegrator.hpp>
#include <polyfem/utils/Types.hpp>

#include <Eigen/Sparse>

extern "C"
{
    using namespace polyfem;

    void setup_sim(void* void_state_ptr, int argc, char **argv)
    {

    }

    double get_objective(void* void_state_ptr, double* x, unsigned int ndofs)
    {
        State* state_ptr = static_cast<State*>(void_state_ptr);
        if (!state_ptr)
        {
            assert(false);
        }
        Eigen::VectorXd eigen_x(ndofs);
        for (int i = 0; i < ndofs; ++i)
        {
            eigen_x(i) = x[i];
        }
        return state_ptr->solve_data.nl_problem->value(eigen_x);
    }

    void get_gradient(void* void_state_ptr, double* x, double* grad, unsigned int ndofs)
    {
        State* state_ptr = static_cast<State*>(void_state_ptr);
        if (!state_ptr)
        {
            assert(false);
        }
        Eigen::VectorXd eigen_x(ndofs);
        for (int i = 0; i < ndofs; ++i)
        {
            eigen_x(i) = x[i];
        }
        Eigen::VectorXd eigen_grad(ndofs);
        state_ptr->solve_data.nl_problem->gradient(eigen_x, eigen_grad);

        for (int i = 0; i < ndofs; ++i)
        {
            grad[i] = eigen_grad(i);
        }
    }

    void get_hessian(void* void_state_ptr, double* x, unsigned int ndofs, unsigned int* rows, unsigned int* cols, double* vals, unsigned int* nvals)
    {
        State* state_ptr = static_cast<State*>(void_state_ptr);
        if (!state_ptr)
        {
            assert(false);
        }
        Eigen::VectorXd eigen_x(ndofs);
        for (int i = 0; i < ndofs; ++i)
        {
            eigen_x(i) = x[i];
        }
        StiffnessMatrix hessian;
        state_ptr->solve_data.nl_problem->hessian(eigen_x, hessian);
        
        *nvals = hessian.nonZeros();

        rows = static_cast<unsigned int*>(malloc(sizeof(unsigned int) * (*nvals)));
        cols = static_cast<unsigned int*>(malloc(sizeof(unsigned int) * (*nvals)));
        vals = static_cast<double*>(malloc(sizeof(double) * (*nvals)));
        if (!(rows && cols && vals))
        {
            assert(fales);
        }

        unsigned int i = 0;
        for (int k = 0; k < hessian.outerSize(); ++k)
        {
            for (StiffnessMatrix::InnerIterator it(hessian, k); it; ++it)
            {
                rows[i] = it.row();
                cols[i] = it.col();
                vals[i] = it.value();
                ++i;
            }
        }
    }

    void get_constraints(void* void_state_ptr, double* x, unsigned int ndofs, double* constraints, unsigned int* nconstraints)
    {

    }

    void make_timestep(void* void_state_ptr, double* x, unsigned int ndofs, int t, double dt)
    {
        State* state_ptr = static_cast<State*>(void_state_ptr);
        if (!state_ptr)
        {
            assert(false);
        }
        Eigen::VectorXd eigen_x(ndofs);
        for (int i = 0; i < ndofs; ++i)
        {
            eigen_x(i) = x[i];
        }

        state_ptr->save_timestep(dt * t, t, 0, dt, eigen_x, Eigen::MatrixXd()); // no pressure

        {
            state_ptr->solve_data.time_integrator->update_quantities(eigen_x);

            state_ptr->solve_data.nl_problem->update_quantities((t + 1) * dt, eigen_x);

            state_ptr->solve_data.update_dt();
            state_ptr->solve_data.update_barrier_stiffness(eigen_x);
        }
    }
    
}