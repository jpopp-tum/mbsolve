/*
 * mbsolve: Framework for solving the Maxwell-Bloch/-Lioville equations
 *
 * Copyright (c) 2016, Computational Photonics Group, Technical University of
 * Munich.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 */

#ifndef MBSOLVE_SOLVER_OPENMP_2LVL_PC_H
#define MBSOLVE_SOLVER_OPENMP_2LVL_PC_H

#include <copy_list_entry_int.hpp>
#include <solver.hpp>

namespace mbsolve {

/**
 * \defgroup MBSOLVE_SOLVER_OPENMP solver-openmp
 * Different solvers that use OpenMP for parallelization.
 */

/**
 *
 * \ingroup MBSOLVE_SOLVER_OPENMP
 */
class sim_constants_2lvl
{
public:
    /* electromagnetic constants */
    real M_CE;
    real M_CH;
    real M_CP;
    real sigma;

    /* quantum mechanical constants */
    real w12;
    real d12;
    real tau1;
    real gamma12;

    /* simulation constants */
    real d_x_inv;
    real d_t;

    /* initialization constants */
    real dm11_init;
    real dm22_init;
};

class sim_source
{
public:
    source::type type;
    unsigned int x_idx;
    unsigned int data_base_idx;
};

/**
 * OpenMP solver for 2-lvl systems using the predictor corrector approach.
 * \ingroup MBSOLVE_SOLVER_OPENMP
 */
class solver_openmp_2lvl_pc : public solver_int
{
public:
    solver_openmp_2lvl_pc(std::shared_ptr<const device> dev,
                          std::shared_ptr<scenario> scen);

    ~solver_openmp_2lvl_pc();

    const std::string& get_name() const;

    void run() const;

private:

    /* TODO: rule of three. make copy constructor etc. private?
     * or implement correctly
     */

    real *m_dm11;
    real *m_dm12r;
    real *m_dm12i;
    real *m_dm22;

    real *m_h;
    real *m_e;

    real *m_result_scratch;

    real *m_source_data;

    unsigned int *m_mat_indices;

    std::vector<sim_constants_2lvl> m_sim_consts;

    std::vector<sim_source> m_sim_sources;

    std::vector<copy_list_entry> m_copy_list;
};

}

#endif
