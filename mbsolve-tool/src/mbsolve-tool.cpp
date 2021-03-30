/*
 * mbsolve: An open-source solver tool for the Maxwell-Bloch equations.
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

/**
 * \defgroup MBSOLVE_TOOL mbsolve-tool
 * Runs different simulation setups.
 */

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <regex>
#include <string>
#include <vector>
#include <cxxopts.hpp>
#include <mbsolve/lib.hpp>

#if USE_CPU
#include <mbsolve/solver-cpu/solver_cpu_loader.hpp>
#endif

#if USE_HDF5
#include <mbsolve/writer-hdf5/writer_hdf5.hpp>
#endif

static std::string device_file;
static std::string output_file;
static std::string scenario_file;
static std::string solver_method;
static std::string writer_method;
static mbsolve::real sim_endtime;
static unsigned int num_gridpoints;

static void
parse_args(int argc, char** argv)
{
    /* set up program name and help one-liner */
    cxxopts::Options options(
        argv[0],
        "mbsolve: An open-source solver "
        "tool for the Maxwell-Bloch equations.");

    /* set up allowed options */
    options.add_options()("h,help", "Print usage")(
        "d,device",
        "Specify device",
        cxxopts::value<std::string>(device_file))(
        "m,method",
        "Specify solver method",
        cxxopts::value<std::string>(solver_method))(
        "o,output",
        "Specify output file",
        cxxopts::value<std::string>(output_file))(
        "s,scenario",
        "Specify scenario",
        cxxopts::value<std::string>(scenario_file))(
        "w,writer",
        "Specify writer",
        cxxopts::value<std::string>(writer_method))(
        "e,endtime",
        "Specify simulation end time",
        cxxopts::value<mbsolve::real>(sim_endtime))(
        "g,gridpoints",
        "Specify number of spatial grid points",
        cxxopts::value<unsigned int>(num_gridpoints));

    try {
        auto result = options.parse(argc, argv);

        /* help request */
        if (result.count("help")) {
            std::cout << options.help();
            exit(0);
        }

        /* mandatory arguments */
        if (!result.count("device")) {
            throw std::invalid_argument("No device specified");
        } else {
            std::cout << "Using device " << device_file << std::endl;
        }

        if (!result.count("method")) {
            throw std::invalid_argument("No method specified");
        }

        if (!result.count("writer")) {
            throw std::invalid_argument("No writer specified");
        }
    } catch (const std::exception& ex) {
        std::cout << options.help();
        std::cout << ex.what() << "." << std::endl;
        exit(1);
    }
}

class wallclock
{
private:
    typedef std::chrono::duration<double> duration;

    std::chrono::time_point<std::chrono::steady_clock> start;

public:
    wallclock() { tic(); }

    void tic() { start = std::chrono::steady_clock::now(); }

    double toc()
    {
        duration elapsed = std::chrono::steady_clock::now() - start;
        return elapsed.count();
    }
};

/**
 * mbsolve-tool main function.
 * \ingroup MBSOLVE_TOOL
 *
 * Specify the simulation setup using the -d parameter (the available setups
 * are described below) and the solver method using the -m parameter.
 * For the complete list of parameters, run mbsolve-tool -h.
 */
int
main(int argc, char** argv)
{
    /* parse command line arguments */
    parse_args(argc, argv);

    try {
#if USE_CPU
        mbsolve::solver_cpu_loader cpu_load;
#endif

#if USE_HDF5
        mbsolve::writer_hdf5_loader hdf5_load;
#endif

        wallclock clock;
        double total_time = 0;

        std::shared_ptr<mbsolve::device> dev;
        std::shared_ptr<mbsolve::scenario> scen;

        if (device_file == "song2005") {
            /**
             * The song2005 setup features a three-level active region which
             * is excited by a sech pulse. For details see literature:
             * Song X. et al., Propagation of a Few-Cycle Laser Pulse
             * in a V-Type Three-Level System, Optics and Spectroscopy,
             * Vol. 99, No. 4, 2005, pp. 517–521
             * https://doi.org/10.1134/1.2113361
             */

            /* set up quantum mechanical description */
            std::vector<mbsolve::real> energies = {
                0,
                2.3717 * mbsolve::HBAR * 1e15,
                2.4165 * mbsolve::HBAR * 1e15
            };

            mbsolve::qm_operator H(energies);

            std::vector<std::complex<mbsolve::real> > dipoles = {
                mbsolve::E0 * 9.2374e-11,
                mbsolve::E0 * 9.2374e-11 * sqrt(2),
                0
            };

            mbsolve::qm_operator u({ 0, 0, 0 }, dipoles);

            mbsolve::real rate = 1e10;
            std::vector<std::vector<mbsolve::real> > scattering_rates = {
                { 0, rate, rate }, { rate, 0, rate }, { rate, rate, 0 }
            };

            auto relax_sop =
                std::make_shared<mbsolve::qm_lindblad_relaxation>(
                    scattering_rates);

            auto qm = std::make_shared<mbsolve::qm_description>(
                6e24, 0, H, u, relax_sop);

            auto mat_vac = std::make_shared<mbsolve::material>("Vacuum");
            mbsolve::material::add_to_library(mat_vac);
            auto mat_ar = std::make_shared<mbsolve::material>("AR_Song", qm);
            mbsolve::material::add_to_library(mat_ar);

            dev = std::make_shared<mbsolve::device>("Song");
            dev->add_region(std::make_shared<mbsolve::region>(
                "Active region (single point)", mat_ar, 0, 0));

            /* default settings */
            if (num_gridpoints == 0) {
                num_gridpoints = 1;
            }
            if (sim_endtime < 1e-21) {
                sim_endtime = 80e-15;
            }

            mbsolve::qm_operator rho_init({ 1, 0, 0 });

            /* Song basic scenario */
            auto ic_d = std::make_shared<mbsolve::ic_density_const>(rho_init);
            auto ic_e = std::make_shared<mbsolve::ic_field_const>(0.0);
            auto ic_m = std::make_shared<mbsolve::ic_field_const>(0.0);

            scen = std::make_shared<mbsolve::scenario>(
                "Basic", 1, sim_endtime, ic_d, ic_e, ic_m, 10000);

            auto sech_pulse = std::make_shared<mbsolve::sech_pulse>(
                "sech",
                0.0,
                mbsolve::source::hard_source,
                3.5471e9,
                3.8118e14,
                17.248,
                1.76 / 5e-15,
                -mbsolve::PI / 2);
            scen->add_source(sech_pulse);

            scen->add_record(std::make_shared<mbsolve::record>(
                "d11", mbsolve::record::type::density, 1, 1, 0, 0));
            scen->add_record(std::make_shared<mbsolve::record>(
                "d22", mbsolve::record::type::density, 2, 2, 0, 0));
            scen->add_record(std::make_shared<mbsolve::record>(
                "d33", mbsolve::record::type::density, 3, 3, 0, 0));
            scen->add_record(std::make_shared<mbsolve::record>("e", 0, 0.0));

        } else if (device_file == "ziolkowski1995") {
            /**
             * The ziolkowski1995 setup is a self induced transparency (SIT)
             * setup that consists of a two-level active region embedded in
             * two vacuum section. For details see literature:
             * https://doi.org/10.1103/PhysRevA.52.3082
             */

            /* set up quantum mechanical description */
            auto qm = std::make_shared<mbsolve::qm_desc_2lvl>(
                1e24, 0, 2 * mbsolve::PI * 2e14, 6.24e-11, 1.0e10, 1.0e10);

            /* materials */
            auto mat_vac = std::make_shared<mbsolve::material>("Vacuum");
            auto mat_ar =
                std::make_shared<mbsolve::material>("AR_Ziolkowski", qm);
            mbsolve::material::add_to_library(mat_vac);
            mbsolve::material::add_to_library(mat_ar);

            /* set up device */
            dev = std::make_shared<mbsolve::device>("Ziolkowski");
            dev->add_region(std::make_shared<mbsolve::region>(
                "Vacuum left", mat_vac, 0, 7.5e-6));
            dev->add_region(std::make_shared<mbsolve::region>(
                "Active region", mat_ar, 7.5e-6, 142.5e-6));
            dev->add_region(std::make_shared<mbsolve::region>(
                "Vacuum right", mat_vac, 142.5e-6, 150e-6));

            /* default settings */
            if (num_gridpoints == 0) {
                num_gridpoints = 32768;
            }
            if (sim_endtime < 1e-21) {
                sim_endtime = 200e-15;
            }

            mbsolve::qm_operator rho_init({ 1, 0, 0 });

            /* Ziolkowski basic scenario */
            auto ic_d = std::make_shared<mbsolve::ic_density_const>(rho_init);
            auto ic_e = std::make_shared<mbsolve::ic_field_const>(0.0);

            scen = std::make_shared<mbsolve::scenario>(
                "Basic", num_gridpoints, sim_endtime, ic_d, ic_e);

            auto sech_pulse = std::make_shared<mbsolve::sech_pulse>
                //("sech", 0.0, mbsolve::source::hard_source, 4.2186e9/2,
                // 2e14,
                ("sech",
                 0.0,
                 mbsolve::source::hard_source,
                 4.2186e9,
                 2e14,
                 10,
                 2e14);
            scen->add_source(sech_pulse);

            scen->add_record(
                std::make_shared<mbsolve::record>("inv12", 2.5e-15));
            scen->add_record(std::make_shared<mbsolve::record>("e", 2.5e-15));

        } else if (device_file == "andreasen2009") {
            /**
             * The andreasen2009 setup is an atomic superfluorescence
             * setup that consists of a two-level active region embedded in
             * two vacuum section. For details see literature:
             * https://doi.org/10.1103/PhysRevA.77.023810
             */

            mbsolve::real num_carrier_cell = 3e4;

            /* set up quantum mechanical description */
            auto qm = std::make_shared<mbsolve::qm_desc_2lvl>(
                2.13e19,
                num_carrier_cell,
                2 * mbsolve::PI * 4.77e14,
                6.875e-11,
                1.32e7,
                1.0e10);

            /* materials */
            auto mat_vac = std::make_shared<mbsolve::material>("Vacuum");
            auto mat_ar =
                std::make_shared<mbsolve::material>("AR_Andreasen", qm);
            mbsolve::material::add_to_library(mat_vac);
            mbsolve::material::add_to_library(mat_ar);

            /* set up device */
            dev = std::make_shared<mbsolve::device>("Andreasen");
            dev->add_region(std::make_shared<mbsolve::region>(
                "Vacuum left", mat_vac, 0, 70e-6));
            dev->add_region(std::make_shared<mbsolve::region>(
                "Active region", mat_ar, 70e-6, 7.07e-3));
            dev->add_region(std::make_shared<mbsolve::region>(
                "Vacuum right", mat_vac, 7.07e-3, 7.14e-3));

            /* default settings */
            if (num_gridpoints == 0) {
                num_gridpoints = 102000;
            }
            if (sim_endtime < 1e-21) {
                sim_endtime = 1e-9;
            }

            /* Andreasen basic scenario */
            auto ic_d = std::make_shared<mbsolve::ic_density_random_2lvl>(
                num_carrier_cell);
            auto ic_e = std::make_shared<mbsolve::ic_field_const>(0.0);

            scen = std::make_shared<mbsolve::scenario>(
                "Basic", num_gridpoints, sim_endtime, ic_d, ic_e);
            scen->set_boundary_flag(1);
            
            scen->add_record(
                std::make_shared<mbsolve::record>("inv12", 1e-13));
            scen->add_record(std::make_shared<mbsolve::record>("e", 1e-13));

        } else if (device_file == "tzenov2016") {
            /**
             * The tzenov2016 setup is a model of an actual quantum cascade
             * laser frequency comb.
             * For details see literature:
             * https://doi.org/10.1364/OE.24.023232
             */

            /* quantum mechanical description of active region */
            mbsolve::qm_operator H(
                { 0.10103 * mbsolve::E0,
                  0.09677 * mbsolve::E0,
                  0.09720 * mbsolve::E0,
                  0.08129 * mbsolve::E0,
                  0.07633 * mbsolve::E0 },
                { 0.0,
                  1.2329e-3 * mbsolve::E0,
                  -1.3447e-3 * mbsolve::E0,
                  0.0,
                  0.0,
                  0.0,
                  0.0,
                  0.0,
                  0.0,
                  0.0 });

            mbsolve::qm_operator u(
                { 0.0, 0.0, 0.0, 0.0, 0.0 },
                { 0.0,
                  0.0,
                  0.0,
                  0.0,
                  0.0,
                  mbsolve::E0 * 4e-9,
                  0.0,
                  0.0,
                  0.0,
                  0.0 });

            mbsolve::real deph_r1 = 1.0 / (1.6e-12);
            mbsolve::real deph_r2 = 1.0 / (4.8e-12);

            std::vector<std::vector<mbsolve::real> > scattering_rates = {
                { deph_r1, 0.4947e12, 0.0974e12, 0.8116e12, 1.0410e12 },
                { 0.8245e12, deph_r1, 0.1358e12, 0.6621e12, 1.1240e12 },
                { 0.0229e12, 0.0469e12, deph_r1, 0.0794e12, 0.0357e12 },
                { 0.0047e12, 0.0029e12, 0.1252e12, deph_r2, 0.2810e12 },
                { 0.0049e12, 0.0049e12, 0.1101e12, 0.4949e12, deph_r2 }
            };

            auto relax_sop =
                std::make_shared<mbsolve::qm_lindblad_relaxation>(
                    scattering_rates);

            auto qm = std::make_shared<mbsolve::qm_description>(
                5.6e21, 6.84e5, H, u, relax_sop);

            auto mat_ar =
                std::make_shared<mbsolve::material>("AR", qm, 17, 0.9, 1400);
            mbsolve::material::add_to_library(mat_ar);

            /* set up device with semi-transparent mirror boundary cond. */
            auto bc =
                std::make_shared<mbsolve::bc_field_reflectivity>(0.8, 0.8);
            dev = std::make_shared<mbsolve::device>("tzenov2016", bc);
            dev->add_region(std::make_shared<mbsolve::region>(
                "Active region", mat_ar, 0, 5e-3));

            /* default settings */
            if (num_gridpoints == 0) {
                num_gridpoints = 8192;
            }
            if (sim_endtime < 1e-21) {
                sim_endtime = 0.2e-9;
            }

            /* initial density matrix */
            mbsolve::qm_operator rho_init({ 0.0, 0.0, 1.0, 0.0, 0.0 });

            /* basic scenario */
            scen = std::make_shared<mbsolve::scenario>(
                "basic", num_gridpoints, sim_endtime, rho_init);
            scen->add_record(std::make_shared<mbsolve::record>(
                "e0", mbsolve::record::electric, 1, 1, 0.0, 0.0));
            for (int i = 1; i <= 5; i++) {
                scen->add_record(std::make_shared<mbsolve::record>(
                    "d" + std::to_string(i) + std::to_string(i),
                    mbsolve::record::type::density,
                    i,
                    i,
                    0.0,
                    0.1e-3));
            }

        } else if (device_file == "tzenov2018-cpml") {
            /**
             * The tzenov2018-cpml setup consists of an absorber region
             * embedded in two gain regions. Each region is modeled as a
             * two-level system. The results show that short pulses are
             * generated due to colliding pulse mode-locking (CPML).
             * For details see literature:
             * https://doi.org/10.1088/1367-2630/aac12a
             */

            /* set up quantum mechanical descriptions */
            auto qm_gain = std::make_shared<mbsolve::qm_desc_2lvl>(
                5e21,
                0,
                2 * mbsolve::PI * 3.4e12,
                2e-9,
                1.0 / 10e-12,
                1.0 / 200e-15,
                1.0);

            auto qm_absorber = std::make_shared<mbsolve::qm_desc_2lvl>(
                1e21,
                0,
                2 * mbsolve::PI * 3.4e12,
                6e-9,
                1.0 / 3e-12,
                1.0 / 160e-15);

            /* materials */
            auto mat_absorber = std::make_shared<mbsolve::material>(
                "Absorber", qm_absorber, 12.96, 1.0, 500);
            auto mat_gain = std::make_shared<mbsolve::material>(
                "Gain", qm_gain, 12.96, 1.0, 500);
            mbsolve::material::add_to_library(mat_absorber);
            mbsolve::material::add_to_library(mat_gain);

            /* set up device */
            dev = std::make_shared<mbsolve::device>("tzenov-cpml");
            dev->add_region(std::make_shared<mbsolve::region>(
                "Gain R", mat_gain, 0, 0.5e-3));
            dev->add_region(std::make_shared<mbsolve::region>(
                "Absorber", mat_absorber, 0.5e-3, 0.625e-3));
            dev->add_region(std::make_shared<mbsolve::region>(
                "Gain L", mat_gain, 0.625e-3, 1.125e-3));

            /* default settings */
            if (num_gridpoints == 0) {
                num_gridpoints = 8192;
            }
            if (sim_endtime < 1e-21) {
                sim_endtime = 2e-9;
            }

            mbsolve::qm_operator rho_init({ 0.5, 0.5 }, { 0.001 });

            /* basic scenario */
            scen = std::make_shared<mbsolve::scenario>(
                "basic", num_gridpoints, sim_endtime, rho_init);

            scen->add_record(
                std::make_shared<mbsolve::record>("inv12", 1e-12));
            scen->add_record(std::make_shared<mbsolve::record>("e", 1e-12));
            scen->add_record(std::make_shared<mbsolve::record>(
                "e0", mbsolve::record::electric, 1, 1, 0.0, 0.0));
            scen->add_record(std::make_shared<mbsolve::record>(
                "h0", mbsolve::record::magnetic, 1, 1, 0.0, 1.373e-7));

        } else if (device_file.rfind("marskar2011", 0) == 0) {
            /**
             * The marskar2011 setup is a 6-level anharmonic ladder
             * system. In the scenario, the interaction with a few-cycle
             * Gaussian pulse is simulated. For details see literature:
             * https://doi.org/10.1364/OE.19.016784
             *
             * This setup has been generalized to an arbitrary number of
             * levels N >= 2, which can be specified by setting the
             * device_file to "marskar2011-Nlvl". If the expression cannot
             * be parsed or the original device_file = "marskar2011" is
             * used, N defaults to 6.
             *
             */

            /* number of levels, defaults to 6 */
            int N = 6;

            /* parse number of levels*/
            std::string num_lvl = device_file;
            num_lvl = std::regex_replace(
                num_lvl, std::regex("^marskar2011[-]?"), "");
            num_lvl = std::regex_replace(num_lvl, std::regex("lvl$"), "");

            try {
                N = std::stoi(num_lvl);
                if (N < 2) {
                    throw std::invalid_argument("N must be >= 2!");
                }
            } catch (const std::exception& e) {
                std::cout << "Warning: Could not determine number of levels."
                          << " Assuming N = 6." << std::endl;
                N = 6;
            }

            /* set up quantum mechanical description */
            std::vector<mbsolve::real> energies(N, 0.0);
            for (int i = 1; i < N; i++) {
                energies[i] = energies[i - 1] +
                    (1.0 - 0.1 * (i - 3)) * mbsolve::HBAR * 2 * mbsolve::PI *
                        1e13;
            }

            mbsolve::qm_operator H(energies);

            mbsolve::real dipole = 1e-29;
            std::vector<std::complex<mbsolve::real> > dipoles;
            for (int i = 1; i < N; i++) {
                for (int j = 0; j < i - 1; j++) {
                    dipoles.push_back(0.0);
                }
                dipoles.push_back(dipole);
            }

            mbsolve::qm_operator u(
                std::vector<mbsolve::real>(N, 0.0), dipoles);

            /* scattering rates are temperature dependent (use T = 600 K)
             */
            mbsolve::real T = 600;

            /* pure dephasing rate */
            mbsolve::real pure_deph = 1e12;

            std::vector<std::vector<mbsolve::real> > scattering_rates;
            for (int i = 0; i < N; i++) {
                std::vector<mbsolve::real> v;
                for (int j = 0; j < N; j++) {
                    mbsolve::real rate;
                    if (i == j) {
                        rate = pure_deph;
                    } else if (i == j - 1) {
                        rate = 1e12 / (1 + 0.1 * i);
                    } else if (i == j + 1) {
                        rate = 1e12 / (1 + 0.1 * j) *
                            exp((energies[j] - energies[i]) /
                                (mbsolve::KB * T));
                    } else {
                        rate = 0.0;
                    }
                    v.push_back(rate);
                }
                scattering_rates.push_back(v);
            }

            auto relax_sop =
                std::make_shared<mbsolve::qm_lindblad_relaxation>(
                    scattering_rates);

            auto qm = std::make_shared<mbsolve::qm_description>(
                1e25, 0, H, u, relax_sop);

            /* materials */
            auto mat_vac = std::make_shared<mbsolve::material>("Vacuum");
            mbsolve::material::add_to_library(mat_vac);
            auto mat_ar = std::make_shared<mbsolve::material>("AR_Marsk", qm);
            mbsolve::material::add_to_library(mat_ar);

            /* set up device */
            dev = std::make_shared<mbsolve::device>(
                "marskar2011-" + std::to_string(N) + "lvl");
            dev->add_region(std::make_shared<mbsolve::region>(
                "Active region", mat_ar, 0, 1e-3));

            /* default settings */
            if (num_gridpoints == 0) {
                num_gridpoints = 8192;
            }
            if (sim_endtime < 1e-21) {
                sim_endtime = 2e-12;
            }

            /* determine equilibrium populations:
             *  - the scattering rates form a tridiagonal N x N matrix
             *  - we solve for the homogeneous solution = equilibrium
             *  - solving yields relations between neighboring populations
             *    of the form rho_22 = (gamma_21/gamma_12) rho_11
             */
            mbsolve::real total = 0.0;
            for (int i = 0; i < N; i++) {
                /* the ratios are chained, e.g.,
                 * rho_33 = (gamma_32/gamma_23) (gamma_21/gamma_12) rho_11
                 */
                mbsolve::real prod = 1;
                for (int j = 0; j < i; j++) {
                    prod *= scattering_rates[j + 1][j] /
                        scattering_rates[j][j + 1];
                }
                total += prod;
            }

            std::vector<mbsolve::real> populations_eq;
            /* rho_11 is normalized w.r.t. to the sum of products */
            populations_eq.push_back(1 / total);
            /* the remaining populations are calculated */
            for (int i = 1; i < N; i++) {
                populations_eq.push_back(
                    populations_eq[i - 1] * scattering_rates[i][i - 1] /
                    scattering_rates[i - 1][i]);
            }

            mbsolve::qm_operator rho_init(populations_eq);

            /* Marskar basic scenario */
            auto ic_d = std::make_shared<mbsolve::ic_density_const>(rho_init);
            auto ic_e = std::make_shared<mbsolve::ic_field_const>(0.0);

            scen = std::make_shared<mbsolve::scenario>(
                "basic", num_gridpoints, sim_endtime, ic_d, ic_e);

            /* input pulse */
            mbsolve::real tau = 100e-15;
            auto pulse = std::make_shared<mbsolve::gaussian_pulse>(
                "gaussian",
                0.0, /* position */
                mbsolve::source::hard_source,
                5e8, /* amplitude */
                1e13, /* frequency */
                3.0 * tau, /* phase: 3*tau */
                tau /* tau */
            );
            scen->add_source(pulse);

            /* select data to be recorded */
            mbsolve::real sample_time = 0.0;
            mbsolve::real sample_pos = 0.0;
            for (int i = 1; i <= N; i++) {
                scen->add_record(std::make_shared<mbsolve::record>(
                    "d" + std::to_string(i) + std::to_string(i),
                    mbsolve::record::type::density,
                    i,
                    i,
                    sample_time,
                    sample_pos));
            }
            scen->add_record(std::make_shared<mbsolve::record>(
                "e",
                mbsolve::record::type::electric,
                0,
                0,
                sample_time,
                sample_pos));
        } else {
            throw std::invalid_argument("Specified device not found!");
        }

        /* tic */
        clock.tic();

        std::shared_ptr<mbsolve::writer> writer =
            mbsolve::writer::create_instance(writer_method);

        std::shared_ptr<mbsolve::solver> solver =
            mbsolve::solver::create_instance(solver_method, dev, scen);

        /* toc */
        double setup_time = clock.toc();
        std::cout << "Time required (setup): " << setup_time << std::endl;
        total_time += setup_time;

        std::cout << solver->get_name() << std::endl;

        /* tic */
        clock.tic();

        /* execute solver */
        solver->run();

        /* toc */
        double run_time = clock.toc();
        std::cout << "Time required (run): " << run_time << std::endl;
        total_time += run_time;

        /* grid point updates per second */
        double gpups = 1e-6 / run_time * scen->get_num_gridpoints() *
            scen->get_endtime() / scen->get_timestep_size();
        std::cout << "Performance: " << gpups << " MGPU/s" << std::endl;

        /* tic */
        clock.tic();

        /* write results */
        std::string default_file = dev->get_name() + "_" + scen->get_name() +
            "." + writer->get_extension();
        writer->write(
            output_file.empty() ? default_file : output_file,
            solver->get_results(),
            dev,
            scen);

        /* toc */
        double write_time = clock.toc();
        std::cout << "Time required (write): " << write_time << std::endl;
        total_time += write_time;

        std::cout << "Time required (total): " << total_time << std::endl;
    } catch (const std::exception& ex) {
        std::cout << "Error: " << ex.what() << std::endl;
        exit(1);
    }

    exit(0);
}
