/**
 * DEM Particle Simulator – Main Driver
 * HPSC / COMP528 Group Assignment 1
 *
 * ALL output is plain CSV in results/.  No matplotlib needed on cluster.
 * After a run, 'make zip' packages results/ for local download.
 * Plots are generated locally: cd scripts/local_plots && python3 plot_all.py
 *
 * Module pipeline (each is a Fortran 90 subroutine):
 *   init_particles → zero_forces → add_gravity →
 *   compute_particle_contacts → compute_wall_contacts →
 *   integrate_particles → compute_kinetic_energy → (CSV output)
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <cstdlib>
#include "dem_params.h"
#include "dem_io.h"
#include "timer.h"

// ── Fortran 90 interface ─────────────────────────────────────────────────────
extern "C" {
    void init_particles_(
        double* x, double* y, double* z,
        double* vx, double* vy, double* vz,
        double* mass, double* radius,
        int* n, int* mode,
        double* lx, double* ly, double* lz,
        unsigned int* seed);

    void zero_forces_(double* fx, double* fy, double* fz, int* n);

    void add_gravity_(
        double* fx, double* fy, double* fz,
        double* mass, int* n,
        double* gx, double* gy, double* gz);

    void compute_particle_contacts_(
        const double* x, const double* y, const double* z,
        const double* vx,const double* vy,const double* vz,
        double* fx, double* fy, double* fz,
        const double* radius, int* n,
        double* kn, double* gamma_n,
        int* ncontacts);

    void compute_wall_contacts_(
        const double* x, const double* y, const double* z,
        const double* vx,const double* vy,const double* vz,
        double* fx, double* fy, double* fz,
        const double* radius, int* n,
        double* kn, double* gamma_n,
        double* lx, double* ly, double* lz);

    void integrate_particles_(
        double* x, double* y, double* z,
        double* vx,double* vy,double* vz,
        const double* fx,const double* fy,const double* fz,
        const double* mass, int* n, double* dt);

    void compute_kinetic_energy_(
        const double* vx,const double* vy,const double* vz,
        const double* mass, int* n, double* ke);

    void set_omp_threads_(int* nthreads);
    int  get_omp_threads_();
}

// ── Utility ──────────────────────────────────────────────────────────────────
static void ensure_dir(const std::string& d) {
    std::string cmd = "mkdir -p " + d;
    std::system(cmd.c_str());
}

static void print_banner(const SimParams& p, int threads) {
    std::cout << "\n═════════════════════════════════════════\n";
    std::cout <<   "  DEM Simulator  COMP528 2026 Assign 1   \n";
    std::cout <<   "═════════════════════════════════════════\n";
    std::cout << "  Particles    : " << p.n_particles << "\n";
    std::cout << "  OMP Threads  : " << threads       << "\n";
    std::cout << "  dt           : " << p.dt          << " s\n";
    std::cout << "  T_final      : " << p.t_final     << " s\n";
    std::cout << "  Init mode    : " << p.init_mode   << "\n";
    std::cout << "  kn / gamma_n : " << p.kn << " / " << p.gamma_n << "\n";
    std::cout << "  Domain       : [0," << p.lx << "]^3\n\n";
}

// ════════════════════════════════════════════════════════════════════════════
int main(int argc, char* argv[]) {

    // ── 1. Load parameters ───────────────────────────────────────────────────
    SimParams p = default_params();
    if (argc > 1) load_params(argv[1], p);

    // ── 2. OpenMP ────────────────────────────────────────────────────────────
    set_omp_threads_(&p.n_threads);
    int actual_threads = get_omp_threads_();
    print_banner(p, actual_threads);

    int N = p.n_particles;

    // ── 3. Particle arrays ───────────────────────────────────────────────────
    std::vector<double> x(N),  y(N),  z(N);
    std::vector<double> vx(N), vy(N), vz(N);
    std::vector<double> fx(N), fy(N), fz(N);
    std::vector<double> mass(N), radius(N);

    // ── 4. Initialise particles ──────────────────────────────────────────────
    {
        unsigned int seed = p.random_seed;
        init_particles_(x.data(), y.data(), z.data(),
                        vx.data(), vy.data(), vz.data(),
                        mass.data(), radius.data(),
                        &N, &p.init_mode,
                        &p.lx, &p.ly, &p.lz, &seed);
    }

    // ── 5. Open CSV output files ─────────────────────────────────────────────
    ensure_dir("results");

    // kinetic_energy.csv   – one row per output step
    std::ofstream f_ke   ("results/kinetic_energy.csv");
    f_ke << "time,kinetic_energy,n_contacts\n";

    // trajectory.csv       – particle positions/velocities
    std::ofstream f_traj ("results/trajectory.csv");
    f_traj << "time,pid,x,y,z,vx,vy,vz,speed\n";

    // energy_balance.csv   – KE + PE + total
    std::ofstream f_en   ("results/energy_balance.csv");
    f_en << "time,kinetic_energy,potential_energy,total_energy\n";

    // contacts.csv         – contact count time series
    std::ofstream f_con  ("results/contacts.csv");
    f_con << "time,n_contacts\n";

    // per_step_timing.csv  – accumulated wall-time at each output step
    std::ofstream f_wt   ("results/per_step_timing.csv");
    f_wt << "step,time,wall_s\n";

    // ── 6. Timers ────────────────────────────────────────────────────────────
    TimerSet timers;
    timers.add("total");
    timers.add("forces");
    timers.add("integration");
    timers.add("io");

    // ── 7. Derived constants ─────────────────────────────────────────────────
    int    n_steps   = static_cast<int>(std::round(p.t_final / p.dt));
    int    out_every = std::max(1, static_cast<int>(p.output_interval / p.dt));
    int    max_traj  = std::min(N, p.max_traj_particles);
    double g_mag     = std::sqrt(p.gx*p.gx + p.gy*p.gy + p.gz*p.gz);

    std::cout << "  n_steps   = " << n_steps   << "\n";
    std::cout << "  out_every = " << out_every  << " steps\n";
    std::cout << "  traj pids = 0.." << max_traj-1 << "\n\n";

    // ── 8. Main loop ─────────────────────────────────────────────────────────
    timers.start("total");

    for (int step = 0; step <= n_steps; ++step) {

        double time = step * p.dt;

        // ─── a) Zero forces ──────────────────────────────────────────────────
        timers.start("forces");
        zero_forces_(fx.data(), fy.data(), fz.data(), &N);

        // ─── b) Gravity ──────────────────────────────────────────────────────
        add_gravity_(fx.data(), fy.data(), fz.data(),
                     mass.data(), &N,
                     &p.gx, &p.gy, &p.gz);

        // ─── c) Particle–particle contacts ──────────────────────────────────
        int ncontacts = 0;
        compute_particle_contacts_(x.data(), y.data(), z.data(),
                                   vx.data(), vy.data(), vz.data(),
                                   fx.data(), fy.data(), fz.data(),
                                   radius.data(), &N,
                                   &p.kn, &p.gamma_n, &ncontacts);

        // ─── d) Wall contacts ────────────────────────────────────────────────
        compute_wall_contacts_(x.data(), y.data(), z.data(),
                               vx.data(), vy.data(), vz.data(),
                               fx.data(), fy.data(), fz.data(),
                               radius.data(), &N,
                               &p.kn, &p.gamma_n,
                               &p.lx, &p.ly, &p.lz);
        timers.stop("forces");

        // ─── e) Integration ──────────────────────────────────────────────────
        timers.start("integration");
        integrate_particles_(x.data(), y.data(), z.data(),
                             vx.data(), vy.data(), vz.data(),
                             fx.data(), fy.data(), fz.data(),
                             mass.data(), &N, &p.dt);
        timers.stop("integration");

        // ─── f) Output (every out_every steps) ──────────────────────────────
        if (step % out_every == 0) {
            timers.start("io");

            // Kinetic energy (Fortran kernel)
            double ke = 0.0;
            compute_kinetic_energy_(vx.data(), vy.data(), vz.data(),
                                    mass.data(), &N, &ke);

            // Potential energy (C++ loop – serial, cheap)
            double pe = 0.0;
            for (int i = 0; i < N; ++i)
                pe += mass[i] * g_mag * z[i];

            // Write CSV rows
            f_ke  << std::fixed << std::setprecision(8)
                  << time << "," << ke << "," << ncontacts << "\n";

            f_con << time << "," << ncontacts << "\n";

            f_en  << time << "," << ke << "," << pe << ","
                  << (ke + pe) << "\n";

            for (int i = 0; i < max_traj; ++i) {
                double spd = std::sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
                f_traj << std::fixed << std::setprecision(8)
                       << time << "," << i    << ","
                       << x[i] << "," << y[i] << "," << z[i] << ","
                       << vx[i]<< "," << vy[i]<< "," << vz[i]<< ","
                       << spd  << "\n";
            }

            double wall = timers.elapsed("total");
            f_wt << step << "," << time << "," << wall << "\n";

            timers.stop("io");

            // ─── Progress bar ────────────────────────────────────────────────
            int report_interval = std::max(out_every, n_steps / 20);
            if (step % report_interval == 0 || step == n_steps) {
                double pct = 100.0 * step / std::max(1, n_steps);
                int bar = (int)(pct / 5);
                std::cout << "  [";
                for (int b=0;b<20;++b) std::cout << (b<bar ? '#' : '-');
                std::cout << "] " << std::setw(5) << std::fixed
                          << std::setprecision(1) << pct << "%"
                          << "  t=" << std::setprecision(4) << time << "s"
                          << "  KE=" << std::scientific << std::setprecision(3) << ke
                          << "  C=" << ncontacts
                          << "  W=" << std::fixed << std::setprecision(1)
                          << wall << "s\n";
            }
        }
    }

    timers.stop("total");

    // ── 9. Close CSV files ───────────────────────────────────────────────────
    f_ke.close(); f_traj.close(); f_en.close();
    f_con.close(); f_wt.close();

    // ── 10. Timing summary ───────────────────────────────────────────────────
    {
        double tot = timers.elapsed("total");
        std::ofstream ft("results/timing.csv");
        ft << "phase,time_s,fraction\n";
        auto row = [&](const std::string& nm, const std::string& k){
            double t = timers.elapsed(k);
            ft << nm << "," << std::fixed << std::setprecision(6)
               << t << "," << std::setprecision(4) << t/tot << "\n";
        };
        row("total",       "total");
        row("forces",      "forces");
        row("integration", "integration");
        row("io",          "io");
    }

    // ── 11. Run metadata ─────────────────────────────────────────────────────
    {
        std::ofstream fm("results/run_metadata.csv");
        fm << "key,value\n"
           << "n_particles,"   << N                        << "\n"
           << "n_threads,"     << actual_threads           << "\n"
           << "dt,"            << p.dt                     << "\n"
           << "t_final,"       << p.t_final                << "\n"
           << "n_steps,"       << n_steps                  << "\n"
           << "kn,"            << p.kn                     << "\n"
           << "gamma_n,"       << p.gamma_n                << "\n"
           << "init_mode,"     << p.init_mode              << "\n"
           << "total_time_s,"  << timers.elapsed("total")  << "\n"
           << "forces_time_s," << timers.elapsed("forces") << "\n"
           << "out_every,"     << out_every                << "\n";
    }

    // ── 12. Final timing banner ──────────────────────────────────────────────
    double tot = timers.elapsed("total");
    std::cout << "\n══════════════════════════════════\n";
    std::cout <<   "       Timing Summary              \n";
    std::cout <<   "══════════════════════════════════\n";
    auto print_row = [&](const std::string& label, const std::string& key){
        double t = timers.elapsed(key);
        std::cout << "║  " << std::left << std::setw(12) << label << " : "
                  << std::right << std::setw(8) << std::fixed
                  << std::setprecision(3) << t << " s ("
                  << std::setw(5) << std::setprecision(1) << 100.0*t/tot
                  << "%)  ║\n";
    };
    print_row("Total",        "total");
    print_row("Forces",       "forces");
    print_row("Integration",  "integration");
    print_row("I/O",          "io");
    std::cout << "═════════════════════════════════\n";
    std::cout << "\nAll CSV results in: results/\n";
    std::cout << "Package for download: make zip\n\n";

    return 0;
}
