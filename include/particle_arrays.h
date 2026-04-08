#pragma once
/**
 * data_structures/particle_arrays.h
 *
 * Documents the Structure-of-Arrays (SoA) layout used for particle data.
 *
 * WHY SoA?
 *   - Better cache line utilisation when looping over a single field
 *     (e.g., all x positions, then all forces).
 *   - Required for efficient OpenMP vectorisation (contiguous memory).
 *   - Fortran subroutines receive simple double* pointers – no struct packing.
 *
 * MEMORY LAYOUT
 * ─────────────────────────────────────────────────────────────────────
 *  Array    Size    Unit    Description
 * ─────────────────────────────────────────────────────────────────────
 *  x[N]     N       m       x-position of each particle
 *  y[N]     N       m       y-position
 *  z[N]     N       m       z-position
 *  vx[N]    N      m/s      x-velocity
 *  vy[N]    N      m/s      y-velocity
 *  vz[N]    N      m/s      z-velocity
 *  fx[N]    N       N       x-force accumulator (zeroed each step)
 *  fy[N]    N       N       y-force accumulator
 *  fz[N]    N       N       z-force accumulator
 *  mass[N]  N       kg      particle mass (constant)
 *  radius[N]N       m       particle radius (constant)
 * ─────────────────────────────────────────────────────────────────────
 *
 * All arrays are allocated in main.cpp as std::vector<double> and their
 * raw pointers (.data()) are passed to Fortran kernels.
 *
 * INDEXING: 0-based in C++, 1-based in Fortran (Fortran views these as
 * 1-indexed; C++ passes the base address so particle i in C++ maps to
 * particle i+1 in Fortran).  The Fortran subroutines declare arrays as
 * real(c_double), intent(...) :: arr(n) and loop from 1 to n.
 */

#include <vector>
#include <cstddef>

/**
 * ParticleArrays – thin wrapper for the SoA particle data.
 * Owns all memory; can be passed around or used as a convenience bundle.
 * The Fortran kernels receive raw pointers from .data() calls.
 */
struct ParticleArrays {
    int N;
    // State
    std::vector<double> x, y, z;       // position (m)
    std::vector<double> vx, vy, vz;    // velocity (m/s)
    // Forces (reset each step)
    std::vector<double> fx, fy, fz;    // force accumulator (N)
    // Material properties (constant)
    std::vector<double> mass;          // kg
    std::vector<double> radius;        // m

    explicit ParticleArrays(int n) : N(n),
        x(n), y(n), z(n),
        vx(n), vy(n), vz(n),
        fx(n), fy(n), fz(n),
        mass(n), radius(n) {}
};

/**
 * ContactPair – used for diagnostics / neighbour-list extension.
 * Not used in the core O(N^2) solver but reserved for future cell-list.
 */
struct ContactPair {
    int    i, j;          // particle indices (0-based)
    double overlap;       // delta_ij (m), > 0 means in contact
    double fn;            // normal force magnitude (N)
};

/**
 * WallForce – per-particle wall force record for diagnostics.
 */
struct WallForce {
    int    pid;            // particle index
    int    wall_id;        // 0=xlo, 1=xhi, 2=ylo, 3=yhi, 4=zlo, 5=zhi
    double fn;             // normal force (N)
};
