#pragma once
#include <string>
#include <fstream>
#include <iostream>

/**
 * All simulation parameters in one struct.
 * Loaded from a simple key=value config file or left at defaults.
 */
struct SimParams {
    // Domain
    double lx = 1.0, ly = 1.0, lz = 1.0;

    // Gravity (m/s²)
    double gx = 0.0, gy = 0.0, gz = -9.81;

    // Contact model (linear spring-dashpot)
    double kn      = 1.0e5;   // Normal stiffness (N/m)
    double gamma_n = 50.0;    // Normal damping   (N·s/m)

    // Time integration
    double dt      = 1.0e-5;
    double t_final = 0.5;

    // Particle defaults
    int    n_particles        = 200;
    double particle_r         = 0.02;
    double particle_m         = 1.0;

    // OpenMP
    int n_threads = 4;

    // Output control
    double output_interval    = 1.0e-3;  // seconds between output writes
    int    max_traj_particles = 500;     // max particles dumped in trajectory.csv

    // Initialisation mode
    //  1 = random cloud    (multi-particle default)
    //  2 = free-fall test  (single particle, z0=0.9*Lz, v=0)
    //  3 = bounce test     (single particle, z0=0.8*Lz, v=0)
    //  4 = const-vel test  (single particle, vz=2.0, g=0)
    int          init_mode   = 1;
    unsigned int random_seed = 42;
};

inline SimParams default_params() { return SimParams{}; }

// Simple key=value parser (lines starting with # are comments)
inline void load_params(const std::string& fname, SimParams& p) {
    std::ifstream f(fname);
    if (!f.is_open()) {
        std::cerr << "[cfg] Cannot open '" << fname << "' – using defaults.\n";
        return;
    }

    auto trim = [](std::string& s) {
        size_t l = s.find_first_not_of(" \t\r\n");
        size_t r = s.find_last_not_of (" \t\r\n");
        s = (l == std::string::npos) ? "" : s.substr(l, r - l + 1);
    };

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string key = line.substr(0, eq);
        std::string val = line.substr(eq + 1);
        trim(key); trim(val);

        try {
            if      (key == "lx")                   p.lx                   = std::stod(val);
            else if (key == "ly")                   p.ly                   = std::stod(val);
            else if (key == "lz")                   p.lz                   = std::stod(val);
            else if (key == "gx")                   p.gx                   = std::stod(val);
            else if (key == "gy")                   p.gy                   = std::stod(val);
            else if (key == "gz")                   p.gz                   = std::stod(val);
            else if (key == "kn")                   p.kn                   = std::stod(val);
            else if (key == "gamma_n")              p.gamma_n              = std::stod(val);
            else if (key == "dt")                   p.dt                   = std::stod(val);
            else if (key == "t_final")              p.t_final              = std::stod(val);
            else if (key == "n_particles")          p.n_particles          = std::stoi(val);
            else if (key == "particle_r")           p.particle_r           = std::stod(val);
            else if (key == "particle_m")           p.particle_m           = std::stod(val);
            else if (key == "n_threads")            p.n_threads            = std::stoi(val);
            else if (key == "output_interval")      p.output_interval      = std::stod(val);
            else if (key == "max_traj_particles")   p.max_traj_particles   = std::stoi(val);
            else if (key == "init_mode")            p.init_mode            = std::stoi(val);
            else if (key == "random_seed")          p.random_seed          = (unsigned)std::stoul(val);
            else std::cerr << "[cfg] Unknown key ignored: '" << key << "'\n";
        } catch (...) {
            std::cerr << "[cfg] Bad value for key '" << key << "': " << val << "\n";
        }
    }
    std::cout << "[cfg] Loaded: " << fname << "\n";
}
