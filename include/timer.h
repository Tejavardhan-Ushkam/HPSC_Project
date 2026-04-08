#pragma once
#include <chrono>
#include <map>
#include <string>
#include <stdexcept>

/**
 * Simple named-timer set using std::chrono.
 */
class TimerSet {
    using Clock = std::chrono::steady_clock;
    using TP    = std::chrono::time_point<Clock>;

    struct Entry {
        TP     start_tp;
        double total = 0.0;
        bool   running = false;
    };
    std::map<std::string, Entry> timers_;

public:
    void add(const std::string& name) { timers_[name]; }

    void start(const std::string& name) {
        auto& e = timers_.at(name);
        e.start_tp = Clock::now();
        e.running  = true;
    }

    void stop(const std::string& name) {
        auto& e = timers_.at(name);
        if (!e.running) return;
        e.total  += std::chrono::duration<double>(Clock::now() - e.start_tp).count();
        e.running = false;
    }

    double elapsed(const std::string& name) const {
        const auto& e = timers_.at(name);
        return e.total;
    }
};
