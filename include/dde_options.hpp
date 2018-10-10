#ifndef DDE_OPTIONS_HPP
#define DDE_OPTIONS_HPP

#include <vector>

extern "C" {

  struct dde_opts_cc {
    double hinit;
    double hmax;
    bool neutral;
    bool track_discontinuities;
    bool interpolation;
    int tracking_level;
    int max_events;
    int max_steps;
    double max_delay;
    int trim_frequency;
  };

  double default_rel_err() {
    return 1.0E-3;
  }

  double default_abs_err() {
    return 1.0E-6;
  }

  struct dde_opts_cc default_dde_opts() {
    struct dde_opts_cc opts;
    opts.hinit = 0.E0;
    opts.hmax = -1.0;
    opts.neutral = false;
    opts.track_discontinuities = true;
    opts.interpolation = false;
    opts.tracking_level = 7;
    opts.max_events = 200;
    opts.max_steps = 10000;
    opts.max_delay = 0.0;
    opts.trim_frequency = 0;
    return opts;
  }
}

#endif
