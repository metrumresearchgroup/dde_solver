#ifndef DDE_OPTIONS_HPP
#define DDE_OPTIONS_HPP

#include <vector>

extern "C" {

  struct dde_opts_cc {
    double hinit;               // initial step size
    double hmax;                // max step size
    bool neutral;               // solve DDE of neutral type
    bool track_discontinuities; // enable tracking of discontinuities
    bool interpolation;         // interpolate solutions
    int tracking_level;         // max discontinuities tracking level
    int max_events;             // Maximum allowable number of event occurrences
    int max_steps;              // Maximum allowable number of integration steps.
    double max_delay;           // upper bound of delay
    int trim_frequency;         // Frequency of trimming solution queue
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
