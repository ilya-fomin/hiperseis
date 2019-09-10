#!/usr/bin/env python

import logging

import numpy as np

import rf

from seismic.receiver_fn.plot_ccp import process_plot_ccp

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

rf_file = r"C:\software\hiperseis\seismic\receiver_fn\DATA\OA-ZRT-R-cleaned.h5"
stream = rf.read_rf(rf_file, 'H5')

outfile_base = "frame_NS_{}.png"
width = 80.0
spacing = 2.0
max_depth = 110.0
stacked_scale = 0.2
channels = 'R'
title = 'OA W->E lines sweeping N->S'
plot_rf = False

# Generate North to South animation along West-East lines
try:
    for i, lat in enumerate(np.linspace(-17.5, -21.5, 9)):
        output_file = outfile_base.format(i)
        start_latlon = (lat, 132.8)
        end_latlon = (lat, 141.2)
        process_plot_ccp(stream, output_file, start_latlon, end_latlon, width, spacing, max_depth,
                         stacked_scale, channels, title, plot_rf, log=log, figsize=(16, 9))
    # end for
except KeyboardInterrupt:
    log.error("Interrupted N-S frame generation")
# end try

outfile_base = "frame_WE_{}.png"
title = 'OA N->S lines sweeping W->E'
try:
    # Generate West to East animation along North-South lines
    for i, lon in enumerate(np.linspace(133.0, 141.0, 17)):
        output_file = outfile_base.format(i)
        start_latlon = (-17.3, lon)
        end_latlon = (-21.7, lon)
        process_plot_ccp(stream, output_file, start_latlon, end_latlon, width, spacing, max_depth,
                         stacked_scale, channels, title, plot_rf, log=log, figsize=(9, 9))
    # end for
except KeyboardInterrupt:
    log.error("Interrupted W-E frame generation")
# end try
