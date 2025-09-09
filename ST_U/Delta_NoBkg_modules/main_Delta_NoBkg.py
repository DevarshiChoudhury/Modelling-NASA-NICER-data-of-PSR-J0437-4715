""" Main module for NICER J0437 <- X-PSI v0.7.0 ST-U. """


import argparse

parser = argparse.ArgumentParser(
    description='''
    Main module for X-PSI ST-U modelling of NICER J0437-4715 event data.

    You can run this module as a script and launch a sampler, optionally
    with a world of MPI processes.

    Alternate usage: mpiexec -n 4 python -m mpi4py %(prog)s [-h] @<config.ini> [--multinest] [--emcee]

    ''',
    fromfile_prefix_chars='@')

_prfx = 'Absolute or relative path to '
parser.add_argument('--matrix-path', type=str, help=_prfx + 'response matrix file.')
parser.add_argument('--event-path', type=str, help=_prfx + 'event list file.')
parser.add_argument('--arf-path', type=str, help=_prfx + 'ARF file.')
parser.add_argument('--rsp-path', type=str, help=_prfx + 'RSP file.')
parser.add_argument('--channels-path', type=str, help=_prfx + 'channel bounds file.')
parser.add_argument('--attenuation-path', type=str, help=_prfx + 'attenuation file.')
parser.add_argument('--atmosphere-path', type=str, help=_prfx + 'atmosphere file.')

parser.add_argument('--multinest', action='store_true',
                    help='Launch MultiNest sampler. Takes precedence.')
parser.add_argument('--emcee', action='store_true',
                    help='Launch emcee sampler.')

if __name__ == '__main__':
    args = parser.parse_args()
else:
    args = parser.parse_args(['@./config_Delta_NoBkg.ini'])

# check if interactive input needed
if not args.matrix_path:
    args.matrix_path = input('Specify the response matrix path: ')

if not args.event_path:
    event_path = input('Specify the event file path: ')

if not args.arf_path:
    args.arf_path = input('Specify the ARF file path: ')

if not args.rsp_path:
    args.rsp_path = input('Specify the RSP file path: ')

if not args.channels_path:
    args.channels_path = input('Specify the channel energy file path: ')

if not args.attenuation_path:
    args.attenuation_path = input('Specify the attenuation file path: ')

if not args.atmosphere_path:
    args.atmosphere_path = input('Specify the atmosphere file path: ')

import numpy as np
import math

import xpsi

print('Rank reporting: %d' % xpsi._rank)

from xpsi.global_imports import gravradius

from .CustomInstrument import CustomInstrument
from .CustomInterstellar import CustomInterstellar
from .CustomSignal import CustomSignal
from .CustomPrior_H import CustomPrior
from .CustomPhotosphere_H import CustomPhotosphere

minCH = 30
maxCH = 300
minEI = 2
maxEI = 1800

try:
    counts = np.loadtxt(args.matrix_path, dtype=np.double)
except IOError:
    data = xpsi.Data.phase_bin__event_list(args.event_path,
                                           channels=np.arange(minCH, maxCH),
                                           phases=np.linspace(0.0, 1.0, 33),
                                           phase_column=0,
                                           channel_column=1,
                                           skiprows=3,
                                           dtype=np.double,
                                           first=0,
                                           last=maxCH-minCH-1,
                                           exposure_time=2319594.0)

    np.savetxt(args.matrix_path, data.counts)
else:
    data = xpsi.Data(counts,
                     channels=np.arange(minCH, maxCH),
                     phases=np.linspace(0.0, 1.0, 33),
                     first=0,
                     last=maxCH-minCH-1,
                     exposure_time=2319594.0)

NICER = CustomInstrument.from_SWG(bounds = dict(beta = (None, None)),
                                  values = {},
                                  ARF = args.arf_path,
                                  RSP = args.rsp_path,
                                  max_input = maxEI,
                                  min_input = minEI,
                                  channel_edges = args.channels_path)

interstellar = CustomInterstellar.from_SWG(args.attenuation_path,
                                           bounds = dict(column_density = (0.0,5.0)))

signal = CustomSignal(data = data,
                      instrument = NICER,
                      interstellar = interstellar,
                      cache = True,
                      workspace_intervals = 1000,
                      epsrel = 1.0e-8,
                      epsilon = 1.0e-3,
                      sigmas = 10.0)

bounds = dict(mass = (None, None), # Gaussian prior
              radius = (3.0*gravradius(1.0), 16.0),
              cos_inclination = (None, None)) # Gaussian prior

spacetime = xpsi.Spacetime(bounds, dict(frequency = 173.68795,
                                        distance = 0.01)) # fixed dummy distance

bounds = dict(super_colatitude = (0.001, math.pi - 0.001),
              super_radius = (0.001, math.pi/2.0 - 0.001),
              phase_shift = (-0.25, 0.75),
              super_temperature = (5.1, 6.8))

primary = xpsi.HotRegion(bounds=bounds,
                            values={},
                            symmetry=True,
                            omit=False,
                            cede=False,
                            concentric=False,
                            sqrt_num_cells=32,
                            min_sqrt_num_cells=16,
                            max_sqrt_num_cells=64,
                            num_leaves=100,
                            num_rays=512,
                            is_secondary=False,
                            image_order_limit=3, # up to tertiary
                            atm_ext="Num4D",
                            prefix='p')

bounds = dict(super_colatitude = (0.001, math.pi - 0.001),
              super_radius = (0.001, math.pi/2.0 - 0.001),
              phase_shift = (-0.25, 0.75),
              super_temperature = (5.1, 6.8))

secondary = xpsi.HotRegion(bounds=bounds,
                            values={},
                            symmetry=True,
                            omit=False,
                            cede=False,
                            concentric=False,
                            sqrt_num_cells=32,
                            min_sqrt_num_cells=16,
                            max_sqrt_num_cells=64,
                            num_leaves=100,
                            num_rays=512,
                            is_antiphased=True,
                            image_order_limit=3,
                            atm_ext="Num4D",
                            prefix='s')

from xpsi import HotRegions

hot = HotRegions((primary, secondary))

photosphere = CustomPhotosphere(hot = hot, elsewhere = None,
                                values=dict(mode_frequency = spacetime['frequency']))
photosphere.hot_atmosphere = args.atmosphere_path

star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

prior = CustomPrior()

likelihood = xpsi.Likelihood(star = star, signals = signal,
                             num_energies = 128,
                             threads = 1,
                             externally_updated = True,
                             prior = prior)

# likelihood.prior = prior

p = [0.146458802770893226E+01,
     0.102924085054358176E+02,
    -0.737823824628629055E+00,
     0.130093166545462435E-01,
     0.203641608954142794E+00,
     0.330910964933071206E+00,
     0.615719910407865800E+01,
    -0.931518131300512198E-01,
     0.312920098722429929E+00,
     0.894553798177937787E-01,
     0.623599843093417672E+01,
     0.342918838503153296E+02,
     0.178524587594019735E-01]

#print (likelihood(p,reinitialise=True))
#exit()
likelihood.check(None, [-38447.2695774], 1.0e-5,
                 physical_points=[p])

if __name__ == '__main__': # sample from the posterior
    # transform relevant input information below to conmmand line arguments
    # and config file arguments

    wrapped_params = [0] * len(likelihood)
    wrapped_params[likelihood.index('p__phase_shift')] = 1
    wrapped_params[likelihood.index('s__phase_shift')] = 1

    runtime_params = {'resume': False,
                      'importance_nested_sampling': False,
                      'multimodal': False,
                      'n_clustering_params': None,
                      'outputfiles_basename': './stu_run_DELTA_data_set_nlive4000_eff0.3_noCONST_noMM_noIS_tol-1',
                      'n_iter_before_update': 100,
                      'n_live_points': 4000,
                      'sampling_efficiency': 0.3,
                      'const_efficiency_mode': False,
                      'wrapped_params': wrapped_params,
                      'evidence_tolerance': 0.1,
                      'max_iter': -1,
                      'verbose': True}

    xpsi.Sample.nested(likelihood, prior, **runtime_params)
