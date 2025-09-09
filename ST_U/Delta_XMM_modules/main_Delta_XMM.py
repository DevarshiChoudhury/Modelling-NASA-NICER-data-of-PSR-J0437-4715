""" Main module for NICER x XMM J0437 <- X-PSI v0.7 ST-U. """


import argparse

parser = argparse.ArgumentParser(
    description='''
    Main module for X-PSI ST-U modelling of NICER & XMM J0740+6620 event data.

    You can run this module as a script and launch a sampler, optionally
    with a world of MPI processes.

    Alternate usage: mpiexec -n 4 python -m mpi4py %(prog)s [-h] @<config.ini> [--multinest] [--emcee]

    ''',
    fromfile_prefix_chars='@')

_prfx = 'Absolute or relative path to '
_instr = 'NICER-XTI '
parser.add_argument('--NICER-matrix-path', type=str, help=_prfx + _instr + 'channel-phase count matrix. This path is written to if the file does not exist.')
parser.add_argument('--NICER-event-path', type=str, help=_prfx + _instr + 'event list file.')
parser.add_argument('--NICER-arf-path', type=str, help=_prfx + _instr + 'ARF file.')
parser.add_argument('--NICER-rsp-path', type=str, help=_prfx + _instr + 'RSP file.')
parser.add_argument('--NICER-channels-path', type=str, help=_prfx + _instr + 'channel bounds file.')

_instr = 'XMM-MOS1 '
parser.add_argument('--MOS1-spectrum-path', type=str, help=_prfx + _instr + 'spectrum file.')
#parser.add_argument('--MOS1-arf-path', type=str, help=_prfx + _instr + 'ARF file.')
parser.add_argument('--MOS1-rsp-path', type=str, help=_prfx + _instr + 'RSP file.')
parser.add_argument('--MOS1-channels-path', type=str, help=_prfx + _instr + 'channel bounds file.')
parser.add_argument('--MOS1-background-path', type=str, help=_prfx + _instr + 'background spectrum file.')

_instr = 'XMM-MOS2 '
parser.add_argument('--MOS2-spectrum-path', type=str, help=_prfx + _instr + 'spectrum file.')
#parser.add_argument('--MOS2-arf-path', type=str, help=_prfx + _instr + 'ARF file.')
parser.add_argument('--MOS2-rsp-path', type=str, help=_prfx + _instr + 'RSP file.')
parser.add_argument('--MOS2-channels-path', type=str, help=_prfx + _instr + 'channel bounds file.')
parser.add_argument('--MOS2-background-path', type=str, help=_prfx + _instr + 'background spectrum file.')

parser.add_argument('--attenuation-path', type=str, help=_prfx + 'attenuation file.')
parser.add_argument('--atmosphere-path', type=str, help=_prfx + 'atmosphere file.')

parser.add_argument('--multinest', action='store_true',
                    help='Launch MultiNest sampler. Takes precedence.')
parser.add_argument('--emcee', action='store_true',
                    help='Launch emcee sampler.')
parser.add_argument('--resume', action='store_true',
                    help='Resume sampling?')

parser.add_argument('--NICER', action='store_true',
                    help='Model NICER event data.')
parser.add_argument('--XMM', action='store_true',
                    help='Model XMM event data.')

if __name__ == '__main__':
    args = parser.parse_args()
else:
    args = parser.parse_args(['@./config_Delta.ini','--NICER','--XMM'])

# check if interactive input needed
# NICER
if not args.NICER_matrix_path:
    args.NICER_matrix_path = input('Specify the NICER response matrix path: ')

if not args.NICER_event_path:
    args.NICER_event_path = input('Specify the NICER event file path: ')

if not args.NICER_arf_path:
    args.NICER_arf_path = input('Specify the NICER ARF file path: ')

if not args.NICER_rsp_path:
    args.NICER_rmf_path = input('Specify the NICER RSP file path: ')

if not args.NICER_channels_path:
    args.NICER_channels_path = input('Specify the NICER channel energy file path: ')

# MOS1
if not args.MOS1_spectrum_path:
    args.MOS1_spectrum_path = input('Specify the XMM-MOS1 spectrum file path: ')

if not args.MOS1_rsp_path:
    args.MOS1_rsp_path = input('Specify the XMM-MOS1 RSP file path: ')

if not args.MOS1_channels_path:
    args.MOS1_channels_path = input('Specify the XMM-MOS1 channel energy file path: ')

if not args.MOS1_background_path:
    args.MOS1_background_path = input('Specify the XMM-MOS1 background file path: ')

# MOS2
if not args.MOS2_spectrum_path:
    args.MOS2_spectrum_path = input('Specify the XMM-MOS2 spectrum file path: ')

if not args.MOS2_rsp_path:
    args.MOS2_rsp_path = input('Specify the XMM-MOS2 RSP file path: ')

if not args.MOS2_channels_path:
    args.MOS2_channels_path = input('Specify the XMM-MOS2 channel energy file path: ')

if not args.MOS2_background_path:
    args.MOS2_background_path = input('Specify the XMM-MOS2 background file path: ')

# shared
if not args.attenuation_path:
    args.attenuation_path = input('Specify the attenuation file path: ')

if not args.atmosphere_path:
    args.atmosphere_path = input('Specify the atmosphere file path: ')

import numpy as np
import math

import xpsi
from xpsi.Parameter import Derive

print('Rank reporting: %d' % xpsi._rank)

from xpsi.global_imports import gravradius

from .CustomInstrument_NICER_XMM import CustomInstrument
from .CustomSignal import CustomSignal
from .CustomInterstellar import CustomInterstellar
from .CustomPhotosphere_H import CustomPhotosphere
from .CustomPrior_H_NICER_XMM import CustomPrior

class namespace():
    pass

interstellar = CustomInterstellar.from_SWG(args.attenuation_path,
                                    bounds = dict(column_density = (0.0,10.0)))

signals = [[],]

alpha_bounds = dict(alpha = (0.1, 1.9))

minCH_NICER = 30
maxCH_NICER = 300
minIN_NICER = 2 #0
maxIN_NICER = 1800 #1500
exposure_time_NICER = 2319594.0 #1313137.0


minCH_XMM = 30 #included
maxCH_XMM = 300 #excluded
maxIN_XMM = 700 #excluded
minIN_XMM = 9 #included

exposure_time_MOS1 = 1.824928E+05
exposure_time_MOS2 = 1.828327E+05
exposure_time_BKG_MOS1 = 1.824978E+05
exposure_time_BKG_MOS2 = 1.828314E+05

#BACKSCAL_SOURCE/BACKSCAL_BKG:
BACKSCAl_MOS1 =  1.0/1.109E+01
BACKSCAL_MOS2 =  1.0/1.170E+01

#-------#
# NICER #
#-------#
if args.NICER:
    NICER = namespace()

    try:
        counts = np.loadtxt(args.NICER_matrix_path, dtype=np.double)
    except IOError:
        NICER.data = xpsi.Data.phase_bin__event_list(args.NICER_event_path,
                                              channels=np.arange(minCH_NICER, maxCH_NICER),
                                              phases=np.linspace(0.0, 1.0, 33),
                                              channel_column=1,
                                              phase_column=0,
                                              skiprows=3,
                                              dtype=np.double,
                                              first=0,
                                              last=maxCH_NICER-minCH_NICER-1,
                                              exposure_time=exposure_time_NICER)

        np.savetxt(args.NICER_matrix_path, NICER.data.counts)
    else:
        NICER.data = xpsi.Data(counts,
                               channels=np.arange(minCH_NICER, maxCH_NICER),
                               phases=np.linspace(0.0, 1.0, 33),
                               first=0,
                               last=maxCH_NICER-minCH_NICER-1,
                               exposure_time=exposure_time_NICER)

    NICER.instrument = CustomInstrument.NICER_XTI(bounds = alpha_bounds,
                                          values = {},
                                          ARF = args.NICER_arf_path,
                                          RSP = args.NICER_rsp_path,
                                          max_input = maxIN_NICER,
                                          max_channel = maxCH_NICER,
                                          min_input = minIN_NICER,
                                          min_channel = minCH_NICER,
                                          channel_edges = args.NICER_channels_path,
                                          prefix = 'XTI')

    NICER.signal = CustomSignal(data = NICER.data,
                                  instrument = NICER.instrument,
                                  interstellar = interstellar,
                                  cache = False,
                                  workspace_intervals = 1000,
                                  epsrel = 1.0e-8,
                                  epsilon = 1.0e-3,
                                  sigmas = 10.0)

    signals[0].append(NICER.signal)

#-----#
# XMM #
#-----#
if args.XMM:

    #----------#
    # XMM-MOS1 #
    #----------#
    MOS1 = namespace()


    MOS1.instrument = CustomInstrument.XMM_MOS1(bounds = alpha_bounds,
                                            values = {},
                                            #ARF = args.MOS1_arf_path,
                                            RSP = args.MOS1_rsp_path,
                                            max_input = maxIN_XMM,
                                            max_channel = maxCH_XMM,
                                            min_input = minIN_XMM, # skip intervals 7 & 8
                                            min_channel = minCH_XMM,
                                            channel_edges = args.MOS1_channels_path,
                                            prefix = 'MOS1')
    spectrum_data = np.loadtxt(args.MOS1_spectrum_path,
                          skiprows=3,
                          usecols=1,
                          dtype=np.double)[minCH_XMM:maxCH_XMM]
    
    MOS1.data = xpsi.Data(spectrum_data.reshape(-1,1),
                          channels=MOS1.instrument.channels,
                          phases=np.array([0.0, 1.0]),
                          first=0,
                          last=len(MOS1.instrument.channels) - 1,
                          exposure_time=exposure_time_MOS1) # 1.795957421875e4

    spectrum = np.loadtxt(args.MOS1_background_path,
                          skiprows=3,
                          usecols=1,
                          dtype=np.double)[minCH_XMM:maxCH_XMM]
                          
    support = np.zeros((len(spectrum), 2), dtype=np.double)
    support[:,0] = spectrum - 4.0 * np.sqrt(spectrum)
    support[support[:,0] < 0.0, 0] = 0.0
    support[:,1] = spectrum + 4.0 * np.sqrt(spectrum)

    for i in range(support.shape[0]):
        if support[i,1] == 0.0:
            for j in range(i, support.shape[0]):
                if support[j,1] > 0.0:
                    support[i,0] = support[j,1]
                    break

    support *= BACKSCAl_MOS1 * (MOS1.data.exposure_time / exposure_time_BKG_MOS1) # BACKSCAL x exposure ratio # 1.109E+01  1.824978E+05 #1.074 1.57623e6

    support /= MOS1.data.exposure_time # need count rate, so divide by exposure time

    MOS1.signal = CustomSignal(data = MOS1.data,
                               instrument = MOS1.instrument,
                               interstellar = interstellar,
                               support = support,
                               cache = False,
                               workspace_intervals = 1000,
                               epsrel = 1.0e-8,
                               epsilon = 1.0e-3,
                               sigmas = 10.0)

    signals[0].append(MOS1.signal)

    #----------#
    # XMM-MOS2 #
    #----------#
    MOS2 = namespace()
    
    class derive(Derive):
        global MOS1
    
        def __init__(self):
            pass
    
        def __call__(self, boundto, caller=None):
            return MOS1.instrument['alpha']
            
    MOS2.instrument = CustomInstrument.XMM_MOS2(bounds = {'alpha': None},
                                            values = {'alpha': derive()},
                                            RSP = args.MOS2_rsp_path,
                                            max_input = maxIN_XMM,
                                            max_channel = maxCH_XMM,
                                            min_input = minIN_XMM, # skip intervals 7 & 8
                                            min_channel = minCH_XMM,
                                            channel_edges = args.MOS2_channels_path,
                                            prefix = 'MOS2')
    spectrum_data = np.loadtxt(args.MOS2_spectrum_path,
                          skiprows=3,
                          usecols=1,
                          dtype=np.double)[minCH_XMM:maxCH_XMM]
    
    MOS2.data = xpsi.Data(spectrum_data.reshape(-1,1),
                          channels=MOS2.instrument.channels,
                          phases=np.array([0.0, 1.0]),
                          first=0,
                          last=len(MOS2.instrument.channels) - 1,
                          exposure_time=exposure_time_MOS2) #1.828327E+05  1.8680734375e4

    spectrum = np.loadtxt(args.MOS2_background_path,
                          skiprows=3,
                          usecols=1,
                          dtype=np.double)[minCH_XMM:maxCH_XMM]

    support = np.zeros((len(spectrum), 2), dtype=np.double)
    support[:,0] = spectrum - 4.0 * np.sqrt(spectrum)
    support[support[:,0] < 0.0, 0] = 0.0
    support[:,1] = spectrum + 4.0 * np.sqrt(spectrum)

    for i in range(support.shape[0]):
        if support[i,1] == 0.0:
            for j in range(i, support.shape[0]):
                if support[j,1] > 0.0:
                    support[i,0] = support[j,1]
                    break

    support *= BACKSCAL_MOS2* (MOS2.data.exposure_time / exposure_time_BKG_MOS2) # BACKSCAL x exposure ratio  1.170E+01   1.828314E+05 #1.260 1.51256e6

    support /= MOS2.data.exposure_time # need count rate, so divide by exposure time

    MOS2.signal = CustomSignal(data = MOS2.data,
                               instrument = MOS2.instrument,
                               interstellar = interstellar,
                               support = support,
                               cache = False,
                               workspace_intervals = 1000,
                               epsrel = 1.0e-8,
                               epsilon = 1.0e-3,
                               sigmas = 10.0)

    signals[0].append(MOS2.signal)

#-------#
# J0437 #
#-------#

bounds = dict(mass = (None, None),
              radius = (3.0*gravradius(1.0), 16.0),
              cos_inclination = (None, None))

spacetime = xpsi.Spacetime(bounds, dict(frequency =  173.68795,
                                        distance = 0.15669))

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
                            image_order_limit=3,
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

likelihood = xpsi.Likelihood(star = star, signals = signals,
                             num_energies = 128,
                             threads = 1,
                             externally_updated = True,
                             prior = prior)
print (likelihood)
# remember: get a point to check the likelihood between before using cluster
# or supercomp time
if args.XMM and args.NICER:
    p =[1.4949909256216776, 6.701010757850769, -0.7348265225366337, 0.7484694308107467, 0.49418992687323, 0.5253861962180668, 5.121164939167022, 0.48877554050290617, 1.888875034679965, 0.11794663371491185, 5.37472766374477, 0.9451158417416874, 1.1745594609294474, 0.9378078736228744]
    #likelihood.check(None, [-242075.635], 1.0e-6,physical_points=[p])
    #-70754.69408
    #-2.42075635e+05
    #exit()
    
elif args.XMM:
    p =[1.4949909256216776, 6.701010757850769, -0.7348265225366337, 0.7484694308107467, 0.49418992687323, 0.5253861962180668, 5.121164939167022, 0.48877554050290617, 1.888875034679965, 0.11794663371491185, 5.37472766374477, 0.9378078736228744, 1.1745594609294474]
    #likelihood.check(None, [-7164.70853764], 1.0e-6,
    #                 physical_points=[p])
elif args.NICER:
    p =[1.4949909256216776, 6.701010757850769, -0.7348265225366337, 0.7484694308107467, 0.49418992687323, 0.5253861962180668, 5.121164939167022, 0.48877554050290617, 1.888875034679965, 0.11794663371491185, 5.37472766374477, 0.9451158417416874, 1.1745594609294474]
    #print (likelihood(p,reinitialise=True))
    #exit()
    likelihood.check(None, [-63589.9855457], 1.0e-6,
                     physical_points=[p])


if __name__ != '__main__':
    names = ['mass', 'radius', 'cos_inclination', 
             'p__phase_shift', 'p__super_colatitude', 'p__super_radius', 'p__super_temperature',
             's__phase_shift', 's__super_colatitude', 's__super_radius', 's__super_temperature',
             'XTI__alpha',
             'column_density',
             'MOS1__alpha',]

    # names of derived variables of interest
    names += ['compactness', 'p__phase_shift_shifted', 's__phase_shift_shifted']

    bounds = {'mass': (2.131 - 10.0 * 0.0724, 2.131 + 10.0 * 0.0724, ),
              'radius': (3.0 * gravradius(1.0), 16.0),
              'cos_inclination': (0.0442 - 10.0 * 0.00308, 0.0442 + 10.0 * 0.00308),
              'p__super_colatitude': (0.001, math.pi - 0.001),
              'p__super_radius': (0.001, math.pi/2.0 - 0.001),
              'p__super_temperature': (5.1, 6.8),
              's__super_colatitude': (0.001, math.pi - 0.001),
              's__super_radius': (0.001, math.pi/2.0 - 0.001),
              's__super_temperature': (5.1, 6.8),
              'column_density': (0.0, 10.0),
              'p__phase_shift': (-0.25,0.75),
              's__phase_shift': (-0.25,0.75),
              'p__phase_shift_shifted': (-0.5,0.5),
              's__phase_shift_shifted': (-0.5,0.5),
              'XTI__alpha': (0.1, 1.9),
              'MOS1__alpha': (0.1, 1.9),
              'compactness': (gravradius(1.0)/16.0, 1.0/3.0)}

    # TeX compatible labels
    labels = {'mass': r"M\;\mathrm{[M}_{\odot}\mathrm{]}",
              'radius': r"R_{\mathrm{eq}}\;\mathrm{[km]}",
              'cos_inclination': r"\cos(i)",
              'p__super_colatitude': r"\Theta_{p}\;\mathrm{[rad]}",
              'p__super_radius': r"\zeta_{p}\;\mathrm{[rad]}",
              'p__super_temperature': r"\mathrm{log10}(T_{p}\;[\mathrm{K}])",
              's__super_colatitude': r"\Theta_{s}\;\mathrm{[rad]}",
              's__super_radius': r"\zeta_{s}\;\mathrm{[rad]}",
              's__super_temperature': r"\mathrm{log10}(T_{s}\;[\mathrm{K}])",
              'column_density': r"N_{\mathrm{H}}\;\mathrm{[10^{20}\;cm^{-2}]}",
              'XTI__alpha': r"\alpha_{\rm XTI}",
              'MOS1__alpha': r"\alpha_{\rm mos1}",
              'p__phase_shift': r"\phi_{p}\;\mathrm{[cycles]}",
              's__phase_shift': r"\phi_{s}\;\mathrm{[cycles]}",
              'compactness': r"M/R_{\mathrm{eq}}",
              'p__phase_shift_shifted': r"\phi_{p}\;\mathrm{[cycles]}",
              's__phase_shift_shifted': r"\phi_{s}\;\mathrm{[cycles]}",}

else:
    wrapped_params = [0] * len(likelihood)
    wrapped_params[likelihood.index('p__phase_shift')] = 1
    wrapped_params[likelihood.index('s__phase_shift')] = 1

    runtime_params = {'resume': False,
                      'importance_nested_sampling': False,
                      'multimodal': False,
                      'n_clustering_params': None,
                      'outputfiles_basename': './stu_NICERxXMM_run1_nlive4000_eff0.3_noCONST_noMM_noIS_tol-1',
                      'n_iter_before_update': 100,
                      'n_live_points': 4000,
                      'sampling_efficiency': 0.3,
                      'const_efficiency_mode': False,
                      'wrapped_params': wrapped_params,
                      'evidence_tolerance': 0.1,
                      'max_iter': -1,
                      'verbose': True}

    xpsi.Sample.nested(likelihood, prior, **runtime_params)
