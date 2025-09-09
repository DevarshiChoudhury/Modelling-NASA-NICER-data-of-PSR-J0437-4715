## Initialization

import xpsi
from xpsi.global_imports import gravradius

## Uncomment the following and adjust paths if needed
# import sys
# sys.path.append('./_auto_modules')

PATH_dir = './ST_U_outputs/'
PATH_3C50_NoBKG = '3C50_NoBKG/init/nlive2000_expf1.2_noCONST_noMM_tol0.1'
PATH_3C50_BKG_lower_1sigma = '3C50_BKG_lower_1sigma/init/nlive2000_expf1.2_noCONST_noMM_tol0.1'
PATH_3C50_BKG_lower_2sigma = '3C50_BKG_lower_2sigma/init/nlive2000_expf1.2_noCONST_noMM_tol0.1'
PATH_3C50_BKG_lower_3sigma = '3C50_BKG_lower_3sigma/init/nlive2000_expf1.2_noCONST_noMM_tol0.1'
PATH_3C50_BKG_lower_smooth_3sigma = '3C50_BKG_lower_smooth_3sigma/init/nlive2000_expf1.2_noCONST_noMM_tol0.1'
PATH_3C50_BKG_lower_4sigma = '3C50_BKG_lower_4sigma/init/nlive2000_expf1.2_noCONST_noMM_tol0.1'
PATH_3C50_BKG_DeltaAGN_lower_smooth_3sigma = '3C50_BKG_DeltaAGN_lower_smooth_3sigma/init/nlive2000_expf1.2_noCONST_noMM_tol0.1'
PATH_3C50_BKG_DeltaAGN_smooth_3sigma = '3C50_BKG_DeltaAGN_smooth_3sigma/init/nlive2000_expf1.2_noCONST_noMM_tol0.1'


## Importing all relevant modules for post-processing

from _auto_modules import main_3C50_NoBKG
print('\n3C50_NoBKG log-likelihood value: ', \
    main_3C50_NoBKG.likelihood(main_3C50_NoBKG.p_maxL, reinitialise=True))
print("3C50_NoBKG module loaded \n")

from _auto_modules import main_3C50_BKG_lower_1sigma
print('\n3C50_BKG_lower_1sigma log-likelihood value: ', \
    main_3C50_BKG_lower_1sigma.likelihood(main_3C50_BKG_lower_1sigma.p_maxL, reinitialise=True))
print("3C50_BKG_lower_1sigma module loaded \n")

from _auto_modules import main_3C50_BKG_lower_2sigma
print('\n3C50_BKG_lower_2sigma log-likelihood value: ', \
    main_3C50_BKG_lower_2sigma.likelihood(main_3C50_BKG_lower_2sigma.p_maxL, reinitialise=True))
print("3C50_BKG_lower_2sigma module loaded \n")

from _auto_modules import main_3C50_BKG_lower_3sigma
print('\n3C50_BKG_lower_3sigma log-likelihood value: ', \
    main_3C50_BKG_lower_3sigma.likelihood(main_3C50_BKG_lower_3sigma.p_maxL, reinitialise=True))
print("3C50_BKG_lower_3sigma module loaded \n")

from _auto_modules import main_3C50_BKG_lower_smooth_3sigma
print('\n3C50_BKG_lower_smooth_3sigma log-likelihood value: ', \
    main_3C50_BKG_lower_smooth_3sigma.likelihood(main_3C50_BKG_lower_smooth_3sigma.p_maxL, reinitialise=True))
print("3C50_BKG_lower_smooth_3sigma module loaded \n")

from _auto_modules import main_3C50_BKG_lower_4sigma
print('\n3C50_BKG_lower_4sigma log-likelihood value: ', \
    main_3C50_BKG_lower_4sigma.likelihood(main_3C50_BKG_lower_4sigma.p_maxL, reinitialise=True))
print("3C50_BKG_lower_4sigma module loaded \n")

from _auto_modules import main_3C50_BKG_DeltaAGN_lower_smooth_3sigma
print('\n3C50_BKG_DeltaAGN_lower_smooth_3sigma log-likelihood value: ', \
    main_3C50_BKG_DeltaAGN_lower_smooth_3sigma.likelihood(main_3C50_BKG_DeltaAGN_lower_smooth_3sigma.p_maxL, reinitialise=True))
print("3C50_BKG_DeltaAGN_lower_smooth_3sigma module loaded \n")

from _auto_modules import main_3C50_BKG_DeltaAGN_smooth_3sigma
print('\n3C50_BKG_DeltaAGN_smooth_3sigma log-likelihood value: ', \
    main_3C50_BKG_DeltaAGN_smooth_3sigma.likelihood(main_3C50_BKG_DeltaAGN_smooth_3sigma.p_maxL, reinitialise=True))
print("3C50_BKG_DeltaAGN_smooth_3sigma module loaded \n")

print('All modules loaded \n\nLoading metadata...')


## Specifying parameter bounds and labels for post-processing

# Additional parameters specified in CustomPrior.transfom to be plotted
names = ['compactness', 'N_H_e20_units']

STU_bounds = {'mass': main_3C50_NoBKG.likelihood.get_param('mass').bounds,
              'radius': main_3C50_NoBKG.likelihood.get_param('radius').bounds,
              'cos_inclination': main_3C50_NoBKG.likelihood.get_param('cos_inclination').bounds,
              'p__super_colatitude': main_3C50_NoBKG.likelihood.get_param('p__super_colatitude').bounds,
              'p__super_radius': main_3C50_NoBKG.likelihood.get_param('p__super_radius').bounds,
              'p__super_temperature': main_3C50_NoBKG.likelihood.get_param('p__super_temperature').bounds,
              's__super_colatitude': main_3C50_NoBKG.likelihood.get_param('s__super_colatitude').bounds,
              's__super_radius': main_3C50_NoBKG.likelihood.get_param('s__super_radius').bounds,
              's__super_temperature': main_3C50_NoBKG.likelihood.get_param('s__super_temperature').bounds,
              'XTI__energy_independent_effective_area_scaling_factor': main_3C50_NoBKG.likelihood.get_param('XTI__energy_independent_effective_area_scaling_factor').bounds,
              'neutral_hydrogen_column_density': main_3C50_NoBKG.likelihood.get_param('neutral_hydrogen_column_density').bounds,
              'p__phase_shift': main_3C50_NoBKG.likelihood.get_param('p__phase_shift').bounds,
              's__phase_shift': main_3C50_NoBKG.likelihood.get_param('s__phase_shift').bounds,
              'compactness': (gravradius(1.0)/16.0, 1.0/3.0),
              'N_H_e20_units': (0.04, 2.0)}

# TeX compatible labels
STU_labels = {'mass': r"M\;\mathrm{[M}_{\odot}\mathrm{]}",
              'radius': r"R_{\mathrm{eq}}\;\mathrm{[km]}",
              'cos_inclination': r"\mathrm{cos}(i)",
              'p__super_colatitude': r"\Theta_{p}\;\mathrm{[rad]}",
              'p__super_radius': r"\zeta_{p}\;\mathrm{[rad]}",
              'p__super_temperature': r"\mathrm{log10}(T_{p}\;[\mathrm{K}])",
              's__super_colatitude': r"\Theta_{s}\;\mathrm{[rad]}",
              's__super_radius': r"\zeta_{s}\;\mathrm{[rad]}",
              's__super_temperature': r"\mathrm{log10}(T_{s}\;[\mathrm{K}])",
              'XTI__energy_independent_effective_area_scaling_factor': r"\alpha_{XTI}",
              'neutral_hydrogen_column_density': r"N_{\mathrm{H}}\;\mathrm{[0.4 \times 10^{20}\;cm^{-2}]}",
              'p__phase_shift': r"\phi_{p}\;\mathrm{[cycles]}",
              's__phase_shift': r"\phi_{s}\;\mathrm{[cycles]}",
              'compactness': r"M/R_{\mathrm{eq}}",
              'N_H_e20_units': r"N_{\mathrm{H}}\;\mathrm{[10^{20}\;cm^{-2}]}"}


# KDE settings
getdist_kde_settings = {'ignore_rows': 0,
                         'min_weight_ratio': 1.0e-40,
                         'contours': [0.683, 0.954, 0.997],
                         'credible_interval_threshold': 0.001,
                         'range_ND_contour': 0,
                         'range_confidence': 0.001,
                         'fine_bins': 2048,
                         'smooth_scale_1D': -1,
                         'num_bins': 100,
                         'boundary_correction_order': 1,
                         'mult_bias_correction_order': 1,
                         'smooth_scale_2D': -1,
                         'max_corr_2D': 0.99,
                         'fine_bins_2D': 1024,
                         'num_bins_2D': 40}


## Loading MultiNest output files for post-processing

main_3C50_NoBKG.names = main_3C50_NoBKG.likelihood.names
main_3C50_NoBKG.names += names
main_3C50_NoBKG.bounds = STU_bounds
main_3C50_NoBKG.labels = STU_labels
main_3C50_NoBKG.runs = xpsi.Runs.load_runs(ID = r"3C50, No BKG",
                                           run_IDs = ['run1'],
                                           roots = [PATH_3C50_NoBKG],
                                           base_dirs = [PATH_dir],
                                           use_nestcheck=[True],
                                           kde_settings=getdist_kde_settings,
                                           likelihood = main_3C50_NoBKG.likelihood,
                                           names = main_3C50_NoBKG.names,
                                           bounds = main_3C50_NoBKG.bounds,
                                           labels = main_3C50_NoBKG.labels,
                                           implementation = 'multinest',
                                           overwrite_transformed = True)
print("\n 3C50 NoBKG metadata loaded successfully")

main_3C50_BKG_lower_1sigma.names = main_3C50_BKG_lower_1sigma.likelihood.names
main_3C50_BKG_lower_1sigma.names += names
main_3C50_BKG_lower_1sigma.bounds = STU_bounds
main_3C50_BKG_lower_1sigma.labels = STU_labels
main_3C50_BKG_lower_1sigma.runs = xpsi.Runs.load_runs(ID = r"3C50, BKG-1X",
                                                      run_IDs = ['run1'],
                                                      roots = [PATH_3C50_BKG_lower_1sigma],
                                                      base_dirs = [PATH_dir],
                                                      use_nestcheck=[True],
                                                      kde_settings=getdist_kde_settings,
                                                      likelihood = main_3C50_BKG_lower_1sigma.likelihood,
                                                      names = main_3C50_BKG_lower_1sigma.names,
                                                      bounds = main_3C50_BKG_lower_1sigma.bounds,
                                                      labels = main_3C50_BKG_lower_1sigma.labels,
                                                      implementation = 'multinest',
                                                      overwrite_transformed = True)
print("\n 3C50_BKG_lower_1sigma metadata loaded successfully")

main_3C50_BKG_lower_2sigma.names = main_3C50_BKG_lower_2sigma.likelihood.names
main_3C50_BKG_lower_2sigma.names += names
main_3C50_BKG_lower_2sigma.bounds = STU_bounds
main_3C50_BKG_lower_2sigma.labels = STU_labels
main_3C50_BKG_lower_2sigma.runs = xpsi.Runs.load_runs(ID = r"3C50, BKG-2X",
                                                      run_IDs = ['run1'],
                                                      roots = [PATH_3C50_BKG_lower_2sigma],
                                                      base_dirs = [PATH_dir],
                                                      use_nestcheck=[True],
                                                      kde_settings=getdist_kde_settings,
                                                      likelihood = main_3C50_BKG_lower_2sigma.likelihood,
                                                      names = main_3C50_BKG_lower_2sigma.names,
                                                      bounds = main_3C50_BKG_lower_2sigma.bounds,
                                                      labels = main_3C50_BKG_lower_2sigma.labels,
                                                      implementation = 'multinest',
                                                      overwrite_transformed = True)
print("\n 3C50_BKG_lower_2sigma metadata loaded successfully")

main_3C50_BKG_lower_3sigma.names = main_3C50_BKG_lower_3sigma.likelihood.names
main_3C50_BKG_lower_3sigma.names += names
main_3C50_BKG_lower_3sigma.bounds = STU_bounds
main_3C50_BKG_lower_3sigma.labels = STU_labels
main_3C50_BKG_lower_3sigma.runs = xpsi.Runs.load_runs(ID = r"3C50, BKG-3X",
                                                      run_IDs = ['run1'],
                                                      roots = [PATH_3C50_BKG_lower_3sigma],
                                                      base_dirs = [PATH_dir],
                                                      use_nestcheck=[True],
                                                      kde_settings=getdist_kde_settings,
                                                      likelihood = main_3C50_BKG_lower_3sigma.likelihood,
                                                      names = main_3C50_BKG_lower_3sigma.names,
                                                      bounds = main_3C50_BKG_lower_3sigma.bounds,
                                                      labels = main_3C50_BKG_lower_3sigma.labels,
                                                      implementation = 'multinest',
                                                      overwrite_transformed = True)
print("\n 3C50_BKG_lower_3sigma metadata loaded successfully")

main_3C50_BKG_lower_smooth_3sigma.names = main_3C50_BKG_lower_smooth_3sigma.likelihood.names
main_3C50_BKG_lower_smooth_3sigma.names += names
main_3C50_BKG_lower_smooth_3sigma.bounds = STU_bounds
main_3C50_BKG_lower_smooth_3sigma.labels = STU_labels
main_3C50_BKG_lower_smooth_3sigma.runs = xpsi.Runs.load_runs(ID = r"3C50, BKG-3X (smooth)",
                                                             run_IDs = ['run1'],
                                                             roots = [PATH_3C50_BKG_lower_smooth_3sigma],
                                                             base_dirs = [PATH_dir],
                                                             use_nestcheck=[True],
                                                             kde_settings=getdist_kde_settings,
                                                             likelihood = main_3C50_BKG_lower_smooth_3sigma.likelihood,
                                                             names = main_3C50_BKG_lower_smooth_3sigma.names,
                                                             bounds = main_3C50_BKG_lower_smooth_3sigma.bounds,
                                                             labels = main_3C50_BKG_lower_smooth_3sigma.labels,
                                                             implementation = 'multinest',
                                                             overwrite_transformed = True)
print("\n 3C50_BKG_lower_smooth_3sigma metadata loaded successfully")

main_3C50_BKG_lower_4sigma.names = main_3C50_BKG_lower_4sigma.likelihood.names
main_3C50_BKG_lower_4sigma.names += names
main_3C50_BKG_lower_4sigma.bounds = STU_bounds
main_3C50_BKG_lower_4sigma.labels = STU_labels
main_3C50_BKG_lower_4sigma.runs = xpsi.Runs.load_runs(ID = r"3C50, BKG-4X",
                                                      run_IDs = ['run1'],
                                                      roots = [PATH_3C50_BKG_lower_4sigma],
                                                      base_dirs = [PATH_dir],
                                                      use_nestcheck=[True],
                                                      kde_settings=getdist_kde_settings,
                                                      likelihood = main_3C50_BKG_lower_4sigma.likelihood,
                                                      names = main_3C50_BKG_lower_4sigma.names,
                                                      bounds = main_3C50_BKG_lower_4sigma.bounds,
                                                      labels = main_3C50_BKG_lower_4sigma.labels,
                                                      implementation = 'multinest',
                                                      overwrite_transformed = True)
print("\n 3C50_BKG_lower_4sigma metadata loaded successfully")

main_3C50_BKG_DeltaAGN_lower_smooth_3sigma.names = main_3C50_BKG_DeltaAGN_lower_smooth_3sigma.likelihood.names
main_3C50_BKG_DeltaAGN_lower_smooth_3sigma.names += names
main_3C50_BKG_DeltaAGN_lower_smooth_3sigma.bounds = STU_bounds
main_3C50_BKG_DeltaAGN_lower_smooth_3sigma.labels = STU_labels
main_3C50_BKG_DeltaAGN_lower_smooth_3sigma.runs = xpsi.Runs.load_runs(ID = r"3C50, BKG-3X, $\Delta$-AGN$_{\mathrm{low}}$ (smooth)",
                                                                      run_IDs = ['run1'],
                                                                      roots = [PATH_3C50_BKG_DeltaAGN_lower_smooth_3sigma],
                                                                      base_dirs = [PATH_dir],
                                                                      use_nestcheck=[True],
                                                                      kde_settings=getdist_kde_settings,
                                                                      likelihood = main_3C50_BKG_DeltaAGN_lower_smooth_3sigma.likelihood,
                                                                      names = main_3C50_BKG_DeltaAGN_lower_smooth_3sigma.names,
                                                                      bounds = main_3C50_BKG_DeltaAGN_lower_smooth_3sigma.bounds,
                                                                      labels = main_3C50_BKG_DeltaAGN_lower_smooth_3sigma.labels,
                                                                      implementation = 'multinest',
                                                                      overwrite_transformed = True)
print("\n 3C50_BKG_DeltaAGN_lower_smooth_3sigma metadata loaded successfully")

main_3C50_BKG_DeltaAGN_smooth_3sigma.names = main_3C50_BKG_DeltaAGN_smooth_3sigma.likelihood.names
main_3C50_BKG_DeltaAGN_smooth_3sigma.names += names
main_3C50_BKG_DeltaAGN_smooth_3sigma.bounds = STU_bounds
main_3C50_BKG_DeltaAGN_smooth_3sigma.labels = STU_labels
main_3C50_BKG_DeltaAGN_smooth_3sigma.runs = xpsi.Runs.load_runs(ID = r"3C50, BKG-33, $\Delta$-AGN (smooth)",
                                                                run_IDs = ['run1'],
                                                                roots = [PATH_3C50_BKG_DeltaAGN_smooth_3sigma],
                                                                base_dirs = [PATH_dir],
                                                                use_nestcheck=[True],
                                                                kde_settings=getdist_kde_settings,
                                                                likelihood = main_3C50_BKG_DeltaAGN_smooth_3sigma.likelihood,
                                                                names = main_3C50_BKG_DeltaAGN_smooth_3sigma.names,
                                                                bounds = main_3C50_BKG_DeltaAGN_smooth_3sigma.bounds,
                                                                labels = main_3C50_BKG_DeltaAGN_smooth_3sigma.labels,
                                                                implementation = 'multinest',
                                                                overwrite_transformed = True)
print("\n 3C50_BKG_DeltaAGN_smooth_3sigma metadata loaded successfully")

print('\n \nAll metadata loaded successfully')
