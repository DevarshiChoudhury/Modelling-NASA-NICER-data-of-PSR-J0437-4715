## Initialization

import xpsi
from xpsi.global_imports import gravradius

## Uncomment the following and adjust paths if needed
# import sys
# sys.path.append('./_auto_modules')

PATH_dir = './ST_PDT_outputs/'
PATH_outfile = 'module_generator/8klp_NICER-delta_space-weather/res1/nlive8000_expf10.0_noCONST_noMM_tol0.1'

## Importing all relevant modules for post-processing

from Delta_modules import main
print("ST-PDT module loaded \n")
print('Log-likelihood value: ', main.likelihood())

## Specifying parameter bounds and labels for post-processing

# Additional parameters specified in CustomPrior.transfom to be plotted
names = ['compactness', 'N_H_e20_units']

ST_PDT_bounds = {'mass': main.likelihood.get_param('mass').bounds,
              'radius': main.likelihood.get_param('radius').bounds,
              'cos_inclination': main.likelihood.get_param('cos_inclination').bounds,
              'p__super_colatitude': main.likelihood.get_param('p__super_colatitude').bounds,
              'p__super_radius': main.likelihood.get_param('p__super_radius').bounds,
              'p__super_temperature': main.likelihood.get_param('p__super_temperature').bounds,
              's__super_colatitude': main.likelihood.get_param('s__super_colatitude').bounds,
              's__super_radius': main.likelihood.get_param('s__super_radius').bounds,
              's__cede_colatitude': main.likelihood.get_param('s__cede_colatitude').bounds,
              's__cede_radius': main.likelihood.get_param('s__cede_radius').bounds,
              's__cede_azimuth': main.likelihood.get_param('s__cede_azimuth').bounds,
              's__super_temperature': main.likelihood.get_param('s__super_temperature').bounds,
              's__cede_temperature': main.likelihood.get_param('s__cede_temperature').bounds,
              'XTI__energy_independent_effective_area_scaling_factor': main.likelihood.get_param('XTI__energy_independent_effective_area_scaling_factor').bounds,
              'neutral_hydrogen_column_density': main.likelihood.get_param('neutral_hydrogen_column_density').bounds,
              'p__phase_shift': main.likelihood.get_param('p__phase_shift').bounds,
              's__phase_shift': main.likelihood.get_param('s__phase_shift').bounds,
              'compactness': (gravradius(1.0)/16.0, 1.0/3.0),
              'N_H_e20_units': (0.04, 2.0)}

# TeX compatible labels
ST_PDT_labels = {'mass': r"M\;\mathrm{[M}_{\odot}\mathrm{]}",
              'radius': r"R_{\mathrm{eq}}\;\mathrm{[km]}",
              'cos_inclination': r"cos(i)\;\mathrm{[rad]}",
              'p__super_colatitude': r"\Theta_{p}\;\mathrm{[rad]}",
              'p__super_radius': r"\zeta_{p}\;\mathrm{[rad]}",
              'p__super_temperature': r"\mathrm{log10}(\mathcal{T}_{p}\;[\mathrm{K}])",
              's__super_colatitude': r"\Theta_{s}\;\mathrm{[rad]}",
              's__super_radius': r"\psi_{s}\;\mathrm{[rad]}",
              's__cede_radius': r"\zeta_{s};\mathrm{[rad]}",
              's__cede_colatitude': r"\varkappa_{s}\;\mathrm{[rad]}",
              's__cede_azimuth': r"\varphi_{s}\;\mathrm{[rad]}",
              's__super_temperature': r"\mathrm{log10}(\mathcal{T}_{s}\;[\mathrm{K}])",
              's__cede_temperature': r"\mathrm{log10}(T_{s}\;[\mathrm{K}])",
              'XTI__energy_independent_effective_area_scaling_factor': r"\alpha\;\mathrm{[kpc^{-2}]}",
              'neutral_hydrogen_column_density': r"N_{\mathrm{H}}\;\mathrm{[10^{20}\;cm^{-2}]}",
              'p__phase_shift': r"\phi_{p}\;\mathrm{[cycles]}",
              's__phase_shift': r"\phi_{s}\;\mathrm{[cycles]}",
              'compactness': r"M/R_{\mathrm{eq}}",
              'N_H_e20_units': r"N_{\mathrm{H}}\;\mathrm{[10^{20}\;cm^{-2}]}"}

# KDE settings
getdist_kde_settings = {'ignore_rows': 0,
                         'min_weight_ratio': 1.0e-10,
                         'contours': [0.683, 0.954, 0.997],
                         'credible_interval_threshold': 0.001,
                         'range_ND_contour': 0,
                         'range_confidence': 0.001,
                         'fine_bins': 2048,
                         'smooth_scale_1D': 0.4,
                         'num_bins': 100,
                         'boundary_correction_order': 1,
                         'mult_bias_correction_order': 1,
                         'smooth_scale_2D': 0.4,
                         'max_corr_2D': 0.99,
                         'fine_bins_2D': 512,
                         'num_bins_2D': 40}


## Loading MultiNest output files for post-processing

main.names = main.likelihood.names
main.names += names
main.bounds = ST_PDT_bounds
main.labels = ST_PDT_labels
main.runs = xpsi.Runs.load_runs(ID = r"$\Delta$, SW",
                                 run_IDs = ['run1'],
                                 roots = [PATH_outfile],
                                 base_dirs = [PATH_dir],
                                 use_nestcheck=[True],
                                 kde_settings=getdist_kde_settings,
                                 likelihood = main.likelihood,
                                 names = main.names,
                                 bounds = main.bounds,
                                 labels = main.labels,
                                 implementation = 'multinest',
                                 overwrite_transformed = True
                                 )

print("\n Metadata loaded successfully")
