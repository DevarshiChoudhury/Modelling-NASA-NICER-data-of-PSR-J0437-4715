## Initialization
import xpsi
from xpsi.global_imports import gravradius

## Uncomment the following and adjust paths if needed
# import sys
# sys.path.append('./_auto_modules')

PATH_dir = './ST_U_outputs/'
PATH_Delta_NoBkg = 'Delta_NoBKG/stu_run_DELTA_data_set_nlive4000_eff0.3_noCONST_noMM_noIS_tol-1'
PATH_Delta_elsewhere_SW = 'Delta_elsewhere_SW/stu_NICER_run1_nlive4000_eff0.3_noCONST_noMM_noIS_tol-1'
PATH_Delta_XMM = 'Delta_XMM/stu_NICERxXMM_run1_nlive4000_eff0.3_noCONST_noMM_noIS_tol-1'


## Importing all relevant modules for post-processing

from Delta_NoBkg_modules import main_Delta_NoBkg
from Delta_elsewhere_SW_modules import main_Delta_elsewhere_SW
from Delta_XMM_modules import main_Delta_XMM

print('Delta modules loaded \n\nLoading metadata...')

## Specifying parameter bounds and labels for post-processing

# Additional parameters specified in CustomPrior.transfom to be plotted
names = ['compactness']

STU_bounds = {'mass': main_Delta_NoBkg.likelihood.get_param('mass').bounds,
              'radius': main_Delta_NoBkg.likelihood.get_param('radius').bounds,
              'cos_inclination': main_Delta_NoBkg.likelihood.get_param('cos_inclination').bounds,
              'p__super_colatitude': main_Delta_NoBkg.likelihood.get_param('p__super_colatitude').bounds,
              'p__super_radius': main_Delta_NoBkg.likelihood.get_param('p__super_radius').bounds,
              'p__super_temperature': main_Delta_NoBkg.likelihood.get_param('p__super_temperature').bounds,
              's__super_colatitude': main_Delta_NoBkg.likelihood.get_param('s__super_colatitude').bounds,
              's__super_radius': main_Delta_NoBkg.likelihood.get_param('s__super_radius').bounds,
              's__super_temperature': main_Delta_NoBkg.likelihood.get_param('s__super_temperature').bounds,
              'beta': main_Delta_NoBkg.likelihood.get_param('beta').bounds,
              'column_density': main_Delta_NoBkg.likelihood.get_param('column_density').bounds,
              'p__phase_shift': main_Delta_NoBkg.likelihood.get_param('p__phase_shift').bounds,
              's__phase_shift': main_Delta_NoBkg.likelihood.get_param('s__phase_shift').bounds,
              'compactness': (gravradius(1.0)/16.0, 1.0/3.0)}

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
              'beta': r"\beta",
              'column_density': r"N_{\mathrm{H}}\;\mathrm{[10^{20}\;cm^{-2}]}",
              'p__phase_shift': r"\phi_{p}\;\mathrm{[cycles]}",
              's__phase_shift': r"\phi_{s}\;\mathrm{[cycles]}",
              'compactness': r"M/R_{\mathrm{eq}}"}

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

main_Delta_NoBkg.names = main_Delta_NoBkg.likelihood.names
main_Delta_NoBkg.names += names
main_Delta_NoBkg.bounds = STU_bounds
main_Delta_NoBkg.labels = STU_labels
main_Delta_NoBkg.runs = xpsi.Runs.load_runs(ID = r"$\Delta$, No BKG",
                                           run_IDs = ['run1'],
                                           roots = [PATH_Delta_NoBkg],
                                           base_dirs = [PATH_dir],
                                           use_nestcheck=[True],
                                           kde_settings=getdist_kde_settings,
                                           likelihood = main_Delta_NoBkg.likelihood,
                                           names = main_Delta_NoBkg.names,
                                           bounds = main_Delta_NoBkg.bounds,
                                           labels = main_Delta_NoBkg.labels,
                                           implementation = 'multinest',
                                           overwrite_transformed = True)
print("\n Delta NoBKG metadata loaded successfully")

del STU_bounds['beta']
del STU_labels['beta']
STU_bounds['XTI__alpha']= main_Delta_elsewhere_SW.likelihood.get_param('XTI__alpha').bounds
STU_labels['XTI__alpha'] = r"\alpha_{XTI}"
STU_bounds['elsewhere_temperature'] = (5.0, 6.5)
STU_labels['elsewhere_temperature'] = r"\mathrm{log10}(T_{else}\;[\mathrm{K}])"

main_Delta_elsewhere_SW.names = main_Delta_elsewhere_SW.likelihood.names
main_Delta_elsewhere_SW.names += names
main_Delta_elsewhere_SW.bounds = STU_bounds
main_Delta_elsewhere_SW.labels = STU_labels
main_Delta_elsewhere_SW.runs = xpsi.Runs.load_runs(ID = r"$\Delta$, SW, elsewhere",
                                           run_IDs = ['run1'],
                                           roots = [PATH_Delta_elsewhere_SW],
                                           base_dirs = [PATH_dir],
                                           use_nestcheck=[True],
                                           kde_settings=getdist_kde_settings,
                                           likelihood = main_Delta_elsewhere_SW.likelihood,
                                           names = main_Delta_elsewhere_SW.names,
                                           bounds = main_Delta_elsewhere_SW.bounds,
                                           labels = main_Delta_elsewhere_SW.labels,
                                           implementation = 'multinest',
                                           overwrite_transformed = True)
print("\n Delta elsewhere SW metadata loaded successfully")

del STU_bounds['elsewhere_temperature']
del STU_labels['elsewhere_temperature']
STU_bounds['MOS1__alpha'] = main_Delta_XMM.likelihood.get_param('MOS1__alpha').bounds
STU_labels['MOS1__alpha'] = r"\alpha_{MOS}"

main_Delta_XMM.names = main_Delta_XMM.likelihood.names
main_Delta_XMM.names += names
main_Delta_XMM.bounds = STU_bounds
main_Delta_XMM.labels = STU_labels
main_Delta_XMM.runs = xpsi.Runs.load_runs(ID = r"$\Delta$, XMM",
                                           run_IDs = ['run1'],
                                           roots = [PATH_Delta_XMM],
                                           base_dirs = [PATH_dir],
                                           use_nestcheck=[True],
                                           kde_settings=getdist_kde_settings,
                                           likelihood = main_Delta_XMM.likelihood,
                                           names = main_Delta_XMM.names,
                                           bounds = main_Delta_XMM.bounds,
                                           labels = main_Delta_XMM.labels,
                                           implementation = 'multinest',
                                           overwrite_transformed = True)
print("\n Delta XMM metadata loaded successfully")