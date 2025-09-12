################################################################################
## N.B.                                                                       ##
##                                                                            ##
## Here we load all modules and outputs corresponding to the CST+PDT model.   ##
## However, this can lead to very long loading times when imporing this file  ##
## in a postprocessing notebook (~ 10 hours depending on the system).         ##
## The main contributors to the long loading times are modules and outputs    ##
## corresponding to high MultiNest resolution runs (filenames and variables   ##
## involving the nomenclature `hiMN_lowXPSI`).                                ##
##                                                                            ##
## However, the modularity of this script and the postprocessing notebook     ##
## implies that you can selectively load the runs you wish to probe, by       ##
## commenting out the unwanted modules.                                       ##
##                                                                            ##
################################################################################

## Initialization

import xpsi
from xpsi.global_imports import gravradius

## Uncomment the following and adjust paths if needed
# import sys
# sys.path.append('./_auto_modules')

PATH_dir = './CST_PDT_outputs/'
PATH_3C50_NoBKG = '3C50_NoBKG/init/nlive4000_expf3.3_noCONST_noMM_tol0.1'
PATH_3C50_NoBKG_hiMN_lowXPSI = '3C50_NoBKG_hiMN_lowXPSI_res/res1/nlive20000_expf3.3_noCONST_noMM_tol0.1'
PATH_3C50_BKG_lower_smooth_3sigma = '3C50_BKG_lower_smooth_3sigma/init/nlive4000_expf3.3_noCONST_noMM_tol0.1'
PATH_3C50_BKG_AGN_smooth_3sigma = '3C50_BKG_AGN_smooth_3sigma/init/nlive4000_expf3.3_noCONST_noMM_tol0.1'
PATH_3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI = '3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI_res/resume7/nlive20000_expf3.3_noCONST_noMM_tol0.1'
PATH_3C50_BKG_AGN_smooth_3sigma_and_XMM_noMM = '3C50_BKG_AGN_smooth_3sigma_and_XMM/noMM/res1/nlive4000_expf3.3_noCONST_noMM_tol0.1'
PATH_3C50_BKG_AGN_smooth_3sigma_and_XMM_MM = '3C50_BKG_AGN_smooth_3sigma_and_XMM/MM/res1/nlive4000_expf3.3_noCONST_MM_tol0.1'
PATH_3C50_NoBKG_NoRadio = '3c50_NoBKG_NoRadio/res2/nlive4000_expf3.3_noCONST_noMM_tol0.1'


## Importing all relevant modules for post-processing

from _auto_modules import main_3C50_NoBKG
print("3C50_NoBKG (high MN low X-PSI res) log-likelihood value (using high X-PSI res): ",
      main_3C50_NoBKG.likelihood(main_3C50_NoBKG.p_maxL_hiMN_lowXPSI_res, reinitialise=True))
print("\n3C50_NoBKG (low MN high X-PSI res) log-likelihood value (using high X-PSI res): ",
      main_3C50_NoBKG.likelihood(main_3C50_NoBKG.p_maxL_lowMN_hiXPSI_res, reinitialise=True))
print("3C50_NoBKG module loaded \n")

from _auto_modules import main_3C50_NoBKG_hiMN_lowXPSI
print("\n3C50_NoBKG_hiMN_lowXPSI (low MN high X-PSI res) log-likelihood value (using low X-PSI res): ",
      main_3C50_NoBKG_hiMN_lowXPSI.likelihood(main_3C50_NoBKG_hiMN_lowXPSI.p_maxL_lowMN_hiXPSI_res, reinitialise=True))
print("3C50_NoBKG_hiMN_lowXPSI (high MN low X-PSI res) log-likelihood value (using low X-PSI res): ",
      main_3C50_NoBKG_hiMN_lowXPSI.likelihood(main_3C50_NoBKG_hiMN_lowXPSI.p_maxL_hiMN_lowXPSI_res, reinitialise=True))
print("3C50_NoBKG_hiMN_lowXPSI module loaded \n")

from _auto_modules import main_3C50_BKG_lower_smooth_3sigma
print("\n3C50_BKG_lower_smooth_3sigma log-likelihood value: ",
      main_3C50_BKG_lower_smooth_3sigma.likelihood(main_3C50_BKG_lower_smooth_3sigma.p_maxL, reinitialise=True))
print("3C50_BKG_lower_smooth_3sigma module loaded \n")

from _auto_modules import main_3C50_BKG_AGN_smooth_3sigma
print("3C50_BKG_AGN_smooth_3sigma (high MN low X-PSI res) log-likelihood value (using high X-PSI res): ",
      main_3C50_BKG_AGN_smooth_3sigma.likelihood(main_3C50_BKG_AGN_smooth_3sigma.p_maxL_hiMN_lowXPSI_res, reinitialise=True))
print("\n3C50_BKG_AGN_smooth_3sigma (low MN high X-PSI res) log-likelihood value (using high X-PSI res): ",
      main_3C50_BKG_AGN_smooth_3sigma.likelihood(main_3C50_BKG_AGN_smooth_3sigma.p_maxL_lowMN_hiXPSI_res, reinitialise=True))
print("3C50_BKG_AGN_smooth_3sigma module loaded \n")

from _auto_modules import main_3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI
print("\n3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI (low MN high X-PSI res) log-likelihood value (using low X-PSI res): ",
      main_3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI.likelihood(main_3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI.p_maxL_lowMN_hiXPSI_res, reinitialise=True))
print("3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI (high MN low X-PSI res) log-likelihood value (using low X-PSI res): ",
      main_3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI.likelihood(main_3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI.p_maxL_hiMN_lowXPSI_res, reinitialise=True))
print("3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI module loaded \n")

from _auto_modules_NICERxXMM import main as main_XMM
print("\n3C50_BKG_AGN_smooth_3sigma_and_XMM noMM log-likelihood value: ",
      main_XMM.likelihood(main_XMM.p_maxL_noMM, reinitialise=True))
print("3C50_BKG_AGN_smooth_3sigma_and_XMM MM log-likelihood value: ",
      main_XMM.likelihood(main_XMM.p_maxL_MM, reinitialise=True))
print("3C50_BKG_AGN_smooth_3sigma_and_XMM module loaded \n")

from _auto_modules import main_3C50_NoBKG_NoRadio_beta_param
print("\n3C50_NoBKG_NoRadio_beta_param (low MN high X-PSI res) log-likelihood value (using high X-PSI res): ",
      main_3C50_NoBKG_NoRadio_beta_param.likelihood(main_3C50_NoBKG_NoRadio_beta_param.p_maxL, reinitialise=True))
print("3C50_NoBKG_NoRadio module loaded \n")

print('All modules loaded \n\nLoading metadata...')


## Specifying parameter bounds and labels for post-processing

# Additional parameters specified in CustomPrior.transfom to be plotted
names = ['compactness', 'N_H_e20_units', 'p__phase_shift_shifted', 's__phase_shift_shifted']

CST_PDT_bounds = {'mass': main_3C50_NoBKG.likelihood.get_param('mass').bounds,
                 'radius': main_3C50_NoBKG.likelihood.get_param('radius').bounds,
                 'cos_inclination': main_3C50_NoBKG.likelihood.get_param('cos_inclination').bounds,
                 'p__super_colatitude': main_3C50_NoBKG.likelihood.get_param('p__super_colatitude').bounds,
                 'p__super_radius': main_3C50_NoBKG.likelihood.get_param('p__super_radius').bounds,
                 'p__super_temperature': main_3C50_NoBKG.likelihood.get_param('p__super_temperature').bounds,
                 'p__omit_radius': main_3C50_NoBKG.likelihood.get_param('p__omit_radius').bounds,
                 's__super_colatitude': main_3C50_NoBKG.likelihood.get_param('s__super_colatitude').bounds,
                 's__super_radius': main_3C50_NoBKG.likelihood.get_param('s__super_radius').bounds,
                 's__super_temperature': main_3C50_NoBKG.likelihood.get_param('s__super_temperature').bounds,
                 's__cede_colatitude': main_3C50_NoBKG.likelihood.get_param('s__cede_colatitude').bounds,
                 's__cede_radius': main_3C50_NoBKG.likelihood.get_param('s__cede_radius').bounds,
                 's__cede_azimuth': main_3C50_NoBKG.likelihood.get_param('s__cede_azimuth').bounds,
                 's__cede_temperature': main_3C50_NoBKG.likelihood.get_param('s__cede_temperature').bounds,
                 'XTI__energy_independent_effective_area_scaling_factor': main_3C50_NoBKG.likelihood.get_param('XTI__energy_independent_effective_area_scaling_factor').bounds,
                 'neutral_hydrogen_column_density': main_3C50_NoBKG.likelihood.get_param('neutral_hydrogen_column_density').bounds,
                 'p__phase_shift': main_3C50_NoBKG.likelihood.get_param('p__phase_shift').bounds,
                 's__phase_shift': main_3C50_NoBKG.likelihood.get_param('s__phase_shift').bounds,
                 'compactness': (gravradius(1.0)/16.0, 1.0/3.0),
                 'N_H_e20_units': (0.04, 2.0),
                 'p__phase_shift_shifted': (0.0, 1.0),
                 's__phase_shift_shifted': (0.0, 1.0)}

#TeX compatible labels
CST_PDT_labels = {'mass': r"M\;\mathrm{[M}_{\odot}\mathrm{]}",
                 'radius': r"R_{\mathrm{eq}}\;\mathrm{[km]}",
                 'cos_inclination': r"\mathrm{cos}(i)",
                 'p__super_colatitude': r"\Theta_{p}\;\mathrm{[rad]}",
                 'p__super_radius': r"\zeta_{p}\;\mathrm{[rad]}",
                 'p__super_temperature': r"\mathrm{log10}(T_{p}\;[\mathrm{K}])",
                 'p__omit_radius': r"\zeta_{o,p}\;\mathrm{[rad]}",
                 's__super_colatitude': r"\Theta_{s}\;\mathrm{[rad]}",
                 's__super_radius': r"\zeta_{s}\;\mathrm{[rad]}",
                 's__super_temperature': r"\mathrm{log10}(T_{s}\;[\mathrm{K}])",
                 's__cede_colatitude': r"\Theta_{c,s}\;\mathrm{[rad]}",
                 's__cede_radius': r"\zeta_{c,s}\;\mathrm{[rad]}",
                 's__cede_azimuth': r"\chi_{s}\;\mathrm{[rad]}",
                 's__cede_temperature': r"\mathrm{log10}(T_{c,s}\;[\mathrm{K}])",
                 'XTI__energy_independent_effective_area_scaling_factor': r"\alpha_{XTI}",
                 'neutral_hydrogen_column_density': r"N_{\mathrm{H}}\;\mathrm{[0.4 \times 10^{20}\;cm^{-2}]}",
                 'p__phase_shift': r"\phi_{p}\;\mathrm{[cycles]}",
                 's__phase_shift': r"\phi_{s}\;\mathrm{[cycles]}",
                 'compactness': r"M/R_{\mathrm{eq}}",
                 'N_H_e20_units': r"N_{\mathrm{H}}\;\mathrm{[10^{20}\;cm^{-2}]}",
                 'p__phase_shift_shifted': r"\phi_{p}\;\mathrm{[cycles]}",
                 's__phase_shift_shifted': r"\phi_{s}\;\mathrm{[cycles]}",}


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
main_3C50_NoBKG.bounds = CST_PDT_bounds
main_3C50_NoBKG.labels = CST_PDT_labels
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
print("\n 3C50_NoBKG metadata loaded successfully")

main_3C50_NoBKG_hiMN_lowXPSI.names = main_3C50_NoBKG_hiMN_lowXPSI.likelihood.names
main_3C50_NoBKG_hiMN_lowXPSI.names += names
main_3C50_NoBKG_hiMN_lowXPSI.bounds = CST_PDT_bounds
main_3C50_NoBKG_hiMN_lowXPSI.labels = CST_PDT_labels
main_3C50_NoBKG_hiMN_lowXPSI.runs = xpsi.Runs.load_runs(ID = r"3C50, No BKG, 20k-LP",
                                           run_IDs = ['run1'],
                                           roots = [PATH_3C50_NoBKG_hiMN_lowXPSI],
                                           base_dirs = [PATH_dir],
                                           use_nestcheck=[True],
                                           kde_settings=getdist_kde_settings,
                                           likelihood = main_3C50_NoBKG_hiMN_lowXPSI.likelihood,
                                           names = main_3C50_NoBKG_hiMN_lowXPSI.names,
                                           bounds = main_3C50_NoBKG_hiMN_lowXPSI.bounds,
                                           labels = main_3C50_NoBKG_hiMN_lowXPSI.labels,
                                           implementation = 'multinest',
                                           overwrite_transformed = True)
print("\n 3C50_NoBKG_hiMN_lowXPSI metadata loaded successfully")

main_3C50_BKG_lower_smooth_3sigma.names = main_3C50_BKG_lower_smooth_3sigma.likelihood.names
main_3C50_BKG_lower_smooth_3sigma.names += names
main_3C50_BKG_lower_smooth_3sigma.bounds = CST_PDT_bounds
main_3C50_BKG_lower_smooth_3sigma.labels = CST_PDT_labels
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

main_3C50_BKG_AGN_smooth_3sigma.names = main_3C50_BKG_AGN_smooth_3sigma.likelihood.names
main_3C50_BKG_AGN_smooth_3sigma.names += names
main_3C50_BKG_AGN_smooth_3sigma.bounds = CST_PDT_bounds
main_3C50_BKG_AGN_smooth_3sigma.labels = CST_PDT_labels
main_3C50_BKG_AGN_smooth_3sigma.runs = xpsi.Runs.load_runs(ID = r"3C50, BKG-33, AGN (smooth)",
                                                           run_IDs = ['run1'],
                                                           roots = [PATH_3C50_BKG_AGN_smooth_3sigma],
                                                           base_dirs = [PATH_dir],
                                                           use_nestcheck=[True],
                                                           kde_settings=getdist_kde_settings,
                                                           likelihood = main_3C50_BKG_AGN_smooth_3sigma.likelihood,
                                                           names = main_3C50_BKG_AGN_smooth_3sigma.names,
                                                           bounds = main_3C50_BKG_AGN_smooth_3sigma.bounds,
                                                           labels = main_3C50_BKG_AGN_smooth_3sigma.labels,
                                                           implementation = 'multinest',
                                                           overwrite_transformed = True)
print("\n 3C50_BKG_AGN_smooth_3sigma metadata loaded successfully")

main_3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI.names = main_3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI.likelihood.names
main_3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI.names += names
main_3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI.bounds = CST_PDT_bounds
main_3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI.labels = CST_PDT_labels
main_3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI.runs = xpsi.Runs.load_runs(ID = r"3C50, BKG-33, AGN (smooth), 20k-LP",
                                                           run_IDs = ['run1'],
                                                           roots = [PATH_3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI],
                                                           base_dirs = [PATH_dir],
                                                           use_nestcheck=[True],
                                                           kde_settings=getdist_kde_settings,
                                                           likelihood = main_3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI.likelihood,
                                                           names = main_3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI.names,
                                                           bounds = main_3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI.bounds,
                                                           labels = main_3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI.labels,
                                                           implementation = 'multinest',
                                                           overwrite_transformed = True)
print("\n 3C50_BKG_AGN_smooth_3sigma_hiMN_lowXPSI metadata loaded successfully")

CST_PDT_bounds['mass'] = main_3C50_NoBKG_NoRadio_beta_param.likelihood.get_param('mass').bounds
CST_PDT_bounds['cos_inclination'] = main_3C50_NoBKG_NoRadio_beta_param.likelihood.get_param('cos_inclination').bounds
CST_PDT_bounds['XTI__energy_independent_effective_area_scaling_factor'] = main_3C50_NoBKG_NoRadio_beta_param.likelihood.get_param('XTI__energy_independent_effective_area_scaling_factor').bounds
CST_PDT_labels['XTI__energy_independent_effective_area_scaling_factor'] = r"\beta\;\mathrm{[kpc^{-2}]}"

main_3C50_NoBKG_NoRadio_beta_param.names = main_3C50_NoBKG_NoRadio_beta_param.likelihood.names
main_3C50_NoBKG_NoRadio_beta_param.names += names
main_3C50_NoBKG_NoRadio_beta_param.bounds = CST_PDT_bounds
main_3C50_NoBKG_NoRadio_beta_param.labels = CST_PDT_labels
main_3C50_NoBKG_NoRadio_beta_param.runs = xpsi.Runs.load_runs(ID = r"3C50, No BKG, No Radio",
                                           run_IDs = ['run1'],
                                           roots = [PATH_3C50_NoBKG_NoRadio],
                                           base_dirs = [PATH_dir],
                                           use_nestcheck=[True],
                                           kde_settings=getdist_kde_settings,
                                           likelihood = main_3C50_NoBKG_NoRadio_beta_param.likelihood,
                                           names = main_3C50_NoBKG_NoRadio_beta_param.names,
                                           bounds = main_3C50_NoBKG_NoRadio_beta_param.bounds,
                                           labels = main_3C50_NoBKG_NoRadio_beta_param.labels,
                                           implementation = 'multinest',
                                           overwrite_transformed = True)
print("\n 3C50_NoBKG_NoRadio metadata loaded successfully")

CST_PDT_bounds['mass'] = main_XMM.likelihood.get_param('mass').bounds
CST_PDT_bounds['cos_inclination'] = main_XMM.likelihood.get_param('cos_inclination').bounds
CST_PDT_bounds['XTI__energy_independent_effective_area_scaling_factor'] = main_XMM.likelihood.get_param('XTI__energy_independent_effective_area_scaling_factor').bounds
CST_PDT_labels['XTI__energy_independent_effective_area_scaling_factor'] = r"\alpha_{XTI}"
CST_PDT_bounds['MOS1__energy_independent_effective_area_scaling_factor'] = main_XMM.likelihood.get_param('MOS1__energy_independent_effective_area_scaling_factor').bounds
CST_PDT_labels['MOS1__energy_independent_effective_area_scaling_factor'] = r"\alpha_{MOS}"

main_XMM.names = main_XMM.likelihood.names
main_XMM.names += names
main_XMM.bounds = CST_PDT_bounds
main_XMM.labels = CST_PDT_labels
main_XMM.runs = xpsi.Runs.load_runs(ID = r"3C50, BKG-33, AGN (smooth), XMM",
                                run_IDs = ['No MM', 'MM'],
                                roots = [PATH_3C50_BKG_AGN_smooth_3sigma_and_XMM_noMM, PATH_3C50_BKG_AGN_smooth_3sigma_and_XMM_MM],
                                base_dirs = [PATH_dir] * 2,
                                use_nestcheck=[True, False],
                                kde_settings=getdist_kde_settings,
                                likelihood = main_XMM.likelihood,
                                names = main_XMM.names,
                                bounds = main_XMM.bounds,
                                labels = main_XMM.labels,
                                implementation = 'multinest',
                                overwrite_transformed = True)
print("\n 3C50_BKG_AGN_smooth_3sigma_and_XMM metadata loaded successfully")


print('\n \nAll metadata loaded successfully')