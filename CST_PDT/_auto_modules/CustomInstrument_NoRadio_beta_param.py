""" Instrument module for X-PSI ST-U + NSX-H modelling of NICER PSR J0437-4715 event data. """

from __future__ import print_function, division

import numpy as np
import math

import xpsi

# from xpsi.utils import make_verbose

from xpsi import Parameter#, make_verbose

class CustomInstrument(xpsi.Instrument):
    """ XTI, and XTI. """
    def construct_matrix(self):
        """ Implement response matrix parameterisation. """

        # Multiplying beta to response matrix
        beta_d = self['energy_independent_effective_area_scaling_factor'] * 0.01**2  # beta_d = beta * (d kpc)^2
                                                                                     # d fixed to 0.01kpc in config file
        matrix = beta_d * self.matrix
        matrix[matrix < 0.0] = 0.0

        return matrix

    def __call__(self, signal, *args):
        """ Overwrite. """

        matrix = self.construct_matrix()

        self._cached_signal = np.dot(matrix, signal)

        return self._cached_signal

    @classmethod
   #  @make_verbose('Loading XTI response matrix',
   #                'Response matrix loaded')
    def XTI(cls,
              bounds,
              values,
              ARF,
              RMF,
              channel_energies,
              max_input,
              max_channel,
              min_input=0,
              min_channel=0,
              effective_area_scaling_factor=1.0,
              ARF_skiprows=0,
              ARF_low_column=1,
              ARF_high_column=2,
              ARF_area_column=3,
              RMF_skiprows=0,
              RMF_usecol=-1,
              channel_energies_skiprows=0,
              channel_energies_low_column=0,
              **kwargs):
        """ Load XTI instrument response matrix. """

        beta = Parameter('energy_independent_effective_area_scaling_factor',
                          strict_bounds = (0.0,100.0),
                          bounds = bounds.get('energy_independent_effective_area_scaling_factor', None),
                          doc='XTI energy-independent effective area scaling factor (Units of kpc^-2)',
                          symbol = r'$\beta_{\rm XTI}$',
                          value = values.get('energy_independent_effective_area_scaling_factor', None))

        # check the loading assumptions and comment out the exception throw if they are true
        #raise NotImplementedError('Implement the class method for loading the XTI instrument.')

        ARF = np.loadtxt(ARF, dtype=np.double, skiprows=ARF_skiprows)
        RMF = np.loadtxt(RMF, dtype=np.double, skiprows=RMF_skiprows, usecols=RMF_usecol)
        channel_energies = np.loadtxt(channel_energies, dtype=np.double, skiprows=channel_energies_skiprows)

        matrix = np.zeros((channel_energies.shape[0], ARF.shape[0]))

        for i in range(ARF.shape[0]):
           matrix[:,i] = RMF[i*channel_energies.shape[0]:(i+1)*channel_energies.shape[0]]

        max_input = int(max_input)
        if min_input != 0:
           min_input = int(min_input)

        edges = np.zeros(max_input - min_input + 1, dtype=np.double)

        edges[0] = ARF[min_input, ARF_low_column]; edges[1:] = ARF[min_input:max_input, ARF_high_column]

        RSP = np.zeros((max_channel - min_channel,
                       max_input - min_input), dtype=np.double)

        for i in range(RSP.shape[0]):
           RSP[i,:] = matrix[i+min_channel, min_input:max_input] * ARF[min_input:max_input, ARF_area_column] * effective_area_scaling_factor

        channels = np.arange(min_channel, max_channel)

        return cls(RSP,
                   edges,
                   channels,
                   channel_energies[min_channel:max_channel+1,channel_energies_low_column],
                   beta, **kwargs)
