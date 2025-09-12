""" Instrument module for X-PSI CST+PDT modelling of NICER & XMM PSR J0437-4715 event data. """

import numpy as np
import math

import xpsi

from xpsi import Parameter
from xpsi.utils import make_verbose

class CustomInstrument(xpsi.Instrument):
    """ XTI, MOS1, and MOS2. """
    def construct_matrix(self):
        """ Implement response matrix parameterisation. """
        matrix = self['energy_independent_effective_area_scaling_factor'] * self.matrix
        matrix[matrix < 0.0] = 0.0

        return matrix

    def __call__(self, signal, *args):
        """ Overwrite. """

        matrix = self.construct_matrix()

        self._cached_signal = np.dot(matrix, signal)

        return self._cached_signal

    @classmethod
    @make_verbose('Loading XTI response matrix',
                  'Response matrix loaded')
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

        alpha = Parameter('energy_independent_effective_area_scaling_factor',
                          strict_bounds = (0.1,1.9),
                          bounds = bounds.get('energy_independent_effective_area_scaling_factor', None),
                          doc='XTI energy-independent effective area scaling factor',
                          symbol = r'$\alpha_{\rm XTI}$',
                          value = values.get('energy_independent_effective_area_scaling_factor',
                                             1.0 if bounds.get('energy_independent_effective_area_scaling_factor', None) is None else None))

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
                   alpha, **kwargs)

    @classmethod
    @make_verbose('Loading MOS1 response matrix',
                  'Response matrix loaded')
    def MOS1(cls,
              bounds,
              values,
              RSP,
              channel_edges,
              max_input,
              max_channel,
              min_input=0,
              min_channel=0,
              effective_area_scaling_factor=1.0,
              RSP_skiprows=0,
              RSP_usecol=-1,
              channel_edges_skiprows=0,
              **kwargs):
        """ Load MOS1 instrument response matrix. """

        alpha = Parameter('energy_independent_effective_area_scaling_factor',
                          strict_bounds = (0.1,1.9),
                          bounds = bounds.get('energy_independent_effective_area_scaling_factor', None),
                          doc='MOS1 energy-independent effective area scaling factor',
                          symbol = r'$\alpha_{\rm MOS1}$',
                          value = values.get('energy_independent_effective_area_scaling_factor',
                                             1.0 if bounds.get('energy_independent_effective_area_scaling_factor', None) is None else None))

        # check the loading assumptions and comment out the exception throw if they are true
        # raise NotImplementedError('Implement the class method for loading the MOS1 instrument.')

        RSP = np.loadtxt(RSP, dtype=np.double, skiprows=RSP_skiprows, usecols=RSP_usecol)
        channel_edges = np.loadtxt(channel_edges, dtype=np.double, skiprows=channel_edges_skiprows)

        max_input = int(max_input)
        if min_input != 0:
           min_input = int(min_input)

        matrix = np.zeros((channel_edges.shape[0], max_input))

        for i in range(matrix.shape[1]):
           matrix[:,i] = RSP[i*channel_edges.shape[0]:(i+1)*channel_edges.shape[0]]

        edges = np.zeros(max_input - min_input + 1, dtype=np.double)

        edges[0] = channel_edges[min_input, 1]; edges[1:] = channel_edges[min_input:max_input, 2]

        RSP = matrix[min_channel:max_channel,min_input:max_input]

        channels = np.arange(min_channel, max_channel)

        return cls(RSP,
                   edges,
                   channels,
                   channel_edges[min_channel:max_channel+1,1],
                   alpha, **kwargs)

    @classmethod
    @make_verbose('Loading MOS2 response matrix',
                  'Response matrix loaded')
    def MOS2(cls,
              bounds,
              values,
              RSP,
              channel_edges,
              max_input,
              max_channel,
              min_input=0,
              min_channel=0,
              effective_area_scaling_factor=1.0,
              RSP_skiprows=0,
              RSP_usecol=-1,
              channel_edges_skiprows=0,
              **kwargs):
        """ Load MOS2 instrument response matrix. """

        alpha = Parameter('energy_independent_effective_area_scaling_factor',
                          strict_bounds = (0.1,1.9),
                          bounds = bounds.get('energy_independent_effective_area_scaling_factor', None),
                          doc='MOS2 energy-independent effective area scaling factor',
                          symbol = r'$\alpha_{\rm MOS2}$',
                          value = values.get('energy_independent_effective_area_scaling_factor',
                                             1.0 if bounds.get('energy_independent_effective_area_scaling_factor', None) is None else None))

        # check the loading assumptions and comment out the exception throw if they are true
        #raise NotImplementedError('Implement the class method for loading the MOS2 instrument.')

        RSP = np.loadtxt(RSP, dtype=np.double, skiprows=RSP_skiprows, usecols=RSP_usecol)
        channel_edges = np.loadtxt(channel_edges, dtype=np.double, skiprows=channel_edges_skiprows)

        max_input = int(max_input)
        if min_input != 0:
           min_input = int(min_input)

        matrix = np.zeros((channel_edges.shape[0], max_input))

        for i in range(matrix.shape[1]):
           matrix[:,i] = RSP[i*channel_edges.shape[0]:(i+1)*channel_edges.shape[0]]

        edges = np.zeros(max_input - min_input + 1, dtype=np.double)

        edges[0] = channel_edges[min_input, 1]; edges[1:] = channel_edges[min_input:max_input, 2]

        RSP = matrix[min_channel:max_channel,min_input:max_input]

        channels = np.arange(min_channel, max_channel)

        return cls(RSP,
                   edges,
                   channels,
                   channel_edges[min_channel:max_channel+1,1],
                   alpha, **kwargs)
