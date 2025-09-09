import numpy as np
import math

import xpsi

from xpsi import Parameter
from xpsi.utils import make_verbose

nCH_NICER = 1501
nIN_NICER = 3451
nIN_XMM = 2400
nCH_XMM = 2400

class CustomInstrument(xpsi.Instrument):
    """ Methods and attributes specific to the NICER instrument.

    Currently tailored to the NICER light-curve SWG model specification.

    """
    def construct_matrix(self):
        """ Implement response matrix parameterisation. """        

        matrix = self['alpha'] * self.matrix
        matrix[matrix < 0.0] = 0.0
        return matrix

    def __call__(self, signal, *args):
        """ Overwrite. """

        matrix = self.construct_matrix()

        self._cached_signal = np.dot(matrix, signal)

        return self._cached_signal

    @classmethod
    @make_verbose('Loading response matrix',
                  'Response matrix loaded')
    def NICER_XTI(cls,
                  bounds, values,
                  ARF, RSP,
                  max_input,
                  max_channel,
                  min_input=0,
                  min_channel=0,
                  channel_edges=None,
                  **kwargs):
        """ Load NICER XTI instrument response matrix. """
        ARF = np.loadtxt(ARF, dtype=np.double, skiprows=3)
        RSP = np.loadtxt(RSP, dtype=np.double, skiprows=3, usecols=-1)

        if channel_edges:
            channel_edges = np.loadtxt(channel_edges, dtype=np.double, skiprows=3)

        matrix = np.zeros((nCH_NICER,nIN_NICER))

        for i in range(nIN_NICER):
            matrix[:,i] = RSP[i*nCH_NICER:(i+1)*nCH_NICER]

        if min_input != 0:
            min_input = int(min_input)

        max_input = int(max_input)

        edges = np.zeros(ARF[min_input:max_input,3].shape[0]+1, dtype=np.double)

        edges[0] = ARF[min_input,1]; edges[1:] = ARF[min_input:max_input,2]

        RSP = matrix[min_channel:max_channel,min_input:max_input]
        
        channels = np.arange(min_channel, max_channel)

        alpha = Parameter('alpha',
                          strict_bounds = (0.1,1.9),
                          bounds = bounds.get('alpha', None),
                          doc='NICER XTI energy-independent scaling factor',
                          symbol = r'$\alpha_{\rm XTI}$',
                          value = values.get('alpha', None))

        return cls(RSP, edges, channels, channel_edges[min_channel:max_channel+1,1],
                   alpha, **kwargs)


    @classmethod
    @make_verbose('Loading XMM-MOS1 response matrix',
                  'Response matrix loaded')
    def XMM_MOS1(cls,
                 bounds, values,
                 RSP,
                 max_input,
                 max_channel,
                 min_input=0,
                 min_channel=0,
                 channel_edges = None,
                 **kwargs):
        """ Load XMM-MOS1 instrument response matrix. """
        RSP = np.loadtxt(RSP, dtype=np.double, skiprows=3, usecols=-1)

        if channel_edges:
            channel_edges = np.loadtxt(channel_edges, dtype=np.double, skiprows=3)

        if min_input != 0:
            min_input = int(min_input)

        max_input = int(max_input)

        matrix = np.zeros((nCH_XMM, max_input))

        for i in range(matrix.shape[1]):
            matrix[:,i] = RSP[i*nCH_XMM:(i+1)*nCH_XMM]

        edges = np.zeros(max_input - min_input + 1, dtype=np.double)

        edges[0] = channel_edges[min_input,1]; edges[1:] = channel_edges[min_input:max_input,2]

        RSP = matrix[min_channel:max_channel,min_input:max_input]

        channels = np.arange(min_channel, max_channel)

        alpha = Parameter('alpha',
                          strict_bounds = (0.1,1.9),
                          bounds = bounds.get('alpha', None),
                          doc='XMM-MOS1 energy-independent scaling factor',
                          symbol = r'$\alpha_{\rm MOS1}$',
                          value = values.get('alpha', None))

        return cls(RSP, edges, channels,
                   channel_edges[min_channel:max_channel+1,1],
                   alpha, **kwargs)

    @classmethod
    @make_verbose('Loading XMM-MOS2 response matrix',
                  'Response matrix loaded')
    def XMM_MOS2(cls,
                 bounds, values,
                 RSP,
                 max_input,
                 max_channel,
                 min_input=0,
                 min_channel=0,
                 channel_edges = None,
                 **kwargs):
        """ Load XMM-MOS2 instrument response matrix. """
        RSP = np.loadtxt(RSP, dtype=np.double, skiprows=3, usecols=-1)

        if channel_edges:
            channel_edges = np.loadtxt(channel_edges, dtype=np.double, skiprows=3)

        if min_input != 0:
            min_input = int(min_input)

        max_input = int(max_input)

        matrix = np.zeros((nCH_XMM, max_input))

        for i in range(matrix.shape[1]):
            matrix[:,i] = RSP[i*nCH_XMM:(i+1)*nCH_XMM]

        edges = np.zeros(max_input - min_input + 1, dtype=np.double)

        edges[0] = channel_edges[min_input,1]; edges[1:] = channel_edges[min_input:max_input,2]

        RSP = matrix[min_channel:max_channel,min_input:max_input]

        channels = np.arange(min_channel, max_channel)

        alpha = Parameter('alpha',
                          strict_bounds = (0.1,1.9),
                          bounds = bounds.get('alpha', None),
                          doc='XMM-MOS2 energy-independent scaling factor',
                          symbol = r'$\alpha_{\rm MOS2}$',
                          value = values.get('alpha', None))

        return cls(RSP, edges, channels,
                   channel_edges[min_channel:max_channel+1,1],
                   alpha, **kwargs)
