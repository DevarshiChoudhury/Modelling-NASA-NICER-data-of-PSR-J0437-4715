import numpy as np
import math

import xpsi

from xpsi import Parameter
from xpsi.utils import make_verbose

minCH = 30
maxCH = 300
NtotCH = 1501
NtotEI = 3451

class CustomInstrument(xpsi.Instrument):
    """ Methods and attributes specific to the NICER instrument.

    Currently tailored to the NICER light-curve SWG model specification.

    """
    def construct_matrix(self):
        """ Implement response matrix parameterisation. """
        # Multiplying beta to response matrix
        beta_d = self['beta'] * 0.01**2  # beta_d = beta * (d kpc)^2
        matrix = beta_d*self.matrix
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
    def from_SWG(cls,
                 bounds, values,
                 ARF, RSP,
                 max_input, min_input=0,
                 channel_edges = None):
        """ Constructor which converts files into :class:`numpy.ndarray`s.

        :param str ARF: Path to ARF which is compatible with
                                :func:`numpy.loadtxt`.

        :param str RSP: Path to RSP which is compatible with
                                :func:`numpy.loadtxt`.

        """
        ARF = np.loadtxt(ARF, dtype=np.double, skiprows=3)
        RSP = np.loadtxt(RSP, dtype=np.double, skiprows=3, usecols=-1)

        if channel_edges:
            channel_edges = np.loadtxt(channel_edges, dtype=np.double, skiprows=3)

        matrix = np.zeros((NtotCH,NtotEI))

        for i in range(NtotEI):
            matrix[:,i] = RSP[i*NtotCH:(i+1)*NtotCH]

        if min_input != 0:
            min_input = int(min_input)

        max_input = int(max_input)

        edges = np.zeros(ARF[min_input:max_input,4].shape[0]+1, dtype=np.double)

        edges[0] = ARF[min_input,1]; edges[1:] = ARF[min_input:max_input,2]

        RSP = matrix[minCH:maxCH,min_input:max_input]
        M = RSP
        for i in range(np.shape(M)[0]):
            Mi = M[i,:]
            if len(Mi[Mi>0.])==0:
                print ("All zeros: i ",i)
        for j in range(np.shape(M)[1]):
            Mj = M[:,j]
            if len(Mj[Mj>0.])==0:
                print ("All zeros: j ",j)
        

        channels = np.arange(minCH, maxCH)

        beta = Parameter('beta',
                          strict_bounds = (20.0,60.0),
                          bounds = bounds.get('beta', None),
                          doc='Units of kpc^-2',
                          symbol = r'$\beta$',
                          value = values.get('beta', None))

        return cls(RSP, edges, channels, channel_edges[minCH:maxCH+1,1], beta)
