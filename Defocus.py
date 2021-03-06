# This file is part of holoaverage.
# Copyright (c) 2018 Tore Niermann
#
# holoaverage is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# holoaverage is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with holoaverage.  If not, see <http://www.gnu.org/licenses/>.
import numpy as np
from math import sqrt

from FFT import aligned_copy
from Series import DataSet
from Grid import Grid

__all__ = ('calcWaveLength', 'propagate')


def calcWaveLength(voltage_kV):
    """Returns electron wave length in [nm] for given voltage.

    Arguments:
        voltage_kV
            float: Voltage [kV]
    """
    # Constants from CODATA 2006
    PLANCK_CONSTANT_eVs = 4.13566733e-15    # eVs
    SPEED_OF_LIGHT_m_s = 299792458.0        # m/s
    ELECTRON_MASS_keV = 0.510998910e+3      # keV

    return PLANCK_CONSTANT_eVs * SPEED_OF_LIGHT_m_s * 1e+6 / sqrt(voltage_kV * (2.0 * ELECTRON_MASS_keV + voltage_kV))


def propagate(data, distance, voltage):
    """
    Propagates wave by distance.

    :param data: DataSet to propagate
    :param distance: Distance to propagate in nm
    :param voltage: Acceleration voltage in kV. If not given taken from data
    """
    grid = Grid.fromDataSet(data)
    waveLength = calcWaveLength(voltage)
    qy, qx = grid.getRcprGrid()

    tmp = grid.forwardFFT(aligned_copy(data.array))
    preD = np.pi * waveLength * distance
    if grid.hasDiagonalSampling():
        propX = np.exp(1.0j * preD * qx**2)
        propY = np.exp(1.0j * preD * qy**2)
        tmp *= propX * propY
    else:
        tmp *= np.exp(1.0j * preD * (qx**2 + qy**2))
    result = DataSet(data.shape, data.dtype)
    result.attrs.update(data.attrs)
    grid.backwardFFT(tmp, out=result.array)
    return result
