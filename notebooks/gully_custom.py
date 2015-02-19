#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filename: gully_custom.py

import numpy as np
from astropy.io import fits


@np.vectorize
def specType(SpT):
    '''
  Converts between float and letter/number M, L, T and Y spectral types
  (e.g. 14.5 => 'L4.5' and 'T3' => 23).

  Returns float spectral type between -20.0 and 39.9
    '''
    if isinstance(SpT, str) and (SpT[0] in ['G', 'K', 'M', 'L', 'T', 'Y']) \
       and float(SpT[1:]) < 10:
        try:
            return [l + float(SpT[1:]) for m, l in
                    zip(['G', 'K', 'M', 'L', 'T', 'Y'],
                        [-20, -10, 0, 10, 20, 30]) if m == SpT[0]][0]
        except ValueError:
                print "Spectral type must be a float between 0 and 40 or a \
                      string of class G,K, M, L, T or Y."
                return np.NaN
    elif isinstance(SpT, float) or isinstance(SpT, int) and 0.0 <= SpT < 40.0:
        return '{}{}'.format('GKMLTY'[int(SpT // 10)], SpT % 10.)
    else:
        return np.NaN


@np.vectorize
def specTypePlus(SpT):
    '''
    Converts between float and letter/number M, L, T and Y spectral types
    (e.g. 14.5 => 'L4.5' and 'T3' => 23).  Additionally parses luminosity class
    "III" -> "Giant"

    Returns float spectral type between -20.0 and 39.9
    '''

    L_class = ''
    pec_class = ''
    SpT_e = 1.0
    if (SpT != SpT):
        SpT = 'junk'
    if (SpT.find('III') != -1):
        L_class = 'Giant'
        SpT = SpT.strip('III')
    if (SpT.find('V') != -1):
        L_class = 'Dwarf'
        SpT = SpT.strip('V')
    if (SpT.find('(IR)') != -1):
        pec_class = 'IR'
        SpT = SpT.strip('(IR)')
    if (SpT.count('β') == 1):
        L_class = 'β'
        SpT = SpT.strip('β')
    if (SpT.count('δ') == 1):
        L_class = 'δ'
        SpT = SpT.strip('δ')
    if (SpT.count('γ') == 1):
        L_class = 'γ'
        SpT = SpT.strip('γ')
    if (SpT.find('V') != -1):
        L_class = 'Dwarf'
        SpT = SpT.strip('V')
    if (SpT.find('pec') != -1):
        pec_class = 'pec'
        SpT = SpT.strip('pec')
    if (SpT.find('p') != -1):
        pec_class = 'pec'
        SpT = SpT.strip('p')
    if (SpT.find('blue') != -1):
        pec_class = 'blue'
        SpT = SpT.strip('blue')
    if (SpT.find('e') != -1):
        pec_class = 'e'
        SpT = SpT.strip('e')
    if (SpT.find('FGK') != -1):
        SpT_e = 10.0
        SpT = SpT.strip('GK')
    if (SpT.find('+') != -1):
        pec_class = 'dbl'
        SpT = ''
    if (SpT.find('::') != -1):
        SpT_e = 4.0
        SpT = SpT.strip('::')
    if (SpT.find(':') != -1):
        SpT_e = 2.0
        SpT = SpT.strip(':')

    if isinstance(SpT, str) and len(SpT) < 2:
        return np.NaN, SpT_e, L_class, pec_class
    if isinstance(SpT, str) and (SpT[0] in ['G', 'K', 'M', 'L', 'T', 'Y']) \
       and float(SpT[1:]) < 10:
        try:
            num_value = [l + float(SpT[1:]) for m, l in
                         zip(['G', 'K', 'M', 'L', 'T', 'Y'],
                         [-20, -10, 0, 10, 20, 30]) if m == SpT[0]][0]
            return num_value, SpT_e, L_class, pec_class
        except ValueError:
                print "Spectral type must be a float between 0 and 40 or a \
                      string of class G,K, M, L, T or Y."
                return np.NaN, np.NaN, L_class, pec_class
    elif isinstance(SpT, float) or isinstance(SpT, int) and 0.0 <= SpT < 40.0:
        return '{}{}'.format('GKMLTY'[int(SpT // 10)], SpT % 10.)
    else:
        return np.NaN, SpT_e, L_class, pec_class


def read_spec(filename):
    '''
    Read a spectrum from FITS assuming a (3 x N) data dimension, wl, fl, err.

    Parameters
    ----------
    filename : string
    name of the fits file with the data

    Returns
    -------
    wavelength : np.ndarray
    wavelength (in Ang)
    flux : np.ndarray
    flux (in erg/s/cm**2)
    error : np.ndarray
    error (in erg/s/cm**2)
    '''

    sp = fits.open(filename)
    header = sp[0].header
    wavelength = sp[0].data[0, :]
    flux = sp[0].data[1, :]
    error = sp[0].data[2, :]
    return wavelength, flux, error, header


def calculate_index_single_numerator(wavelength, flux, df):
    '''
    Calculate the index for a single numerator.
    Parameters
    ----------
    wavelength: np.ndarray
    array of wavelength values
    flux: np.ndarray
    array of flux values
    df: DataFrame
    a DataFrame with a single row with N-, and D- start and end.

    Returns
    -------
    index: np.float
    the value of the index
    unc: np.float
    the uncertainty in the value of the index
    '''
    Nids = (df.N_Start < wavelength) & (wavelength < df.N_End)
    Dids = (df.D_Start < wavelength) & (wavelength < df.D_End)

    if (Nids.sum == 0) or (Dids.sum == 0):
        return np.NaN, np.NaN

    index_out = np.mean(flux[Nids]) / np.mean(flux[Dids])
    sigma_D = np.std(flux[Dids] / len(Dids))
    index_unc = np.sqrt((index_out - np.mean(flux[Nids]) /
                        (np.mean(flux[Dids]) + sigma_D)) ** 2)

    return index_out, index_unc
