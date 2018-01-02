#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pytest
from goodvibes import GoodVibes as GV
from conftest import datapath


@pytest.mark.parametrize("path, QH, temp, E, ZPE, H, TS, TqhS, G, qhG", [
    # Grimme, 298.15K
    ('Al_298K.out', 'grimme', 298.15, -242.328708, 0.000000, -242.326347, 0.017670, 0.017670, -242.344018, -242.344018),
    ('Al_400K.out', 'grimme', 298.15, -242.328708, 0.000000, -242.326347, 0.017670, 0.017670, -242.344018, -242.344018),
    ('allene.out', 'grimme', 298.15, -116.569605, 0.053913, -116.510916, 0.027618, 0.027621, -116.538534, -116.538537),
    ('CuCN.out', 'grimme', 298.15, -289.005463, 0.006594, -288.994307, 0.025953, 0.025956, -289.020260, -289.020264),
    ('ethane.out', 'grimme', 298.15, -79.830421, 0.075238, -79.750770, 0.027523, 0.027525, -79.778293, -79.778295),
    ('ethane_spc.out', 'grimme', 298.15, -79.830421, 0.075238, -79.750770, 0.027523, 0.027525, -79.778293, -79.778295),
    ('ethane_TZ.out', 'grimme', 298.15, -79.858399, None, None, None, None, None, None),
    ('H2O.out', 'grimme', 298.15, -76.368128, 0.020772, -76.343577, 0.021458, 0.021458, -76.365035, -76.365035),
    ('HCN_singlet.out', 'grimme', 298.15, -93.358851, 0.015978, -93.339373, 0.022896, 0.022896, -93.362269, -93.362269),
    ('HCN_triplet.out', 'grimme', 298.15, -93.153787, 0.012567, -93.137780, 0.024070, 0.024070, -93.161850, -93.161850),
    ('methylaniline.out', 'grimme', 298.15, -326.664901, 0.142118, -326.514489, 0.039668, 0.039535, -326.554157, -326.554024),
    # Grimme, 100.0K
    ('Al_298K.out', 'grimme', 100.0, -242.328708, 0.000000, -242.327916, 0.005062, 0.005062, -242.332978, -242.332978),
    ('Al_400K.out', 'grimme', 100.0, -242.328708, 0.000000, -242.327916, 0.005062, 0.005062, -242.332978, -242.332978),
    ('allene.out', 'grimme', 100.0, -116.569605, 0.053913, -116.514408, 0.007423, 0.007423, -116.521831, -116.521831),
    ('CuCN.out', 'grimme', 100.0, -289.005463, 0.006594, -288.997568, 0.006944, 0.006946, -289.004512, -289.004514),
    ('ethane.out', 'grimme', 100.0, -79.830421, 0.075238, -79.753900, 0.007558, 0.007559, -79.761458, -79.761459),
    ('ethane_spc.out', 'grimme', 100.0, -79.830421, 0.075238, -79.753900, 0.007558, 0.007559, -79.761458, -79.761459),
    ('ethane_TZ.out', 'grimme', 100.0, -79.858399, None, None, None, None, None, None),
    ('H2O.out', 'grimme', 100.0, -76.368128, 0.020772, -76.346089, 0.005812, 0.005812, -76.351901, -76.351901),
    ('HCN_singlet.out', 'grimme', 100.0, -93.358851, 0.015978, -93.341765, 0.006385, 0.006385, -93.348150, -93.348150),
    ('HCN_triplet.out', 'grimme', 100.0, -93.153787, 0.012567, -93.140111, 0.006803, 0.006803, -93.146915, -93.146915),
    ('methylaniline.out', 'grimme', 100.0, -326.664901, 0.142118, -326.521226, 0.009864, 0.009905, -326.531090, -326.531131),
    # Truhlar, 298.15K
    ('Al_298K.out', 'truhlar', 298.15, -242.328708, 0.000000, -242.326347, 0.017670, 0.017670, -242.344018, -242.344018),
    ('Al_400K.out', 'truhlar', 298.15, -242.328708, 0.000000, -242.326347, 0.017670, 0.017670, -242.344018, -242.344018),
    ('allene.out', 'truhlar', 298.15, -116.569605, 0.053913, -116.510916, 0.027618, 0.027618, -116.538534, -116.538534),
    ('CuCN.out', 'truhlar', 298.15, -289.005463, 0.006594, -288.994307, 0.025953, 0.025953, -289.020260, -289.020260),
    ('ethane.out', 'truhlar', 298.15, -79.830421, 0.075238, -79.750770, 0.027523, 0.027523, -79.778293, -79.778293),
    ('ethane_spc.out', 'truhlar', 298.15, -79.830421, 0.075238, -79.750770, 0.027523, 0.027523, -79.778293, -79.778293),
    ('ethane_TZ.out',  'truhlar', 298.15, -79.858399, None, None, None, None, None, None),
    ('H2O.out', 'truhlar', 298.15, -76.368128, 0.020772, -76.343577, 0.021458, 0.021458, -76.365035, -76.365035),
    ('HCN_singlet.out', 'truhlar', 298.15, -93.358851, 0.015978, -93.339373, 0.022896, 0.022896, -93.362269, -93.362269),
    ('HCN_triplet.out', 'truhlar', 298.15, -93.153787, 0.012567, -93.137780, 0.024070, 0.024070, -93.161850, -93.161850),
    ('methylaniline.out', 'truhlar', 298.15, -326.664901, 0.142118, -326.514489, 0.039668, 0.039668, -326.554157, -326.554157),
    # Truhlar, 100.0K
    ('Al_298K.out', 'truhlar', 100.0, -242.328708, 0.000000, -242.327916, 0.005062, 0.005062, -242.332978, -242.332978),
    ('Al_400K.out', 'truhlar', 100.0, -242.328708, 0.000000, -242.327916, 0.005062, 0.005062, -242.332978, -242.332978),
    ('allene.out', 'truhlar', 100.0, -116.569605, 0.053913, -116.514408, 0.007423, 0.007423, -116.521831, -116.521831),
    ('CuCN.out', 'truhlar', 100.0, -289.005463, 0.006594, -288.997568, 0.006944, 0.006944, -289.004512, -289.004512),
    ('ethane.out', 'truhlar', 100.0, -79.830421, 0.075238, -79.753900, 0.007558, 0.007558, -79.761458, -79.761458),
    ('ethane_spc.out', 'truhlar', 100.0, -79.830421, 0.075238, -79.753900, 0.007558, 0.007558, -79.761458, -79.761458),
    ('ethane_TZ.out', 'truhlar', 100.0, -79.858399, None, None, None, None, None, None),
    ('H2O.out', 'truhlar', 100.0, -76.368128, 0.020772, -76.346089, 0.005812, 0.005812, -76.351901, -76.351901),
    ('HCN_singlet.out', 'truhlar', 100.0, -93.358851, 0.015978, -93.341765, 0.006385, 0.006385, -93.348150, -93.348150),
    ('HCN_triplet.out', 'truhlar', 100.0, -93.153787, 0.012567, -93.140111, 0.006803, 0.006803, -93.146915, -93.146915),
    ('methylaniline.out', 'truhlar', 100.0, -326.664901, 0.142118, -326.521226, 0.009864, 0.009864, -326.531090, -326.531090),
])
def test_all(path, QH, temp, E, ZPE, H, TS, TqhS, G, qhG):
    # Defaults, no temp interval, no conc interval
    path = datapath(path)
    conc = GV.ATMOS / (GV.GAS_CONSTANT * temp)
    freq_cutoff, freq_scale_factor, solv, spc = 100.0, 1.0, 'none', False
    bbe = GV.calc_bbe(path, QH, freq_cutoff, temp, conc, freq_scale_factor, solv, spc)
    precision = 6 # if temp == 298.15 else 4e-4
    assert E == round(bbe.scf_energy, precision)
    if hasattr(bbe, "gibbs_free_energy"):
        assert ZPE == round(bbe.zpe, precision)
        assert H == round(bbe.enthalpy, precision)
        assert TS == round(temp * bbe.entropy, precision)
        assert TqhS == round(temp * bbe.qh_entropy, precision)
        assert G == round(bbe.gibbs_free_energy, precision)
        assert qhG == round(bbe.qh_gibbs_free_energy, precision)


@pytest.mark.parametrize("QH, E, ZPE, H, TS, TqhS, G, qhG", [
    ('grimme',  -242.328708, 0.000000, -242.327125, 0.011221, 0.011221, -242.338346, -242.338346),
    ('truhlar', -242.328708, 0.000000, -242.327125, 0.011221, 0.011221, -242.338346, -242.338346)
])
def test_temperature_corrections(QH, E, ZPE, H, TS, TqhS, G, qhG):
    temp = 200
    conc = GV.ATMOS / (GV.GAS_CONSTANT * temp)
    freq_cutoff, freq_scale_factor, solv, spc = 100.0, 1.0, 'none', False
    bbe298 = GV.calc_bbe(datapath('Al_298K.out'), QH, freq_cutoff, temp, conc, freq_scale_factor, solv, spc)
    bbe400 = GV.calc_bbe(datapath('Al_400K.out'), QH, freq_cutoff, temp, conc, freq_scale_factor, solv, spc)
    precision = 6
    assert E == round(bbe298.scf_energy, precision) == round(bbe400.scf_energy, precision)
    assert ZPE == round(bbe298.zpe, precision) == round(bbe400.zpe, precision)
    assert H == round(bbe298.enthalpy, precision) == round(bbe400.enthalpy, precision)
    assert TS == round(temp * bbe298.entropy, precision) == round(temp * bbe400.entropy, precision)
    assert TqhS == round(temp * bbe298.qh_entropy, precision) == round(temp * bbe400.qh_entropy, precision)
    assert G == round(bbe298.gibbs_free_energy, precision) == round(bbe400.gibbs_free_energy, precision)
    assert qhG == round(bbe298.qh_gibbs_free_energy, precision) == round(bbe400.qh_gibbs_free_energy, precision)


@pytest.mark.parametrize("spc, E_spc, E, ZPE, H, TS, TqhS, GT, qhGT", [
    (False,        None, -79.830421, 0.075238, -79.750770, 0.027523, 0.027525, -79.778293, -79.778295),
    ('link', -79.830421, -79.830421, 0.075238, -79.750770, 0.027523, 0.027525, -79.778293, -79.778295),
    ('spc',  -79.858399, -79.830421, 0.075238, -79.778748, 0.027523, 0.027525, -79.806271, -79.806273),
    ('TZ',   -79.858399, -79.830421, 0.075238, -79.778748, 0.027523, 0.027525, -79.806271, -79.806273)
])
def test_single_point_correction(spc, E_spc, E, ZPE, H, TS, TqhS, GT, qhGT):
    temp = 298.15
    conc = GV.ATMOS / (GV.GAS_CONSTANT * temp)
    QH, freq_cutoff, freq_scale_factor, solv = 'grimme', 100.0, 1.0, 'none'
    precision = 6

    bbe = GV.calc_bbe(datapath('ethane.out'), QH, freq_cutoff, temp, conc, freq_scale_factor, solv, spc)
    if E_spc:
        assert E_spc == round(bbe.sp_energy, precision)
    assert E == round(bbe.scf_energy, precision)
    assert ZPE == round(bbe.zpe, precision)
    assert H == round(bbe.enthalpy, precision)
    assert TS == round(temp * bbe.entropy, precision)
    assert TqhS == round(temp * bbe.qh_entropy, precision)
    assert GT == round(bbe.gibbs_free_energy, precision)
    assert qhGT == round(bbe.qh_gibbs_free_energy, precision)


@pytest.mark.parametrize("filename, freq_scale_factor, zpe", [
    ('ethane.out', 0.977, 0.073508)
])
def test_scaling_factor_search(filename, freq_scale_factor, zpe):
    temp = 298.15
    conc = GV.ATMOS / (GV.GAS_CONSTANT * temp)
    QH, freq_cutoff, solv, spc = 'grimme', 100.0, 'none', False
    precision = 6

    bbe = GV.calc_bbe(datapath('ethane.out'), QH, freq_cutoff, temp, conc, freq_scale_factor, solv, spc)
    assert zpe == round(bbe.zpe, precision)
