# =========================================================================
# VERSION INFO
#	Last Modified 	--- 2022-05-22
#	Version 		--- 1.0
# -------------------------------------------------------------------------
# INTRODUCTION
#	- Calculate the pure-tone attenuation coefficient due to the atmospheric
#       absorption based on the standard ISO 9613-1.
# -------------------------------------------------------------------------
# REFERENCES
#	[1] ISO Technical Committees: Noise. Acoustics — Attenuation of sound
#		during propagation outdoors — Part 1: Calculation of the absorption
#		of sound by the atmosphere: ISO 9613-1:1993[S]. Geneva:
#		International Organization for Standardization, 1993.
#	[2]	National Physical Laboratory. NPL Acoustics: Calculation of
#		absorption of sound by the atmosphere[EB/OL]. [2018-08-08].
#		http://resource.npl.co.uk/acoustics/techguides/absorption/.
# -------------------------------------------------------------------------
# INPUT
#	freq		- Frequency, in Hertz.
# OPTIONS
#	humidity 	- relative humidity in percentage
#				- dafault: 60
#	temperature	- Ambient atmospheric temperature, in Celcius
#				- default: 25
#	pressure 	- the atmospheric pressure, in kilopscals
#				- default: 101.325
# NOTES
#	- the dimension of all inputs must be compatible
# -------------------------------------------------------------------------
# OUTPUT
#	alpha_Np	- pure-tone sound attenuation coeff. in Neper per meter, for
#					atmospheric absorption
#	alpha_dB	- pure-tone sound attenutaion oeffi. in dB per meter, for
#					atmospheric absorption
# =========================================================================
import numpy as np


def AbsorpAttenCoef(freq, temperature=20, pressure=101.325, humidity=70):
    # 后三个为可选参数及默认值，单位分别为Celsius,kPa,%
    T0 = 293.15  # reference air temperature in Kelvins (20 Celsius)
    T01 = 273.16  # the triple-point isotherm temperature (+0.01 Celsius)
    T = temperature + 273.15  # Ambient atmospheric temperature in Kelvins

    pr = 101.325  # reference ambient atmospheric pressure in kilopascals

    C = -6.8346 * (T01 / T) ** 1.261 + 4.6151
    psat = pr * 10 ** C
    h = humidity * (psat / pr) * (pressure / pr)

    f_rO = pressure / pr * (24 + 4.04 * 10 ** 4 * h * (0.02 + h) / (0.391 + h))
    f_rN = pressure / pr * (T / T0) ** (-1 / 2) * (9 + 280 * h * np.exp(-4.17 * ((T / T0) ** (-1 / 3) - 1)))
    alpha_Np = freq ** 2 * (1.84 * 10 ** (-11) * pr / pressure * (T / T0) ** (-1 / 2) + (T / T0) ** (-5 / 2) * (
            0.01275 * np.exp(-2239.1 / T) / (f_rO + freq ** 2 / f_rO) + 0.1068 * np.exp(-3352.0 / T) / (
            f_rN + freq ** 2 / f_rN)))
    alpha_dB = 20 / np.log(10) * alpha_Np

    return alpha_Np, alpha_dB
