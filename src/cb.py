import cbsyst as cb
import numpy as np

K_values = {
    "K1": 1.0e-14,
    "K2": 4.8e-14,
    "KB": 4.8e-14,
    "KW": 1.0e-14,
    "KS": 1.0e-14,
    "KF": 1.0e-14,
    "K0": 1.0e-14,
    "KP1": 1.0e-14,
    "KP2": 1.0e-14,
    "KP3": 1.0e-14,
    "KSi": 1.0e-14,
    "KspA": 1.0e-14,
    "KspC": 1.0e-14,
}

system = cb.CBsys(pHtot=8.25, CO2=9.3, T_in=24.5, S_in=40.7, Ks=K_values)
print(system)
