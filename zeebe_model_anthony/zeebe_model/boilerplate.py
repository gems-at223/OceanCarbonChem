"""
Functions for interfacing with Zeebe's microenvironment model
"""

import os
import numpy as np
import pandas as pd
from glob import glob

# import pkg_resources as pkgrs
# resource_dir = pkgrs.resource_filename('zeebe_model', 'resources')

resource_dir = os.path.join(os.path.split(__file__)[0], "resources")


# data import
def import_modelrun(folder="."):
    """
    Imports all .sv4 files in the directory as a dict.
    Concentrations are in µM.
    """
    sv4s = glob(folder + "/*.sv4")
    if len(sv4s) == 0:
        raise ValueError(f"No output (.sv4) files in foler {folder}")

    data = {
        os.path.basename(f).replace(".sv4", ""): np.genfromtxt(f)
        for f in sv4s
        if "par" not in f
    }
    data["pH"] = -np.log10(data["h"] * 1e-6)
    data["dic"] = data["co2"] + data["co3"] + data["hco3"]
    data["alk"] = data["hco3"] + 2 * data["co3"] + data["boh4"] + data["oh"]

    with open(folder + "/par.sv4") as f:
        data["par"] = f.read()

    return data


def parse_modelrun(data):
    """
    Parse model run into a pandas DataFrame and a dict of metadata.

    Parameters
    ----------
    data : str or dict
        Either a dictionary output by import_modelrun or a path
        to a folder that can be parsed by import_modelrun.

    Returns
    -------
    tuple : (pd.DataFrame, dict)
    """
    if isinstance(data, str):
        data = import_modelrun(data)

    n = len(data["r"])
    meta = {}
    d = {}
    for k, v in data.items():
        if len(v) == n:
            d[k] = v
        else:
            meta[k] = v

    return pd.DataFrame.from_dict(d).set_index("r"), meta


def c_run(path):
    # open('./a.out', 'a').close()
    os.system("gcc " + path + " -lm; ./a.out")


def make_runfile(params, outpath="./py_run.c", template=None, itmax=100, slowc=0.3):
    if isinstance(template, str):
        with open(template, "r") as f:
            template = f.read()

    for k, v in params.items():
        template = template.replace("**" + k + "**", "{:.9e}".format(v))

    template = template.replace("**ITMAX**", "{:.0f}".format(itmax))
    template = template.replace("**SLOWC**", "{:.2f}".format(slowc))

    with open(outpath, "w") as f:
        f.write(template)


def cp_nrutil(tpath="./py_run/"):
    if not os.path.exists(tpath + "/nrutil.c"):
        os.system(
            f"cp {os.path.join(resource_dir, 'nrutil.c')} {os.path.join(tpath, 'nrutil.c')}"
        )


# parameter handling
def make_params(
    RADIUS,
    CO3UPT,
    CO2UPT,
    HCO3UPT,
    PHBULK,
    DICBULK,
    UALKBULK,
    BORMULT,
    SYMCO2UPT,
    SYMHCO3UPT,
    SYMTCUPT,
    VMAX,
    SYMDIST,
    REDS,
    SALINITY,
    TEMP,
):
    params = {
        "RADIUS": RADIUS,  # Radius of foram
        "CO3UPT": CO3UPT,  # foram CO3 uptake
        "CO2UPT": CO2UPT,  # foram CO2 uptake
        "HCO3UPT": HCO3UPT,  # foram HCO3 uptake
        "PHBULK": PHBULK,  # pH of bulk solution
        "DICBULK": DICBULK,  # DIC of bulk solution
        "UALKBULK": UALKBULK,  # Alkalinity of bulk solution
        "SYMCO2UPT": SYMCO2UPT,  # symbiont CO2 uptake
        "SYMHCO3UPT": SYMHCO3UPT,  # symbiont CO3 uptake
        "SYMTCUPT": SYMTCUPT,  # symbiont total C uptake
        "VMAX": VMAX,  # Michaelis-Menten Vmax
        "SYMDIST": SYMDIST,  # Distance of symbionts from shell surface
        "REDS": REDS,  # Redfield ratio of symbionts
        "SALINITY": SALINITY,  # Salinity... obv
        "TEMP": TEMP,
    }  # Temp ºC
    params["BORTBULK"] = BORMULT * (
        416 * params["SALINITY"] / 35
    )  # boron total in solution, adjusted for salinity

    return params


# run model
def run(
    params,
    modelname="solvde42_py_run.c",
    tpath="./py_run/",
    template=None,
    itmax=400,
    slowc=0.3,
):
    """
    Runs the model with the given parameter dict.

    Parameters
    ----------
    params : dict
        dict of params to pass to the model
    itmax : int
        the maximum number of iterations run by the model
    slowc : float
        A parameter that controls the step-size taken during
        iteration. If you're getting a lot of oscillations in
        the boundary layer, try decreasing this number. IF you
        decrease this number, you should also increase itmax
        as convergence will be slower.
    """
    if not os.path.exists(tpath):
        os.mkdir(tpath)

    if template is None:
        template = os.path.join(resource_dir, "solvde42_py_temp.c")
    make_runfile(
        params, os.path.join(tpath, modelname), template, itmax=itmax, slowc=slowc
    )

    cp_nrutil(tpath)

    curdir = os.getcwd()

    os.chdir(tpath)
    c_run(modelname)
    os.chdir(curdir)

    return parse_modelrun(tpath)
