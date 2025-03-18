import marimo

__generated_with = "0.11.13"
app = marimo.App(width="medium")


@app.cell
def _():
    import os
    import matplotlib

    import matplotlib.pyplot as plt
    import marimo as mo
    import numpy as np
    import pandas as pd
    import pyfraglib.scores as pfs

    matplotlib.style.use("default")
    return matplotlib, mo, np, os, pd, pfs, plt


@app.cell
def _(os, pd):
    DATA_DIR = "/mnt/ramses/scratch/dschuet7/nf_out/"
    INFO_DIR = "/mnt/ramses/projects/uk-lymphoma-cfdna/PCNSL/"

    _sample_name = "DED005_BL"
    _wps_filepath = os.path.join(DATA_DIR, _sample_name, f"wps_{_sample_name}.csv")
    DED005_BL_wps_df = pd.read_csv(_wps_filepath, low_memory=False)
    return DATA_DIR, DED005_BL_wps_df, INFO_DIR


@app.cell
def _(DED005_BL_wps_df, pfs):
    pfs.score_line_plot(DED005_BL_wps_df, name="DED005_BL_MYD88", out_dir="./", region=(38_180_125-100, 38_182_819+100),
                        exclude_chroms=["1", "2"] + [str(i) for i in range(4, 23)] + ["X", "Y", "M"],
                        score="ratio_span_total", log_transform=False)
    return


@app.cell
def _(DATA_DIR, os, pd):
    _sample_name = "DED006_BL"
    _wps_filepath = os.path.join(DATA_DIR, _sample_name, f"wps_{_sample_name}.csv")
    DED006_BL_wps_df = pd.read_csv(_wps_filepath, low_memory=False)
    return (DED006_BL_wps_df,)


@app.cell
def _(DED006_BL_wps_df, pfs):
    pfs.score_line_plot(DED006_BL_wps_df, name="DED006_BL_MYD88", out_dir="./", region=(38_180_125-100, 38_182_819+100),
                        exclude_chroms=["1", "2"] + [str(i) for i in range(4, 23)] + ["X", "Y", "M"],
                        score="ratio_span_total", log_transform=False)
    return


@app.cell
def _(DATA_DIR, os, pd):
    _sample_name = "ctrl1"
    _wps_filepath = os.path.join(DATA_DIR, _sample_name, f"wps_{_sample_name}.csv")
    ctrl1_wps_df = pd.read_csv(_wps_filepath, low_memory=False)
    return (ctrl1_wps_df,)


@app.cell
def _(ctrl1_wps_df, pfs):
    pfs.score_line_plot(ctrl1_wps_df, name="ctrl1_MYD88", out_dir="./", region=(38_180_125-100, 38_182_819+100),
                        exclude_chroms=["1", "2"] + [str(i) for i in range(4, 23)] + ["X", "Y", "M"],
                        score="ratio_span_total", log_transform=False)
    return


@app.cell
def _(DATA_DIR, os, pd):
    _sample_name = "ctrl2"
    _wps_filepath = os.path.join(DATA_DIR, _sample_name, f"wps_{_sample_name}.csv")
    ctrl2_wps_df = pd.read_csv(_wps_filepath, low_memory=False)
    return (ctrl2_wps_df,)


@app.cell
def _(ctrl2_wps_df, pfs):
    pfs.score_line_plot(ctrl2_wps_df, name="ctrl2_MYD88", out_dir="./", region=(38_180_125-100, 38_182_819+100),
                        exclude_chroms=["1", "2"] + [str(i) for i in range(4, 23)] + ["X", "Y", "M"],
                        score="ratio_span_total", log_transform=False)
    return


if __name__ == "__main__":
    app.run()
