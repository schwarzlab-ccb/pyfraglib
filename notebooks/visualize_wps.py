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
def _():
    DATA_DIR = "/mnt/ramses/scratch/dschuet7/nf_out/"
    INFO_DIR = "/mnt/ramses/projects/uk-lymphoma-cfdna/PCNSL/"
    return DATA_DIR, INFO_DIR


@app.cell
def _(DATA_DIR, os, pd, pfs):
    for _sample_name in ["ctrl1", "ctrl2", "ctrl3", "ctrl4", "ctrl5", "ctrl6", "ctrl7", "ctrl8", "ctrl9", "ctrl10",
                         "DED005_BL", "DED006_BL", "DED006DRR049_BL", "DED006DRR049_BLrr",
                         "DED006DRR049_M3", "DED006DRR049_W3", "DED006DRR049_W6", "DED011_BL",
                         "DED017DRR058_BL", "DED017DRR058_M3", "DED017DRR058_W3", "DED017DRR058_W6",
                         "DED021_BL", "DED022_BL", "DED024_BL", "DED024_W6", "DED025_BL", "DED025_M1",
                         "DED029_BL", "DED031_BL", "DED034_BL", "DED044_BL", "DED045_BL", "DED048_BL",
                         "DED049_BL", "DED050_BL", "DED055_BL", "DED057_BL", "DED069_BL", "DED071_BL",
                         "DED072_BL", "DED078_BL", "DED078_M3", "DED082_BL", "DED082_M3", "DED084_BL",
                         "DED084_M3", "DED086_BL", "DED086_BLrr", "DED086_M1", "DED086_M3", "DED090DRR077_BL",
                         "DED090DRR077_BLrr", "DED090DRR077_M1", "DED090DRR077_M3", "DED092_BL", "DED092_M1",
                         "DED092_M3", "DED094_BL", "DED094_M1", "DED094_M6", "DED096DRR098_BL", "DED096DRR098_M3",
                         "DED096DRR098_W3", "DED105_BL", "DED105_M1", "DED105_M3", "DED106_BL", "DED106_M1",
                         "DED108_BL", "DED108_M1", "DED108_M3", "DED109_BL", "DED109_CSF", "DED109_M1",
                         "DED109_M3", "DED110_BL", "DED110_M1", "DED110_M3", "DED117_BL", "DED117_W4", "DED118_BL",
                         "DED118_CSF", "DED118_M3", "DED118_W4", "DED119_BL", "DED119_M1", "DED120_BL",
                         "DED121_BL", "DED123_BL", "DED123_CSF", "DED123_M1", "DED123_M3", "DED126_BL",
                         "DED126_CSF", "DED126_M1", "DED126_M3", "DED128_BL", "DED128_M1", "DED128_M3",
                         "DED129_BL", "DED129_M1", "DED131_BL", "DED131_CSF", "DED131_M3", "DED131_W3",
                         "DED131rr_BL", "DED131rr_M3", "DED131rr_W3", "DED132_BL", "DED132_CSF", "DED132_M1",
                         "DED132_M3", "DED133_BL", "DED133_M1", "DED134_BL", "DED134_M1", "DED134_M3", "DED135_BL",
                         "DED136_BL", "DED136_W3", "DED137_BL", "DED137_CSF", "DED138_BL", "DED139_BL",
                         "DED139_W3", "DED140_BL", "DED140_M1", "DED140_M3", "DED140DRR098_BL", "DED140DRR098_M1",
                         "DED142_BL", "DED143_BL", "DED143_M1", "DED143_M3", "DED143_W6", "DED144_BL",
                         "DED145_M1", "DED145_M3", "DED146_BL", "DED146_W2", "DED146_W6", "DED147_BL",
                         "DED147_W3", "DED147_W6", "DED148_BL", "DED148_M1", "DED148_M3", "DED149_BL",
                         "DED149_M1", "DED149_M3", "DED150_BL", "DED150_M1", "DED151_BL", "DRR017_BL",
                         "DRR017_W3", "DRR017_W6", "DRR071_BL", "DRR071_M1", "DRR074_BL", "DRR074_M3",
                         "DRR074_W3", "DRR074_W6", "DRR076_BL"]:
        _wps_filepath = os.path.join(DATA_DIR, _sample_name, f"wps_{_sample_name}.csv")
        _wps_df = pd.read_csv(_wps_filepath, low_memory=False)
        try:
            pfs.score_line_plot(_wps_df, name=f"{_sample_name}_MYD88",
                                out_dir="./", region=(38_180_125-100, 38_182_819+100),
                                exclude_chroms=["1", "2"] + [str(i) for i in range(4, 23)] + ["X", "Y", "M"],
                                score="ratio_span_total", log_transform=False)
        except Exception as _:
            print(f"{_wps_filepath} failed")
    return


if __name__ == "__main__":
    app.run()
