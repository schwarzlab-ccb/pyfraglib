import marimo

__generated_with = "0.10.14"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _():
    import json
    import os

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    from sksurv.nonparametric import kaplan_meier_estimator
    from lifelines import KaplanMeierFitter
    return (
        KaplanMeierFitter,
        json,
        kaplan_meier_estimator,
        mo,
        np,
        os,
        pd,
        plt,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
        # `pyfraglib` applied to the PCNSL cohort
        In this notebook, we will start to investigate the readouts of `pyfraglib` with regards to their clinical applicability.

        First, we will analyze the fitted parameters of our Gaussian mixture model in the context
        of the primary central nervous system lymphoma (PCNSL) cohort from AG Borchmann.
        The feature vectors available to us are ($n = \\text{number of Gaussians}$):

        $(\\sigma_1, \\sigma_2, ..., \\sigma_n)$
        $(\\mu_1, \\mu_2, ..., \\mu_n)$
        $(\\pi_1, \\pi_2, ..., \\pi_n)$

        Please refer to `pyfraglib` for more explanation on these parameters in the context of the fitted GMM.

        In our workflow, we are `ssh`-mounting the RAMSES HPC cluster. Data is then loaded from the mount point and visually explored.
        """
    )
    return


@app.cell(hide_code=True)
def _(os, pd):
    DATA_DIR = "/mnt/ramses/scratch/dschuet7/nf_out/"
    INFO_DIR = "/mnt/ramses/projects/uk-lymphoma-cfdna/PCNSL/"
    sample_sheet = pd.read_csv(os.path.join(INFO_DIR, "sample_sheet.csv"))
    clin_info_sheet = pd.read_excel(os.path.join(INFO_DIR, "clinical_annotations.xlsx"))

    clin_info_sheet_column_metadata = clin_info_sheet.columns
    clin_info_sheet.columns = clin_info_sheet.iloc[0]
    clin_info_sheet = clin_info_sheet[1:].reset_index(drop=True)
    clin_info_sheet = clin_info_sheet.iloc[:, 1:]
    return (
        DATA_DIR,
        INFO_DIR,
        clin_info_sheet,
        clin_info_sheet_column_metadata,
        sample_sheet,
    )


@app.cell(hide_code=True)
def _(DATA_DIR, json, os, pd, sample_sheet):
    data: dict[str, object] = {}
    for id, timepoint, timepoint_hom in zip(sample_sheet["sample_id"], sample_sheet["time_point"], sample_sheet["time_point_hom"]):
        sample_name = f"{id}_{timepoint}"
        model_params_filepath = os.path.join(DATA_DIR, sample_name, f"{sample_name}_gmm_frags_len.json")
        with open(model_params_filepath, "r") as f:
            json_data = json.load(f)
            assert json_data
            assert sample_name not in data, f"{sample_name} is already in {data.keys()}"

            json_data["time_point"] = timepoint
            json_data["time_point_hom"] = timepoint_hom
            data[sample_name] = json_data

    model_params = pd.DataFrame(data).T
    return (
        data,
        f,
        id,
        json_data,
        model_params,
        model_params_filepath,
        sample_name,
        timepoint,
        timepoint_hom,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
        ## Investigating fitted model parameters
        We are interested in how well our fitted GMMs approximate the fragment length distributions they are supposed to model. That's not an easy question to answer, but we have a few metrics at hand that `pyfraglib` calculated for us. Since the sample size is huge and all of these metrics are somewhat sensitive to that fact (i.e. they detect even minor deviations of the empirical distribution from the fit if too many samples are considered), they should be taken with a grain of salt.

        The first plot down below shows the final value of the objective function per sample on a $\\log_e$ scale. That's useful because we are not optimizing the negative log-likelihood but

        \[
            argmin_{\\sigma_i,\\mu_i,\\pi_i}{\\frac{\\text{NLL}_{\\sigma_i,\\mu_i,\\pi_i}}{\\text{NLL}_{\\sigma_{ini},\\mu_{ini},\\pi_{ini}}}}
        \]

        where $x_{ini}$ is the initial guess of the respective model parameter.
        Thus, we get an estimate of how far away from the initial guesses the fitted parameters / final value of the objective function ends up.
        """
    )
    return


@app.cell(hide_code=True)
def _(model_params, np, plt):
    color_map = {"BL": "red", "BLrr": "darkred", "c1": "blue", "c2": "darkblue", "end": "yellow", "CSF": "green"}
    color_code = [color_map[tp] for tp in model_params["time_point_hom"]]

    _, axis = plt.subplots(figsize=(14,5))
    axis.bar(x=model_params.index,
             height=np.log(model_params["objective_value"].astype(float)),
             color=color_code)
    axis.tick_params(axis="x", rotation=90, labelsize=6)
    axis.set_title(r"$\log_e(\text{obj-func})$")
    handles = [plt.Rectangle((0, 0), 1, 1, color=color_map[label]) for label in color_map]
    axis.legend(handles, color_map.keys(), title="Timepoints")
    axis
    return axis, color_code, color_map, handles


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Next, we are plotting the fitted parameters:""")
    return


@app.cell
def _(clin_info_sheet, clin_info_sheet_column_metadata):
    for idx, item in enumerate(clin_info_sheet_column_metadata):
        print(f"{idx}: {item}")
    clin_info_sheet.iloc[:, 70:]
    return idx, item


@app.cell
def _(KaplanMeierFitter, clin_info_sheet, kaplan_meier_estimator, np, plt):
    days_diff = clin_info_sheet.loc[:, "last_FU"] - clin_info_sheet.loc[:, "ED_date"]
    days_diff = np.array(days_diff)
    times = [d.days for d in days_diff]
    events = clin_info_sheet.loc[:, "survival"] == 1

    def km_plot_scikit(events, times):
        surv_time, surv_prob = kaplan_meier_estimator(events, times)
        _, _axis = plt.subplots()
        _axis.step(surv_time, surv_prob, where="post")
        _axis.set_ylim((0.0, 1.05))
        return _axis

    def km_plot_lifelines(events, times):
        kmf = KaplanMeierFitter()
        kmf.fit(times, events)
        from math import isinf
        print("mOS not reached" if isinf(kmf.median_survival_time_)
              else f"mOS {kmf.median_survival_time_} days")
        
        _axis = kmf.plot_survival_function(at_risk_counts=True, show_censors=True)
        _axis.set_ylim((0.0, 1.05))
        _axis.set_xlabel("Time (days)")
        _axis.set_ylabel("Overall survival probability")
        return _axis

    km_plot_lifelines(events, times)
    return days_diff, events, km_plot_lifelines, km_plot_scikit, times


if __name__ == "__main__":
    app.run()
