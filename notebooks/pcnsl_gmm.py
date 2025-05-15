import marimo

__generated_with = "0.11.13"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _():
    import json
    import os

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns
    import pandas as pd

    from math import isinf
    from sksurv.nonparametric import kaplan_meier_estimator
    from scipy.stats import mannwhitneyu, ttest_ind
    from lifelines import KaplanMeierFitter
    from lifelines.plotting import add_at_risk_counts
    from lifelines.statistics import logrank_test
    return (
        KaplanMeierFitter,
        add_at_risk_counts,
        isinf,
        json,
        kaplan_meier_estimator,
        logrank_test,
        mannwhitneyu,
        mo,
        np,
        os,
        pd,
        plt,
        sns,
        ttest_ind,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
        # `pyfraglib` applied to the PCNSL cohort
        In this notebook, we will investigate the readouts of `pyfraglib` with regards to their clinical applicability.

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
    DATA_DIR = "/mnt/ramses/projects/uk-lymphoma-cfdna/PCNSL/frag_out_26032025/"
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
            argmin_{\\sigma_i,\\mu_i,\\pi_i}{\\frac{\\text{NLL}_{\\sigma_i,\\mu_i,\\pi_i}}{|\\text{NLL}_{\\sigma_{ini},\\mu_{ini},\\pi_{ini}}|}}
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
    handles = [plt.Rectangle((0, 0), 1, 1, color=color_map[label]) for label in color_map]

    _, _axis = plt.subplots(figsize=(14,5))
    _axis.bar(x=model_params.index,
             height=np.log(model_params["objective_value"].astype(float)),
             color=color_code)
    _axis.tick_params(axis="x", rotation=90, labelsize=6)
    _axis.set_ylabel(r"$\log_e(\text{obj-func})$")
    _axis.set_xlabel("Samples")
    _axis.legend(handles, color_map.keys(), title="Timepoints")
    _axis
    return color_code, color_map, handles


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Next, we are plotting the goodness-of-fit parameters and the model parameters:""")
    return


@app.cell(hide_code=True)
def _(color_code, color_map, handles, model_params, np, plt):
    _, _axis = plt.subplots(figsize=(14,5))
    _axis.bar(x=model_params.index,
             height=np.log10(model_params["kolmogorov_smirnov_statistic"].astype(float)),
             color=color_code)
    _axis.tick_params(axis="x", rotation=90, labelsize=6)
    _axis.set_ylabel(r"$\log_{10}(\text{KS Test Statistic})$")
    _axis.set_xlabel("Samples")
    _axis.legend(handles, color_map.keys(), title="Timepoints")
    _axis
    return


@app.cell(hide_code=True)
def _(color_code, color_map, handles, model_params, plt):
    _, _axis = plt.subplots(nrows=3, ncols=1, figsize=(14,10))

    for idx in [0, 1, 2]:
        _axis[idx].bar(x=model_params.index,
                          height=[x[idx] for x in model_params["estimated_means"]],
                          color=color_code)
        _axis[idx].tick_params(axis="x", rotation=90, labelsize=6)
        _axis[idx].set_ylabel(f"$\mu_{idx+1}$")
        _axis[idx].set_xlabel("Samples")
        _axis[idx].legend(handles, color_map.keys(), title="Timepoints")

    plt.tight_layout()
    plt.show()
    return (idx,)


@app.cell(hide_code=True)
def _(color_code, color_map, handles, model_params, plt):
    _, _axis = plt.subplots(nrows=3, ncols=1, figsize=(16,10))

    for _idx in [0, 1, 2]:
        _axis[_idx].bar(x=model_params.index,
                          height=[x[_idx] for x in model_params["estimated_pis"]],
                          color=color_code)
        _axis[_idx].tick_params(axis="x", rotation=90, labelsize=6)
        _axis[_idx].set_ylabel(f"$\pi_{_idx+1}$")
        _axis[_idx].set_xlabel("Samples")
        _axis[_idx].legend(handles, color_map.keys(), title="Timepoints")

    plt.tight_layout()
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
        ## Analysis of clinical cohort characteristics
        Finally, we want to relate the fitted model parameters to patients' clinical data. That is mainly survival and toxicities of chemotherapy. Along the way, we might find additional characteristics that are interesting to predict using our fragmentomics features.
        """
    )
    return


@app.cell(hide_code=True)
def _(
    KaplanMeierFitter,
    clin_info_sheet,
    isinf,
    kaplan_meier_estimator,
    np,
    pd,
    plt,
):
    days_diff = clin_info_sheet.loc[:, "last_FU"] - clin_info_sheet.loc[:, "ED_date"]
    days_diff = np.array(days_diff)
    times = pd.Series([d.days for d in days_diff])
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
        print("mOS not reached" if isinf(kmf.median_survival_time_)
              else f"mOS {kmf.median_survival_time_} days")

        _axis = kmf.plot_survival_function(at_risk_counts=True, show_censors=True)
        _axis.set_ylim((0.0, 1.05))
        _axis.set_xlabel("Time (days)")
        _axis.set_ylabel("Overall survival probability")
        return _axis

    # km_plot_lifelines(events, times)
    return days_diff, events, km_plot_lifelines, km_plot_scikit, times


@app.cell(hide_code=True)
def _(
    KaplanMeierFitter,
    add_at_risk_counts,
    clin_info_sheet,
    events,
    logrank_test,
    plt,
    times,
):
    _, _axis = plt.subplots()
    strata = clin_info_sheet["age"] > 60
    stratum_true_label = "Age > 60"
    stratum_false_label = "Age <= 60"

    kmf_train = KaplanMeierFitter()
    kmf_train.fit(times[strata], events[strata], label=stratum_true_label)
    kmf_train.plot_survival_function(ax=_axis, show_censors=True)

    kmf_val = KaplanMeierFitter()
    kmf_val.fit(times[~strata], events[~strata], label=stratum_false_label)
    kmf_val.plot_survival_function(ax=_axis, show_censors=True)

    logrank_pval = logrank_test(times[strata], times[~strata], events[strata], events[~strata], alpha=.99).p_value

    _axis.set_ylim((0.0, 1.05))
    _axis.set_xlabel("Time (days)")
    _axis.set_ylabel("Overall survival probability")
    _axis.set_title(f"log-rank p={logrank_pval:0.2}")
    add_at_risk_counts(kmf_train, kmf_val, ax=_axis)
    return (
        kmf_train,
        kmf_val,
        logrank_pval,
        strata,
        stratum_false_label,
        stratum_true_label,
    )


@app.cell(hide_code=True)
def _(clin_info_sheet, model_params, pd):
    model_params.index = [name.split("_")[0] for name in model_params.index]
    model_clin_combined: pd.DataFrame = model_params.merge(right=clin_info_sheet, left_index=True, right_on="study_ID", how="inner")
    return (model_clin_combined,)


@app.cell(hide_code=True)
def _(mannwhitneyu, model_clin_combined, pd, plt, sns, ttest_ind):
    _df = model_clin_combined
    PEAK: int = 2
    plot_data: pd.DataFrame = pd.DataFrame(
        [
            ("BL" if tp == "BLrr" else tp,
             "dead" if surv == 1 else "alive",
             "relapse" if rel == 1 else "no relapse",
             pis[PEAK-1],
             sample_id,
             "yes" if csf_protein else "no",
             "yes" if ldh_elevated else "no",
             "yes" if multiple_lesions else "no",
             "low" if cfdna_conc < 13.68 else "high",
             "low" if ctdna_conc < 1.63 else "high")
            for (tp, surv, rel, pis, sample_id, csf_protein, ldh_elevated, multiple_lesions, cfdna_conc, ctdna_conc)
            in zip(_df["time_point_hom"], _df["survival"], _df["relapse"], _df["estimated_pis"], _df["study_ID"],
                   _df["CSF_protein"], _df["LDH_elevated"], _df["multiple_lesions"], _df["cfDNA"], _df["ctDNA"])
            if tp not in ["CSF", "c2"]
        ],
        columns=["time_point", "survival", "relapse", "pi_2", "sample_id", "csf_protein", "ldh_elevated",
                 "multiple_lesions", "cfDNA_concentration", "ctDNA_concentration"]
    )
    _bl = plot_data[plot_data["time_point"] == "BL"]["pi_2"]
    _end = plot_data[plot_data["time_point"] == "end"]["pi_2"]
    mwu_pval = mannwhitneyu(_bl, _end, alternative="two-sided").pvalue
    tt_pval = ttest_ind(_bl, _end, alternative="two-sided").pvalue

    plt.figure(figsize=(8, 6))
    sns.violinplot(data=plot_data, x="time_point", y="pi_2", hue="ctDNA_concentration", order=["BL", "c1", "end"])
    plt.axhline(y=0.1, color="red", alpha=0.6, linestyle="--", linewidth=2)
    plt.xlabel("Timepoints")
    plt.ylabel(f"$\pi_{PEAK}$")
    plt.title(f"Comparison of GMM $\pi_{PEAK}$ across treatment timepoints in PCNSL cohort"
              f"\nMann-Whitney-U BL vs. end p={mwu_pval:0.2}")
    return PEAK, mwu_pval, plot_data, tt_pval


@app.cell
def _(PEAK, pd, plot_data, sns):
    plot_data["time_point"] = pd.Categorical(plot_data["time_point"],
                                             categories=["BL", "c1", "end"],
                                             ordered=True)
    #sns.lineplot(data=plot_data, x="time_point", y="pi_2", units="sample_id", estimator=None,
    #             hue="survival", marker="o")
    _plot = sns.relplot(
        data=plot_data,
        x="time_point",
        y="pi_2",
        units="sample_id",
        hue="relapse",
        style="relapse",
        col="ctDNA_concentration",
        kind="line",
        estimator=None,
        linewidth=2,
        alpha=0.6,
        facet_kws={"sharey": True},
    )
    _plot.set_xlabels("Timepoints")
    _plot.set_ylabels(f"$\pi_{PEAK}$")
    return


@app.cell(hide_code=True)
def _(model_clin_combined, pd, plt, sns):
    _df = model_clin_combined
    _df = _df[_df["time_point_hom"] == "BL"]

    _data: pd.DataFrame = pd.DataFrame(
        [("dead" if surv == 1 else "alive", pis[1]) for (surv, pis) in zip(_df["survival"], _df["estimated_pis"])],
        columns=["survival", "pi_2"]
    )
    plt.figure(figsize=(8, 6))
    sns.violinplot(data=_data, x="survival", y="pi_2")
    plt.xlabel("Patient survival")
    plt.ylabel(r"$\pi_2$")
    plt.title(r"Comparison of GMM $\pi_2$ across BL samples in PCNSL cohort")
    return


@app.cell(hide_code=True)
def _(model_clin_combined, pd, plt, sns):
    _df = model_clin_combined
    _df = _df[_df["time_point_hom"] == "end"]

    _data: pd.DataFrame = pd.DataFrame(
        [("dead" if surv == 1 else "alive", pis[1]) for (surv, pis) in zip(_df["survival"], _df["estimated_pis"])],
        columns=["survival", "pi_2"]
    )
    plt.figure(figsize=(8, 6))
    sns.violinplot(data=_data, x="survival", y="pi_2")
    plt.xlabel("Patient survival")
    plt.ylabel(r"$\pi_2$")
    plt.title(r"Comparison of GMM $\pi_2$ across end-of-treatment samples in PCNSL cohort")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
        # Todo
        - Build quotient of $\\frac{end}{BL}$ $\pi_3$ values. Correlate with clinical outcome.
        - Sankey plots of timepoints and metrics
        - Cutoff for mixture fraction?
        - repeat mixture fraction analysis with means and standard deviations
        - write statistics to JSON and plot them
        - explore WPS and possible ways to use it to predict patient outcome
        - correlate pi's with toxicity metrics (cummulative dose of methotrexate?)
        """
    )
    return


if __name__ == "__main__":
    app.run()
