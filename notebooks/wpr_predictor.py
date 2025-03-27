import marimo

__generated_with = "0.11.13"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _():
    import os

    import matplotlib.pyplot as plt
    import marimo as mo
    import numpy as np
    import pandas as pd
    import polars as pl
    import plotly.express as px
    import seaborn as sns

    from sklearn.decomposition import PCA
    from sklearn.model_selection import train_test_split, cross_val_score
    from sklearn.linear_model import LogisticRegressionCV
    from sklearn.metrics import accuracy_score, roc_auc_score
    from sklearn.preprocessing import StandardScaler
    return (
        LogisticRegressionCV,
        PCA,
        StandardScaler,
        accuracy_score,
        cross_val_score,
        mo,
        np,
        os,
        pd,
        pl,
        plt,
        px,
        roc_auc_score,
        sns,
        train_test_split,
    )


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

    # @NOTE(ds): We need to copy 3 rows of the dataframe where we have patients that have DED _and_
    # DRR samples. One of them is even inconsistently named DEDxxxrr. Yuck.
    row_DED006 = clin_info_sheet.loc[clin_info_sheet["study_ID"] == "DED006DRR049"]
    row_DED006.loc[:, "study_ID"] = "DED006"
    row_DED131 = clin_info_sheet.loc[clin_info_sheet["study_ID"] == "DED131"]
    row_DED131.loc[:, "study_ID"] = "DED131rr"
    row_DED140DRR098 = clin_info_sheet.loc[clin_info_sheet["study_ID"] == "DED140"]
    row_DED140DRR098.loc[:, "study_ID"] = "DED140DRR098"

    clin_info_sheet = pd.concat([clin_info_sheet, row_DED006, row_DED131, row_DED140DRR098], ignore_index=True)
    return (
        DATA_DIR,
        INFO_DIR,
        clin_info_sheet,
        clin_info_sheet_column_metadata,
        row_DED006,
        row_DED131,
        row_DED140DRR098,
        sample_sheet,
    )


@app.cell(hide_code=True)
def _(DATA_DIR, os, pl, sample_sheet):
    def load_wpr_all_samples() -> (pl.Series, pl.DataFrame):
        reference_df: pl.DataFrame | None = None
        dfs: list[pl.DataFrame] = []

        all_samples: list[tuple[str, str, str]] = \
            list(zip([f"ctrl{i}" for i in range(1, 11)], [""]*10, [""]*10))
        all_samples += \
            list(zip(sample_sheet["sample_id"], sample_sheet["time_point"], sample_sheet["time_point_hom"]))

        for id, timepoint, timepoint_hom in all_samples:
            sample_name: str = f"{id}"
            if timepoint != "":
                sample_name = f"{sample_name}_{timepoint}"

            wps_filepath = os.path.join(DATA_DIR, sample_name, f"wps_{sample_name}.csv")

            try:
                if reference_df is None:
                    reference_df = pl.read_csv(wps_filepath, columns=["abs_pos"])
                wps_df = pl.read_csv(wps_filepath, columns=["abs_pos", "spanning_frags", "ending_frags"])
                wps_df = wps_df.with_columns((pl.col("spanning_frags") / (pl.col("spanning_frags") + pl.col("ending_frags"))).alias(sample_name))
                wps_df = wps_df.select(["abs_pos", sample_name])
                dfs.append(wps_df)
            except Exception as e:
                print(f"exception {e} when processing {wps_filepath}")

        wpr_all_samples = reference_df
        for df in dfs:
            wpr_all_samples = wpr_all_samples.join(df, on="abs_pos", how="left")

        wpr_all_samples = wpr_all_samples.drop_nans()
        abs_pos: pl.Series = wpr_all_samples.drop_in_place("abs_pos")
        return (abs_pos, wpr_all_samples)

    if False:
        abs_pos, wpr_all_samples = load_wpr_all_samples()
        wpr_all_samples_abs_pos = wpr_all_samples.with_columns(
            pl.Series("abs_pos", abs_pos)
        )
        wpr_all_samples_abs_pos.write_csv("./cache/wpr_all_samples.csv")
    else:
        wpr_all_samples = pl.read_csv("./cache/wpr_all_samples.csv")
        abs_pos = wpr_all_samples.drop_in_place("abs_pos")
    return (
        abs_pos,
        load_wpr_all_samples,
        wpr_all_samples,
        wpr_all_samples_abs_pos,
    )


@app.cell(hide_code=True)
def _(wpr_all_samples):
    def get_timepoint(n: str, hom=False) -> str:
        s = n.split("_")
        if len(s) == 1:
            return "ctrl"
        elif len(s) == 2:
            tp = s[1]
            if hom:
                match tp:
                    case "BL" | "BLrr": tp = "BL"
                    case "W2" | "W3" | "W4" | "M1": tp = "c1"
                    case "W6": tp = "c2"
                    case "M3" | "M6": tp = "end"
            return tp
        else:
            raise ValueError("weird split")


    timepoint_annos: list[str] = list(map(lambda n: get_timepoint(n, hom=True), wpr_all_samples.columns))
    return get_timepoint, timepoint_annos


@app.cell
def _(PCA, StandardScaler, pl, timepoint_annos, wpr_all_samples):
    def do_pca(wpr_df, timepoint_annos):
        samplename_annos = wpr_df.columns
    
        scaler = StandardScaler(with_mean=True, with_std=True)
        wpr_df = wpr_df.transpose()
        wpr_array_norm = scaler.fit_transform(wpr_df.to_numpy().T).T
        wpr_df_norm = pl.DataFrame(wpr_array_norm, schema=wpr_df.schema)

        pca = PCA(n_components=3)
        wpr_pca = pca.fit_transform(wpr_df_norm)
    
        pca_df = pl.DataFrame(wpr_pca, schema=["PC1", "PC2", "PC3"])
        pca_df = pca_df.with_columns(pl.Series("Patient", samplename_annos),
                                     pl.Series("Timepoint", timepoint_annos))
        return pca_df, pca
    
    pca_df, pca = do_pca(wpr_all_samples, timepoint_annos)
    return do_pca, pca, pca_df


@app.cell
def _(pca, pca_df, plt, sns):
    plt.figure(figsize=(22, 6))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Timepoint", palette="tab10", s=100, alpha=0.5)
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)")
    plt.title("PCA of Genomic Sites WPS Across Patients")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.show()
    return


@app.cell
def _(pca_df, sns):
    sns.kdeplot(data=pca_df, x="PC1", hue="Timepoint", palette="tab10", alpha=0.3, fill=True)
    return


@app.cell
def _(pca_df, px):
    fig = px.scatter_3d(pca_df, x="PC1", y="PC2", z="PC3",  color="Timepoint", color_continuous_scale="viridis",
                        symbol="Timepoint", opacity=0.8)
    fig.update_layout(scene=dict(xaxis_title="PC1", yaxis_title="PC2", zaxis_title="PC3"))
    return (fig,)


@app.cell
def _(pca, sns):
    # @NOTE(ds): Here, we inspect the loadings for principle component 1. They
    # represent how strongly the respective original feature contributes to
    # that particular component. Absolute value represents importance, sign
    # represents direction.
    component_num = 3
    sns.histplot(pca.components_[component_num-1, :])
    return (component_num,)


@app.cell
def _(component_num, np, pca, pd, timepoint_annos, wpr_all_samples):
    mean_df = wpr_all_samples.mean().to_numpy().flatten()
    min_indices = np.argsort(pca.components_[component_num-1])[0:100]
    min_loadings = wpr_all_samples[min_indices, :].mean().to_numpy().flatten()
    max_indices = np.argsort(pca.components_[component_num-1])[::-1][0:100]
    max_loadings = wpr_all_samples[max_indices, :].mean().to_numpy().flatten()
    df = pd.DataFrame({"timepoints": timepoint_annos,
                       "means": mean_df,
                       "min_loadings": min_loadings,
                       "max_loadings": max_loadings})
    return df, max_indices, max_loadings, mean_df, min_indices, min_loadings


@app.cell
def _(df, sns):
    sns.kdeplot(df, x="max_loadings", hue="timepoints")
    return


if __name__ == "__main__":
    app.run()
