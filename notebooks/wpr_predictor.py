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
    from sklearn.model_selection import train_test_split, cross_val_score, cross_val_predict, StratifiedKFold
    from sklearn.linear_model import LogisticRegressionCV
    from sklearn.manifold import MDS
    from sklearn.metrics import accuracy_score, roc_auc_score, pairwise_distances, \
                                classification_report, confusion_matrix, ConfusionMatrixDisplay, \
                                roc_curve, roc_auc_score
    from sklearn.preprocessing import StandardScaler
    return (
        ConfusionMatrixDisplay,
        LogisticRegressionCV,
        MDS,
        PCA,
        StandardScaler,
        StratifiedKFold,
        accuracy_score,
        classification_report,
        confusion_matrix,
        cross_val_predict,
        cross_val_score,
        mo,
        np,
        os,
        pairwise_distances,
        pd,
        pl,
        plt,
        px,
        roc_auc_score,
        roc_curve,
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


@app.cell(hide_code=True)
def _(clin_info_sheet, timepoint_annos, wpr_all_samples):
    # @NOTE(ds): Coding is as follows - relapse 0=no, 1=yes, survival 0=alive, 1=dead
    metadata_df = clin_info_sheet.loc[:, ["study_ID", "relapse", "survival"]]
    metadata_parallel = []
    for sample, timepoint in zip(wpr_all_samples.columns, timepoint_annos):
        if timepoint == "ctrl":
            metadata_parallel.append((timepoint, 0, 0))
        else:
            base = sample.split("_")[0]
            row = metadata_df[metadata_df["study_ID"] == base]
            metadata_parallel.append((timepoint, row["relapse"].iloc[0], row["survival"].iloc[0]))
    return base, metadata_df, metadata_parallel, row, sample, timepoint


@app.cell(hide_code=True)
def _(
    PCA,
    StandardScaler,
    metadata_parallel,
    pl,
    timepoint_annos,
    wpr_all_samples,
):
    def do_pca(wpr_df, timepoint_annos):
        samplename_annos = wpr_df.columns

        scaler = StandardScaler(with_mean=True, with_std=True)
        wpr_df = wpr_df.transpose()

        # @NOTE(ds): Double transpose because we want to scale across a patient, not across a site.
        wpr_array_norm = scaler.fit_transform(wpr_df.to_numpy().T).T
        # wpr_array_norm = scaler.fit_transform(wpr_df.to_numpy())
        wpr_df_norm = pl.DataFrame(wpr_array_norm, schema=wpr_df.schema)

        pca = PCA(n_components=3)
        wpr_pca = pca.fit_transform(wpr_df_norm)

        pca_df = pl.DataFrame(wpr_pca, schema=["PC1", "PC2", "PC3"])
        pca_df = pca_df.with_columns(pl.Series("Patient", samplename_annos),
                                     pl.Series("Timepoint", timepoint_annos),
                                     pl.Series("Survival", list(map(lambda x: "dead" if x[2] == 1 else "alive", metadata_parallel))),
                                     pl.Series("Relapse", list(map(lambda x: "yes" if x[1] == 1 else "no", metadata_parallel))))
        return pca_df, pca

    pca_df, pca = do_pca(wpr_all_samples, timepoint_annos)
    return do_pca, pca, pca_df


@app.cell(hide_code=True)
def _(pca, pca_df, plt, sns):
    plt.figure(figsize=(22, 6))
    # sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Timepoint", style="Timepoint", palette="tab10", s=140, alpha=0.5)
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Timepoint", palette="tab10", s=180, alpha=0.5)
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)")
    plt.title("PCA of Genomic Sites WPS Across Patients")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.show()
    return


@app.cell(hide_code=True)
def _(pca_df, sns):
    sns.kdeplot(data=pca_df, x="PC1", hue="Timepoint", palette="tab10", alpha=0.3, fill=True)
    return


@app.cell(hide_code=True)
def _(pca_df, px):
    fig = px.scatter_3d(pca_df, x="PC1", y="PC2", z="PC3", color="Timepoint", color_continuous_scale="viridis",
                        symbol="Survival", opacity=0.8)
    fig.update_layout(scene=dict(xaxis_title="PC1", yaxis_title="PC2", zaxis_title="PC3"))
    return (fig,)


@app.cell(hide_code=True)
def _(pca, sns):
    # @NOTE(ds): Here, we inspect the loadings for principle component 1. They
    # represent how strongly the respective original feature contributes to
    # that particular component. Absolute value represents importance, sign
    # represents direction.
    component_num = 3
    sns.histplot(pca.components_[component_num-1, :])
    return (component_num,)


@app.cell(hide_code=True)
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


@app.cell(hide_code=True)
def _(df, sns):
    sns.kdeplot(df, x="max_loadings", hue="timepoints")
    return


@app.cell(hide_code=True)
def _(
    MDS,
    StandardScaler,
    metadata_parallel,
    pairwise_distances,
    pl,
    timepoint_annos,
    wpr_all_samples,
):
    def do_mds(wpr_df, timepoint_annos) -> pl.DataFrame():
        samplename_annos = wpr_df.columns
        scaler = StandardScaler(with_mean=True, with_std=True)

        # @NOTE(ds): Input matrix must again be patients x features.
        wpr_df = wpr_df.transpose()

        # @NOTE(ds): Double transpose because we want to scale across a patient, not across a site.
        # > wpr_array_norm = scaler.fit_transform(wpr_df.to_numpy())
        wpr_array_norm = scaler.fit_transform(wpr_df.to_numpy().T).T
        wpr_df_norm = pl.DataFrame(wpr_array_norm, schema=wpr_df.schema)
        dist_matrix = pairwise_distances(wpr_df_norm, metric="euclidean")  # correlation, euclidean, cosine

        mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42, metric=True)
        mds_fit = mds.fit_transform(dist_matrix)

        mds_df = pl.DataFrame(mds_fit, schema=["MDS1", "MDS2"])
        mds_df = mds_df.with_columns(pl.Series("Patient", samplename_annos),
                                     pl.Series("Timepoint", timepoint_annos),
                                     pl.Series("Survival", list(map(lambda x: "dead" if x[2] == 1 else "alive", metadata_parallel))),
                                     pl.Series("Relapse", list(map(lambda x: "yes" if x[1] == 1 else "no", metadata_parallel))))
        return mds_df


    mds_df = do_mds(wpr_all_samples, timepoint_annos)
    return do_mds, mds_df


@app.cell(hide_code=True)
def _(mds_df, sns):
    sns.scatterplot(x="MDS1", y="MDS2", hue="Timepoint", style="Survival", data=mds_df)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Regression
        Can we make any meaningful predictions using the windowed protection ratio?
        """
    )
    return


@app.cell
def _(
    LogisticRegressionCV,
    StandardScaler,
    StratifiedKFold,
    classification_report,
    cross_val_predict,
    pl,
):
    def do_regression(wpr_df, y):
        scaler = StandardScaler(with_mean=True, with_std=True)
        wpr_df = wpr_df.transpose()

        # @NOTE(ds): Double transpose because we want to scale across a patient, not across a site.
        # > wpr_array_norm = scaler.fit_transform(wpr_df.to_numpy())
        wpr_array_norm = scaler.fit_transform(wpr_df.to_numpy().T).T
        wpr_df_norm = pl.DataFrame(wpr_array_norm, schema=wpr_df.schema)

        cv = StratifiedKFold(10, shuffle=True)  # Better than plain KFold for class imbalance
        clf = LogisticRegressionCV(
            Cs=10,  # Grid of inverse regularization strengths
            cv=cv,
            penalty="elasticnet",
            l1_ratios=[0.1, 0.5, 0.9],
            solver="saga",
            scoring="roc_auc",  # "accuracy" didn't perform so well
            max_iter=5000,
            tol=1e5,
            verbose=10,
            n_jobs=14,  # up to 15 should be okay
            refit=True
        )

        print("fitting model...")
        clf.fit(wpr_df_norm, y)

        print("predicting based on model...")
        y_pred = clf.predict(wpr_df_norm)
        print(classification_report(y, y_pred))

        # Get cross-validated predictions of probability for class 1.
        y_proba = cross_val_predict(
            clf, wpr_df_norm, y, cv=cv,
            method="predict_proba",
            n_jobs=12
        )[:, 1]

        return clf, y_pred, y_proba
    return (do_regression,)


@app.cell
def _(do_regression, metadata_parallel, wpr_all_samples):
    # @NOTE(ds): tuple ordering is: (timepoint, relapse, survival)
    #
    # y = list(map(lambda x: "control" if x[0].startswith("ctrl") else "tumor",
    #              metadata_parallel))
    #
    y = list(map(lambda x: "relapse" if x[1] == 1 else "no relapse",
                 metadata_parallel))
    clf_fit, y_pred, y_proba = do_regression(wpr_all_samples, y)
    return clf_fit, y, y_pred, y_proba


@app.cell
def _(ConfusionMatrixDisplay, clf_fit, confusion_matrix, plt, y, y_pred):
    cm = confusion_matrix(y, y_pred)
    disp = ConfusionMatrixDisplay(cm, display_labels=clf_fit.classes_)
    disp.plot()
    plt.show()
    return cm, disp


@app.cell
def _(plt, roc_auc_score, roc_curve, y, y_proba):
    fpr, tpr, thresholds = roc_curve(y, y_proba)
    auc = roc_auc_score(y, y_proba)
    plt.plot(fpr, tpr, label=f'ROC curve (AUC = {auc:.2f})', color='darkorange', lw=2)
    plt.plot([0, 1], [0, 1], color='navy', linestyle='--', lw=2)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC)')
    plt.legend(loc='lower right')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.show()
    return auc, fpr, thresholds, tpr


@app.cell
def _(abs_pos, clf_fit, np):
    # @NOTE(ds): Inspect feature importance.
    coefs = clf_fit.coef_
    top_features = np.argsort(np.abs(coefs), axis=1)[:, -10:]  # top 10 per class
    print(abs_pos[top_features.flatten()])
    return coefs, top_features


if __name__ == "__main__":
    app.run()
