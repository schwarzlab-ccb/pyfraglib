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
    import seaborn as sns

    from sklearn.decomposition import PCA
    from sklearn.model_selection import train_test_split, cross_val_score
    from sklearn.linear_model import LogisticRegressionCV
    from sklearn.metrics import accuracy_score, roc_auc_score
    return (
        LogisticRegressionCV,
        PCA,
        accuracy_score,
        cross_val_score,
        mo,
        np,
        os,
        pd,
        plt,
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
def _(DATA_DIR, os, pd, sample_sheet):
    def load_wps_all_BL_samples() -> pd.DataFrame:
        data: pd.DataFrame | None = None
        index: pd.Series | None = None

        for id, timepoint, timepoint_hom in zip(sample_sheet["sample_id"], sample_sheet["time_point"], sample_sheet["time_point_hom"]):
            if timepoint_hom != "BL":
                continue

            sample_name = f"{id}_{timepoint}"
            wps_filepath = os.path.join(DATA_DIR, sample_name, f"wps_{sample_name}.csv")
            wps_df = pd.read_csv(wps_filepath, low_memory=False)
            if data is None:
                data = pd.DataFrame(wps_df["wps"])
                data.columns = [sample_name]
                index = wps_df["abs_pos"]
            else:
                col_names = list(data.columns)
                col_names.append(sample_name)
                data = pd.concat([data, wps_df["wps"]], axis=1)
                data.columns = col_names

        # @NOTE(ds): We only have a few uninformative features (~0.08%)
        data.index = index
        wps_all_samples = data.loc[(data != 0).any(axis=1)]
        return wps_all_samples

    wps_all_samples = load_wps_all_BL_samples()
    return load_wps_all_BL_samples, wps_all_samples


@app.cell(hide_code=True)
def _(clin_info_sheet, wps_all_samples):
    def transform(s: str) -> str:
        t = s.split("_")[0]
        return t

    metadata_df = clin_info_sheet.loc[:, ["study_ID", "relapse"]]  # survival
    wps_study_ids = wps_all_samples.columns.map(lambda s: transform(s))
    wps_all_samples.columns = wps_study_ids
    metadata_df = metadata_df.set_index("study_ID").loc[wps_study_ids]
    return metadata_df, transform, wps_study_ids


@app.cell(hide_code=True)
def _(metadata_df, np, train_test_split, wps_all_samples):
    # z-score normalization
    mean = wps_all_samples.mean(axis=0).values[np.newaxis, :]
    std = wps_all_samples.std(axis=0).values[np.newaxis, :]
    wps_all_samples_norm = (wps_all_samples - mean) / std

    # train-test split
    patient_ids = metadata_df.index
    train_patients, test_patients = train_test_split(patient_ids, stratify=metadata_df["relapse"],
                                                     test_size=0.2, random_state=42)
    X_train = wps_all_samples_norm[train_patients]
    X_test = wps_all_samples_norm[test_patients]
    y_train = metadata_df.T[train_patients]
    y_test = metadata_df.T[test_patients]
    return (
        X_test,
        X_train,
        mean,
        patient_ids,
        std,
        test_patients,
        train_patients,
        wps_all_samples_norm,
        y_test,
        y_train,
    )


@app.cell(hide_code=True)
def _(plt, wps_all_samples, wps_all_samples_norm):
    _plt, _axis = plt.subplots(2, 1)
    _axis[0].plot(wps_all_samples_norm["DED005"])
    _axis[1].plot(wps_all_samples["DED005"])
    return


@app.cell
def _(LogisticRegressionCV, X_train, y_train):
    lasso = LogisticRegressionCV(
        penalty="l1",
        solver="saga",
        cv=5,
        scoring="roc_auc",
        max_iter=6000,
        tol=3e-4,  # 1e-4
        verbose=3,
        n_jobs=-1,
    )
    lasso.fit(X_train.T, y_train.iloc[0].astype(int))
    return (lasso,)


@app.cell
def _(
    X_test,
    accuracy_score,
    lasso,
    np,
    roc_auc_score,
    wps_all_samples,
    y_test,
):
    y_pred = lasso.predict(X_test.T)
    y_prob = lasso.predict_proba(X_test.T)[:, 1]  # Probability of class 1

    accuracy = accuracy_score(y_test.iloc[0].astype(int), y_pred)
    auc = roc_auc_score(y_test.iloc[0].astype(int), y_prob)
    print(f"✅ Accuracy: {accuracy:.3f}")
    print(f"✅ ROC-AUC: {auc:.3f}")
    feature_importance = np.abs(lasso.coef_).flatten()
    top_sites = np.argsort(feature_importance)[::-1][:10]  # Top 10 sites

    print("\n🔬 Top 10 Most Important Genomic Sites:")
    for idx in top_sites:
        site_name = wps_all_samples.index[idx]
        print(f"{site_name}: Importance {feature_importance[idx]:.4f}")
    return (
        accuracy,
        auc,
        feature_importance,
        idx,
        site_name,
        top_sites,
        y_pred,
        y_prob,
    )


@app.cell
def _(PCA, pd, wps_all_samples_norm):
    pca = PCA(n_components=2)
    wps_all_samples_norm_pca = pca.fit_transform(wps_all_samples_norm.T)
    pca_df = pd.DataFrame(wps_all_samples_norm_pca, columns=["PC1", "PC2"])
    pca_df["Patient"] = wps_all_samples_norm.columns
    return pca, pca_df, wps_all_samples_norm_pca


@app.cell
def _(pca, pca_df, plt, sns):
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Patient", palette="tab10", s=100)
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)")
    plt.title("PCA of Genomic Sites WPS Across Patients")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.show()
    return


if __name__ == "__main__":
    app.run()
