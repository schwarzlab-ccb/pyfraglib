import marimo

__generated_with = "0.11.13"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _():
    import glob

    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    from math import isnan
    from cyvcf2 import VCF  # type: ignore
    return VCF, glob, isnan, pd, plt, sns


@app.cell(hide_code=True)
def _(VCF, boolean, pd):
    def load_vcf_files(
        vcf_files: list[str], filter_vaf: float, top_n: int | None = None,
        collapse_muts: bool = False, blacklist_genes=[]
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        print(f"loading VCF files, blacklisted: {blacklist_genes}")
        mutation_data: list[dict[str, object]] = []
        for vcf_file in vcf_files:
            vcf: VCF = VCF(vcf_file)  # type: ignore
            sample_name: str = vcf_file.split(sep='/')[-1].split(sep='.vcf')[0]
            this_mutation_data: list[dict[str, object]] = []

            for record in vcf:  # type: ignore
                gene: str = record.INFO.get("Gene.refGene", f"{record.CHROM}:{record.POS}")  # type: ignore
                mutation_type: str = record.INFO.get("Func.refGene", "Unknown")  # type: ignore
                vaf: float = float(record.INFO.get("AF_rel", None))  # type: ignore

                if not vaf or vaf*100 < filter_vaf:
                    continue
                if gene in blacklist_genes:
                    continue

                new_mutation_entry: dict[str, object] = {
                    "Sample": sample_name,
                    "Gene_Name": gene,
                    "Position": record.POS,
                    "Gene_Name_Position": f"{gene} ({record.POS})",  # type: ignore
                    "Mutation_Type": mutation_type,
                    "VAF": vaf
                }

                # @NOTE(ds): Nasty. If we want to keep a max. of 1 mutation per gene, we must iterate the
                # mutation entries and replace the entry for the current gene if the current gene's VAF is
                # higher (we want to keep the genes with the highest VAFs).
                if collapse_muts:
                    found: boolean = False
                    for idx, entry in enumerate(this_mutation_data):
                        if entry["Gene_Name"] == new_mutation_entry["Gene_Name"]:
                            found = True
                            if entry["VAF"] < new_mutation_entry["VAF"]:
                                this_mutation_data[idx] = new_mutation_entry
                    if not found:
                        this_mutation_data.append(new_mutation_entry)
                else:
                    this_mutation_data.append(new_mutation_entry)

            this_top_n: int = len(this_mutation_data) if top_n is None else min(top_n, len(this_mutation_data))
            if this_top_n == 0:
                print(f"no mutations found/passing filter for sample {sample_name}")
                continue

            this_mutation_data.sort(key=lambda record: record["VAF"], reverse=True)  # type: ignore
            this_mutation_data = this_mutation_data[:this_top_n]
            print(f"kept {len(this_mutation_data)} mutations for sample {sample_name}")

            mutation_data += this_mutation_data

        df: pd.DataFrame = pd.DataFrame(mutation_data)
        mutation_matrix: pd.DataFrame = df.pivot_table(
            index="Sample", 
            columns="Gene_Name" if collapse_muts else "Gene_Name_Position",
            values="Mutation_Type", 
            aggfunc=lambda x: ", ".join(x) if isinstance(x, list) else x  # type: ignore
        )
        vaf_matrix: pd.DataFrame = df.pivot_table(
            index="Sample", 
            columns="Gene_Name" if collapse_muts else "Gene_Name_Position",
            values="VAF"
        )
        return mutation_matrix, vaf_matrix
    return (load_vcf_files,)


@app.cell(hide_code=True)
def _(isnan, pd, plt, sns):
    def plot_oncoplot(
        mutation_matrix: pd.DataFrame, vaf_matrix: pd.DataFrame, title: str,
        figsize: tuple[int, int] = (25, 25)
    ) -> None:
        print(f"creating Oncoplot of {mutation_matrix.shape} matrix")

        mutation_matrix.fillna("No Mutation", inplace=True)
        unique_mut_types = pd.unique(mutation_matrix.values.ravel())  # type: ignore
        integer_coding: dict[str, int] = {mut_type: i for i, mut_type in enumerate(unique_mut_types)}  # type: ignore
        num_unique_muts = len(integer_coding)
        int_matrix = mutation_matrix.map(lambda mut_type: integer_coding[mut_type])

        plt.figure(figsize=figsize)
        cmap = sns.color_palette("deep", num_unique_muts)
        cmap[integer_coding["No Mutation"]] = (1.0, 1.0, 1.0)
        ax = sns.heatmap(
            int_matrix,
            annot=vaf_matrix.map(lambda x: "" if isnan(x) else f"{round(x*100, 1)}"),  # type: ignore
            annot_kws={"fontsize": 12, "color": "black", "alpha": 0.9},
            fmt="s",
            linewidths=0.5,
            linecolor="gray",
            cbar=True,
            cmap=cmap  # type: ignore
        )

        colorbar = ax.collections[0].colorbar  # type: ignore
        r = colorbar.vmax - colorbar.vmin  # type: ignore
        colorbar.set_ticks([colorbar.vmin + r / num_unique_muts * (0.5 + i)  # type: ignore
                            for i in range(num_unique_muts)])
        colorbar.set_ticklabels(list(integer_coding.keys()))  # type: ignore

        plt.title(title)
        plt.xlabel("Genes")
        plt.ylabel("Samples")
        plt.xticks(rotation=90)
        plt.show()  # type: ignore
    return (plot_oncoplot,)


@app.cell
def _(glob, load_vcf_files):
    vcf_files: list[str] = glob.glob("/mnt/ramses/projects/uk-lymphoma-cfdna/PCNSL/mutations/*.vcf")
    filter_vaf: float = 1.0
    top_n: int = 3
    mutation_matrix, vaf_matrix = load_vcf_files(vcf_files[:], filter_vaf=filter_vaf, top_n=top_n,
                                                 collapse_muts=True, blacklist_genes=["ACTB", "MUC16", "MUC22", "HLA-A", "HLA-B", "HLA-C", "LINC00226:LINC00221"])
    return filter_vaf, mutation_matrix, top_n, vaf_matrix, vcf_files


@app.cell
def _(filter_vaf, mutation_matrix, plot_oncoplot, top_n, vaf_matrix):
    plot_oncoplot(mutation_matrix, vaf_matrix, figsize=(28, 22),
                  title=f"Oncoplot - Top {top_n} mutations per samples (only VAFs > {filter_vaf}%)")
    return


if __name__ == "__main__":
    app.run()
