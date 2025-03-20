import requests

gene_list: list[str] = [
 "ACTB", "ADAMTS12", "ANK3", "APC2", "ARHGEF1", "ARID1A", "ARID1B",
 "ARID5B", "ATM", "ATP6V1B2", "ATRX", "B2M",
 "BAHCC1",  # manually curated, located at chr17:79,373,540
 "BCL10",
 "BCL11A", "BCL2", "BCL2-NA", "BCL6", "BCL6-NA", "BCL7A", "BCOR",
 "BIRC6", "BRAF", "BRCA1", "BTG1", "BTG2", "BTK", "CACNA1E",
 "CADPS", "CARD11", "CCL4", "CCND3", "CD19", "CD274", "CD28",
 "CD58", "CD70", "CD79A", "CD79B", "CD83", "CDC73", "CDH23",
 "CDKN1A", "CDKN2A", "CDKN2B", "CHST2", "CIITA", "CIITA,CIITA-TRANSLOCATIONS",
 "CIITA-TRANSLOCATIONS", "COQ7", "CPEB2", "CPO", "CREBBP", "CRIP1", "CTNNB1",
 "CUX1", "CXCR4", "CXCR5", "CT55", "CYLD", "DAZAP1", "DDX3X", "DNAH7",
 "DNMT3A", "DOCK2", "DPYS", "DTX1", "DUSP2", "EBF1", "EEF1A1",
 "EIF2AK3", "EP300", "ERCC5", "ETS1", "ETV6", "EZH2", "PIEZO1",  # "FAM38A",
 "FAM46C", "FAS", "FAT3", "FAT4", "FCGR2B", "FOXC1", "FOXO1",
 "FOXO3", "FUT5", "FYN", "GNA13", "GNAI2", "GPR98", "GRB2",
 "GRHPR", "GSG2", "GTF2I", "HDAC7", "MROH2B",  # "HEATR7B2",
 "HIST1H1B", "HIST1H1C",
 "HIST1H1D", "HIST1H1E", "HIST1H2AC", "HIST1H2AM", "HIST1H2BC", "HIST1H2BK",
 "HIST1H3B", "HIST2H2AA3", "HIST2H2BE", "HIVEP1", "HLA-A", "HLA-B", "HLA-C",
 "HLA-DMA", "HLA-DMB", "HNF1B", "HNRNPD", "HVCN1", "ID3", "IDH2", "IGH@-BCL6",
 "IGH@-MYC", "IGHJ1", "IGHJ2", "IGHJ3",
 "IGHJ4", "IGHJ5", "IGHJ6",  # "IGHJ6,IGHJ4,IGHJ5,IGH@-MYC",
 "IGHV1-3", "IGHV1-69", "IGHV2-5", "IGHV2-70", "IGHV3-23", "IGHV3-33",
 "IGHV3-7", "IGHV4-34", "IGHV4-39", "IGHV5-51", "IGLL5", "IHGV4-34-INTERGENIC",
 "IKZF3", "IL10RA", "IL16", "IL6", "IRF2BP2", "IRF4", "IRF8", "ITPKB", "JAK1",
 "JAK2", "KDSR", "KLHL14", "KLHL21", "KLHL42", "KLHL6", "KMT2A", "KMT2B",
 "KMT2D", "KRAS", "LAMA2", "LRP1B", "LRP2", "LTB", "LYN", "MALT1",
 "MAP2K1", "MCL1", "MDM2", "MEF2B", "MEF2C", "MFHAS1", "MIR17HG",
 "KMT2C",  # "MLL3",
 "MLLT3", "MPEG1", "MTOR", "MUC16", "MUC22",
 "MYC", "MYC-IGH", "MYD88", "NANOG", "NAV1", "NCKIPSD", "NFKBIA",
 "NFKBIE", "NFKBIZ", "NLRP8", "NOL9", "NOTCH1", "NOTCH2", "NRAS",
 "OSBPL10", "P2RY8", "PARP2", "PAX5", "PCDH10", "PCDH15", "PCLO",
 "PCNXL2",  # "PCNX2,PCNXL2"
 "PDCD1LG2", "PDE4DIP", "PDPK1", "PIK3CA", "PIK3R1",
 "PIM1", "PIM2", "PKD1", "PKHD1", "PLCG1", "POLG", "POT1",
 "POU2AF1", "POU2F2", "PPP1R9B", "PRDM1", "PRKCB", "PRPS1", "PTEN",
 "PTPN23", "PTPN6", "RAD9A", "RB1", "REL", "RELN", "RERE",
 "RHOA", "RRAGC", "S1PR2", "SEL1L3", "SESN1", "SETD1B", "SETD2",
 "SETX", "SF3B1", "SGK1", "SIN3A", "SMEK1", "SOCS1", "SPEN",
 "SPIB", "STAT3", "STAT5B", "STAT6", "SUSD2", "TACC2", "TBL1XR1",
 "TCF3", "TCHH", "TENM4", "TET2", "TFEB", "TFPT", "TLR2",
 "TMEM132B", "TMEM30A", "TMSB4X", "TNFAIP3", "TNFRSF14", "TOX", "TP53",
 "TP53BP1", "TP63", "TP73", "TRIM33", "TRIP12", "TRRAP", "TTN",
 "UBE2A", "UBE2O", "UGGT2", "UNC80", "VAV1", "VMP1", "VPS13B",
 "WDFY3", "WEE1", "XPO1", "YY1", "ZAN", "ZC3H12A", "ZEB2",
 "ZFP36L1", "ZNF292", "ZNF423", "ZNF804A"
]


def get_canonical_tss(gene_symbol: str) -> dict[str, object]:
    url: str = "https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/" \
              f"{gene_symbol}?expand=1"
    headers: dict[str, str] = {"Content-Type": "application/json"}
    response = requests.get(url, headers=headers)
    data: object = response.json()

    transcript: dict[str, object]
    for transcript in data.get("Transcript", []):  # type: ignore
        if transcript.get("is_canonical"):
            return {
                "gene": gene_symbol,
                "start": transcript["start"],
                "end": transcript["end"],
                "canonical_transcript": transcript["id"],
                "contig": transcript["seq_region_name"],
                "strand": transcript["strand"],
                "TSS": transcript["start"]
                if transcript["strand"] == 1 else transcript["end"]
            }
    raise ValueError(f"no canonical transcript found for {gene_symbol}")


with open("loci_covered_by_panel_canonical_tss.bed", "w") as bed_file:
    for gene in gene_list:
        try:
            d: dict[str, object] = get_canonical_tss(gene)
            tss: int = int(d["TSS"])  # type: ignore
            bed_file.write(
                f"{d['contig']}\t{tss-1000}\t{tss+1000}\t{d['gene']}\n"
            )
        except ValueError as e:
            print(e)
