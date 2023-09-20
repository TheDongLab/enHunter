### modified from Ruifeng (/data/neurogen/AMPPD/Aim1_DE/src/eRNA directory)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import ttest_ind

from matplotlib.backends.backend_pdf import PdfPages

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# %%
# gene_name_ls = ["chr10_45472340_45475290_plus","chrX_120452800_120454340_minus"]
gene_name_ls =["chr16_11639850_11640300_minus_chr16_11640470_11641120_plus"]

#gene_names = ["chr16_11611980_11612400","chr16_11612780_11613560","chr16_11613470_11613780","chr16_11613950_11614560"]


# %%
df =  pd.read_csv("../inputs/bidirectional_pairs/eRNA.bidirectional_pairs.class1_2_TNE_900.xls",index_col=0,header=0,sep="\t")
df = df.loc[gene_name_ls,:].T

######
PPMI_cov = pd.read_csv("../run_inout/PPMI_CaseA_CtrlA_cov.tsv", sep="\t", index_col=0, header=0)
PPMI_data = pd.merge(df,PPMI_cov.loc[:,["case_control_other_latest"]],left_index=True,right_index=True)
PPMI_data.loc[:,"dataset"]="PPMI"


PDBF_cov = pd.read_csv("../run_inout/PDBF_CaseA_CtrlA_cov.tsv", sep="\t", index_col=0, header=0)
PDBF_data = pd.merge(df,PDBF_cov.loc[:,["case_control_other_latest"]],left_index=True,right_index=True)
PDBF_data.loc[:,"dataset"]="PDBF"

df = pd.concat([PPMI_data,PDBF_data],axis=0)

# %%
palette ={"Case": "#be1e2d", "Control": "#bcbec0"}
with PdfPages("../output/expression_boxplot_2Combined_TNEs.pdf") as expr_pdf:
    ## loop each single gene
    for gene in gene_name_ls:
        print(gene)
        ppmi_case = df.loc[(df["dataset"]=="PPMI") &(df["case_control_other_latest"]=="Case"),gene].values
        ppmi_ctrl = df.loc[(df["dataset"]=="PPMI") &(df["case_control_other_latest"]=="Control"),gene].values
        t, ppmi_p = ttest_ind(ppmi_case,ppmi_ctrl)

        pdbf_case = df.loc[(df["dataset"]=="PDBF") &(df["case_control_other_latest"]=="Case"),gene].values
        pdbf_ctrl = df.loc[(df["dataset"]=="PDBF") &(df["case_control_other_latest"]=="Control"),gene].values
        t, pdbf_p = ttest_ind(pdbf_case,pdbf_ctrl)

        plt.figure(figsize=(6, 8), constrained_layout=True)
        # sns.boxplot(x="dataset", y=gene, hue="case_control_other_latest",color="g",hue_order=["Case","Control"], data=df, palette=palette,flierprops={"marker": "o"})
        sns.boxplot(x="dataset", y=gene, hue="case_control_other_latest",color="#737373",hue_order=["Case","Control"], data=df, palette=palette)
        sns.despine(offset=10, trim=True)
        plt.ylabel("Gene expression (read counts)",fontname="Arial", fontsize=14)
        plt.xlabel('Dataset', fontsize=14,fontname="Arial")
        plt.xticks(fontsize=12, rotation=0, ha='center',fontname="Arial")
        plt.yticks(fontsize=12,fontname="Arial")
        title_str = gene + "\nPPMI: "+str(ppmi_p)
        title_str += "\nPDBF: "+str(pdbf_p)
        plt.title(title_str,fontname="Arial")
        plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left",frameon=False)
        expr_pdf.savefig()
        plt.close()

