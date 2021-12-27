#!/usr/bin/env python3
import matplotlib as mpl
mpl.use("Qt5Agg")
import matplotlib.pyplot as plt
import gspread
import gspread_dataframe as gsdf

import os
import pandas as pd
import numpy as np
import glob
import re
import scipy.optimize as spo
import scipy.stats as sps
import seaborn as sns

# from datetime import datetime
from funcs.stritfuncs import get_q_dir, get_wwtp_df, sample_name_parser, excel_cols

pd.set_option("display.width", 1000)
pd.options.display.max_rows = 2000
pd.options.display.max_columns = 15

gc = gspread.service_account(filename="cowwid-wave5-e31fa2771e32.json")

##------------------------------------------------------------------------------------------------##
def _parse_assay_counts(droplet_cts):
    """
    Data parser for SARS-CoV-2 Variant ddPCR Drop-off assay.

    takes a pandas Series or tuple of values with the following order:
    (
        TotalDroplets,
        FAM single-positive droplets,
        HEX single-positive droplets,
        Cy5 single-positive droplets,
        Double-positive droplets (either FAM+HEX or Cy5+HEX)
    )

    Right now only works with the S:L452R (Delta) & 69-70del (Omicron) assays
    Ask Lea Caduff for more details
    """
    if type(droplet_cts) == pd.Series:
        tot_ct = droplet_cts.iloc[0] #TotalDroplets
        fam_ct = droplet_cts.iloc[1] #FAM single-positive droplets
        hex_ct = droplet_cts.iloc[2] #HEX single-positive droplets
        cy5_ct = droplet_cts.iloc[3] #Cy5 single-positive droplets
        dbl_ct = droplet_cts.iloc[4] #Double-positive droplets (either FAM+HEX or Cy5+HEX)
    elif type(droplet_cts) == tuple:
        tot_ct = droplet_cts[0] #TotalDroplets
        fam_ct = droplet_cts[1] #FAM single-positive droplets
        hex_ct = droplet_cts[2] #HEX single-positive droplets
        cy5_ct = droplet_cts[3] #Cy5 single-positive droplets
        dbl_ct = droplet_cts[4] #Double-positive droplets (either FAM+HEX or Cy5+HEX)

    if np.isnan(cy5_ct):##DELTA Assay
        mutat_ct = dbl_ct #Mutation-positive droplets
        wldtp_ct = fam_ct #Mutation-negative, SARS-CoV-2 positive droplets ("wildtype")
        unkwn_ct = hex_ct #Unknown-yet-positive droplets; Right now unused
    elif np.isnan(fam_ct):##OMICRON Assay
        mutat_ct = hex_ct #Mutation-positive droplets
        wldtp_ct = dbl_ct #Mutation-negative, SARS-CoV-2 positive droplets ("wildtype")
        unkwn_ct = cy5_ct #Unknown-yet-positive droplets; Right now unused
    else:
        raise TypeError("That boy ain't right")

    return tot_ct, mutat_ct, wldtp_ct, unkwn_ct
##------------------------------------------------------------------------------------------------##
def calc_neg_droplets(droplet_cts):

    tot_ct, mutat_ct, wldtp_ct, unkwn_ct = _parse_assay_counts(droplet_cts)

    pos_ct = mutat_ct + wldtp_ct + unkwn_ct
    neg_ct = tot_ct - pos_ct

    return neg_ct
##------------------------------------------------------------------------------------------------##
def calc_perc_mutation(droplet_cts):

    tot_ct, mutat_ct, wldtp_ct, unkwn_ct = _parse_assay_counts(droplet_cts)

    pos_ct = mutat_ct + wldtp_ct + unkwn_ct
    perc_mut = mutat_ct / pos_ct*100

    return perc_mut
##------------------------------------------------------------------------------------------------##
def r_hat(droplet_cts):
    """
    Function translated from R, orig. author David Dreifuss:

    r_hat <- function(x){
      if(is.null(dim(x))){
        x_0 <- x[1]
        x_1 <- x[2]
        x_2 <- x[3]
      }else{
        x_0 <- x[,1]
        x_1 <- x[,2]
        x_2 <- x[,3]
      }
      n <- x_0 + x_1 + x_2
      l1 <- -log((x_0+x_1)/n)
      l2 <- -log(x_0/n)
      1-(l1/l2)
    }

    r̂ is the estimate of the ratio of λ1 and λ2
    """
    tot_ct, mutat_ct, wldtp_ct, unkwn_ct = _parse_assay_counts(droplet_cts)

    pos_ct = mutat_ct + wldtp_ct + unkwn_ct
    neg_ct = tot_ct - pos_ct

    lambda1 = -np.log((neg_ct + mutat_ct) / tot_ct)
    lambda2 = -np.log(neg_ct / tot_ct)

    return 1 - (lambda1 / lambda2)
##------------------------------------------------------------------------------------------------##
def loglik_trinom_prof(droplet_cts,ratio):
    """
    Function translated from R, orig. author David Dreifuss:

    loglik_trinom_prof <- function(y, ratio){
      if(is.null(dim(y))){
        x_0 <- y[1]
        x_1 <- y[2]
        x_2 <- y[3]
      }else{
        x_0 <- y[,1]
        x_1 <- y[,2]
        x_2 <- y[,3]
      }
      n <-  x_0 + x_1 + x_2
      x_0*log(x_0/n) +
        ifelse(x_1 == 0, yes = 0, no = x_1*log(-x_0/n + exp((1-ratio)*log(x_0/n)))) +
        ifelse(x_2 == 0, yes = 0, no = x_2*log(1 - exp((1-ratio)*log(x_0/n))))
    }
    """

    tot_ct, mutat_ct, wldtp_ct, unkwn_ct = _parse_assay_counts(droplet_cts)

    pos_ct = mutat_ct + wldtp_ct + unkwn_ct
    neg_ct = tot_ct - pos_ct

    neg_ratio = neg_ct / tot_ct
    lambda_hat = -np.log(neg_ratio)

    P0 = neg_ct * -lambda_hat

    if mutat_ct == 0:
        P1 = 0
    else:
        P1 = mutat_ct * np.log(-neg_ratio + np.exp((1-ratio) * -lambda_hat))

    if wldtp_ct == 0:
        P2 = 0
    else:
        P2 = wldtp_ct*np.log(1 - np.exp((1-ratio) * -lambda_hat))

    trinomial = P0 + P1 + P2

    if np.isinf(trinomial):
        print("trinomial infinite result")
        return 1
    else:
        # print(trinomial)
        return trinomial
##------------------------------------------------------------------------------------------------##
def lower_e(droplet_cts, level=0.95):
    """
    Function translated from R, orig. author David Dreifuss:

    lower_e <- function(x, level=0.95){
      if(is.nan(r_hat(x))){
        return(NaN)
      }else if(r_hat(x) == 0){
        return(0)
      }else{
        uniroot(function(y){2*(loglik_trinom_prof(x,y) - loglik_trinom_prof(x,r_hat(x))) + qchisq(level, 1)},
                interval=c(0, r_hat(x)))$root
      }
    }"""
    r_hat_val = r_hat(droplet_cts)
    if np.isnan(r_hat_val):
        print(np.nan)
        return np.nan
    elif r_hat_val == 0:
        print(0)
        return 0
    else:
        lp_rhat = loglik_trinom_prof(droplet_cts, r_hat_val)
        chi_sq = sps.chi2.ppf(level, 1)

        try:
            root = spo.brentq(lambda r: 2*(loglik_trinom_prof(droplet_cts,r) - lp_rhat) + chi_sq,
                a=0.00001,
                b=r_hat_val
            )
            return root
        except ValueError:
            return np.nan


##------------------------------------------------------------------------------------------------##
def upper_e(droplet_cts, level=0.95):
    """
    Function translated from R, orig. author David Dreifuss:

    upper_e <- function(x, level=0.95){
      if (is.nan(r_hat(x))){
        return(NaN)
      }else if (r_hat(x) == 1){
        return(1)
      }else{
        uniroot(function(y){2*(loglik_trinom_prof(x,y) - loglik_trinom_prof(x,r_hat(x))) + qchisq(level, 1)},
                interval=c(r_hat(x), 1), f.lower = qchisq(level, 1))$root
      }
    }
    """
    r_hat_val = r_hat(droplet_cts)
    if np.isnan(r_hat_val):
        print(np.nan)
        return np.nan
    elif r_hat_val == 1:
        # print(1)
        return 1
    else:

        lp_rhat = loglik_trinom_prof(droplet_cts, r_hat_val)
        chi_sq = sps.chi2.ppf(level, 1)

        try:
            root = spo.brentq(lambda r: 2*(loglik_trinom_prof(droplet_cts,r) - lp_rhat) + chi_sq,
                a=r_hat_val,
                b=0.9999
            )
            return root
        except ValueError:
            return np.nan


##------------------------------------------------------------------------------------------------##
if __name__ == '__main__':

    # q_dir = get_q_dir()

    wwtp_df = get_wwtp_df("data/wwtp_info.csv")

    os.chdir("/Volumes/PCR_Cowwid/01_dPCR_data/01_Stilla/NCX_CowwidVariants/04_Delta_Monitoring_Excel_new")
    excel_list = sorted(glob.glob("**/[!~]*.xlsx"))

    variant_df = pd.DataFrame()
    cols_to_use = ["ChamberName", "ChamberID", "Droplets"]
    for excel in excel_list:

        print(excel)

        results_df = pd.read_excel(excel,engine="openpyxl",sheet_name="Results")
        df_cols = [col for col in results_df.columns if re.search('|'.join(cols_to_use), col)]
        results_df = results_df[df_cols].copy()

        chamber_df = results_df["ChamberName"].str.rstrip().str.split("_")
        fn_split = excel.split("_")
        run_date = fn_split[0]
        run_date = f"{run_date[:4]}-{run_date[4:6]}-{run_date[6:]}"

        target = fn_split[2]
        if target in ("Delta","L452R","157-158del"):
            target = "S:L452R (Delta)"
        elif target in ("69-70del","Omicron"):
            target = "69-70del (Alpha, Omicron)"

        sethr = fn_split[-1].split(".")[0]
        if (len(sethr) == 2) & sethr.isdigit():
            run_date = f"{run_date}-{sethr}"
        else:
            run_date = f"{run_date}-00"

        new_df = pd.DataFrame()

        for i in range(len(chamber_df)):
            new_dict = {
                "file_name":excel,
                "run_date":run_date,
                "PCR_ID": pd.NA,
                "target_variant":target,
                "wwtp": pd.NA,
                "sample_date":pd.NA,
                "replicate":"A",
                "dilution":"1x",
                "sample_type": "ww",
                # "volume_mL":68
            }

            string_to_parse = chamber_df.loc[i]

            new_dict["ChamberID"] = string_to_parse[0]
            # print(string_to_parse)
            new_dict = sample_name_parser(string_to_parse, new_dict)
            # print(new_dict)
            sample_date = new_dict["sample_date"]

            line_df = pd.DataFrame.from_dict(new_dict, orient="index").T

            new_df = new_df.append(line_df, ignore_index=True)

        results_df = pd.merge(new_df, results_df, on="ChamberID")
        results_df.columns = results_df.columns.str.replace("NumberOf","")

        results_df["SampleName"] = results_df["wwtp"].str.cat(
            results_df["sample_date"],
            sep="_",
            na_rep=results_df["run_date"]
        )

        channels = results_df.filter(regex="Droplets$").columns.to_list()
        results_df.filter(regex="Droplets$").columns.to_list()
        # print(channels)

        variant_df = pd.concat([variant_df,results_df], axis=0, ignore_index=True)
        # variant_df.drop(columns=variant_df.filter(regex="none_").columns, inplace=True)

    variant_df["wwtp"] = variant_df["wwtp"].apply(
        lambda x: wwtp_df.loc[x]["ARANAME"]
            if x in wwtp_df.index
            else pd.NA
    )

    variant_df.drop(columns=variant_df.filter(regex="_NegativeDroplets$").columns, inplace=True)

    with pd.ExcelWriter("VariantResults.xlsx", date_format="YYYY-MM-DD", engine="openpyxl") as writer:
        variant_df.to_excel(writer, sheet_name="Results", index=True)



##Upload to Google_sheet

    droplet_cols = [
        "wwtp","sample_date","target_variant",
        "TotalDroplets",
        "double_PositiveDroplets",
        "HEX_single_PositiveDroplets",
        "FAM_single_PositiveDroplets",
        "Cy5_single_PositiveDroplets",
    ]

    final_variant_df = variant_df[droplet_cols].copy()
    final_variant_df.dropna(axis=0,subset=["sample_date"],inplace=True)

    gsheet_groupby = final_variant_df.groupby(["wwtp","sample_date","target_variant"])
    final_variant_df = gsheet_groupby.sum(min_count=1)

    calc_cols = [
        "TotalDroplets",
        "FAM_single_PositiveDroplets",
        "HEX_single_PositiveDroplets",
        "Cy5_single_PositiveDroplets",
        "double_PositiveDroplets"
    ]

    calc_df = final_variant_df[calc_cols]

    final_variant_df["NegativeDroplets"] = calc_df.apply(calc_neg_droplets,axis=1)
    final_variant_df["PercMutation"] = calc_df.apply(calc_perc_mutation,axis=1)
    final_variant_df["r̂"] = calc_df.apply(r_hat,axis=1)
    # final_variant_df["LogLikelihoodTrinomProf_R=.99"] = calc_df.apply(loglik_trinom_prof,axis=1,ratio=1)
    final_variant_df["lower_ci"] = calc_df.apply(lower_e,axis=1,level=0.95)*100
    final_variant_df["upper_ci"] = calc_df.apply(upper_e,axis=1,level=0.95)*100


    save_dir = "data/variant_data.csv"

    final_variant_df.to_csv(save_dir)
    ## droplet_cts=calc_df.iloc[93]
    # r_hat(droplet_cts)
    # calc_df.apply(loglik_trinom_prof,axis=1,ratio=r_hat(droplet_cts))
    # #
    # lower_e(droplet_cts)
    # upper_e(droplet_cts)



    gsheet = gc.open("COWWID LAB SHEET")
    worksheet = gsheet.worksheet("variant_monitoring")
    worksheet.clear()

    gsdf.set_with_dataframe(worksheet, final_variant_df, include_index=True)

    # plot_df = final_variant_df.reset_index().copy()
    # plot_df["sample_date"]= pd.to_datetime(plot_df["sample_date"])

    # fig, axs = plt.subplots(2,3, figsize=(14,9), sharex=True, sharey=True)
    # plt.style.use('dark_background')
    # colors = ["k",'darkgray','blue','red','lightgreen','magenta','cyan','orange',"lavender"]
    #
    # for i,wwtp in enumerate(plot_df["wwtp"].unique()):
    #
    #
    #     axis = axs.ravel()[i]
    #     subplot_df = plot_df[plot_df['wwtp'] == wwtp].copy()
    #
    #     axis.set_title(wwtp, fontsize=10)
    #     axis.set_xlabel("",fontsize=10)
    #     axis.set_ylabel("",fontsize=10)
    #
    #     for tick in axis.get_xticklabels():
    #         tick.set_rotation(45)
    #     sns.lineplot(
    #         x=subplot_df["sample_date"],
    #         y=subplot_df["PercMutation"],
    #         hue=subplot_df["target_variant"],
    #         style=subplot_df["target_variant"],
    #         palette="hsv",
    #         markers=True,
    #         alpha=0.8,
    #         ax=axis
    #     )
    #
    #     for variant in subplot_df["target_variant"].unique():
    #         variant_df = subplot_df[subplot_df["target_variant"] == variant]
    #         axis.fill_between(
    #             x=variant_df["sample_date"],
    #             y1=variant_df["upper_ci"],
    #             y2=variant_df["lower_ci"],
    #             fc="white",
    #             lw=0,
    #             alpha=0.3,
    #         )
    #
    # plt.tight_layout()
    # plt.show()
    #
    #
    # os.chdir("/Volumes/KatahdinHD/ResilioSync/DATA/pydata/cowwid/pycode")

# def fix_later():
    # confint_df = variant_df[["NegativeDroplets", "FAM_single_PositiveDroplets","HEX_single_PositiveDroplets","double_PositiveDroplets"]]
    #
    # confint_df.apply(r_hat,axis=1)
    #
    #
    # confint_df.apply(loglik_trinom_prof,axis=1,ratio=1)
    #
    # confint_df.apply(lower_e,axis=1)



    # print(r_hat(val))
    # print(loglik_trinom_prof(val,ratio=0))
    # print(loglik_trinom_prof(val,ratio=r_hat(val)))
# final_df = variant_df[droplet_cols].copy()
    # test_df = variant_df[["SampleName"]+droplet_cols].dropna().reset_index(drop=True)
# with pd.ExcelWriter("VariantResults.xlsx", date_format="YYYY-MM-DD", engine="openpyxl") as writer:
#     final_df.to_excel(writer, sheet_name="Results", index=True)
    # test_df["DoubleNegDroplets"] = test_df["TotalDroplets"] - test_df[droplet_cols[1:]].sum(axis=1)
    # x=test_df[["DoubleNegDroplets","FAM_single_pos._PositiveDroplets","double_pos._PositiveDroplets"]]











##------------------------------------------------------------------------------------------------##


##------------------------------------------------------------------------------------------------##

##------------------------------------------------------------------------------------------------##

##------------------------------------------------------------------------------------------------##
