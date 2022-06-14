#!/usr/bin/env python3
import matplotlib as mpl
mpl.use("Qt5Agg")
# import matplotlib.pyplot as plt

import os
import glob
import re
import pandas as pd
import numpy as np
import scipy.optimize as spo
import scipy.stats as sps
import seaborn as sns
import gspread
import gspread_dataframe as gsdf

from datetime import datetime

pd.set_option("display.width", 1000)
pd.options.display.max_rows = 2000
pd.options.display.max_columns = 15



##------------------------------------------------------------------------------------------------##
def sample_name_parser(string_to_parse, new_dict):
    j = 1
    sample_types = ("ntc","neg","pos","post","std","extcon","extc","ext","ec","exc",
                    "ef","sand","saliva","stock","pex")

    while j <= len(string_to_parse):
        try:
            substr = string_to_parse[j]

            substr = substr.lower().strip()

            if (len(substr) <= 2) & (substr.isdigit()):
                if len(substr) == 1: substr = substr.zfill(2)

                wwtp = substr
                new_dict["wwtp"] = wwtp

            elif substr in sample_types:
                if substr in sample_types[0:2]:
                    substr = "ntc"
                    new_dict["wwtp"] = substr
                elif substr in sample_types[2:5]:
                    substr = "pos"
                    new_dict["wwtp"] = substr
                elif substr in sample_types[5:10]:
                    substr = "exc"
                    new_dict["wwtp"] = substr

                new_dict["sample_type"] = substr

            elif (substr in ("2020","2021","2022")) & (len(substr) == 4):
                year = substr
                month = string_to_parse[j+1]
                dayrep = string_to_parse[j+2]

                if not dayrep.isdigit():
                    day = 00
                    j -= 1

                if dayrep[-1].upper() in ("A","B","C","D","I"):
                    day = dayrep[:-1]
                    rep = dayrep[-1]
                    new_dict["replicate"] = rep.upper()

                elif "to" in dayrep:
                    dayrep = dayrep.split(" ")

                    if len(dayrep) > 1:
                        day = dayrep[0]
                        dil = dayrep[1].split("to")[-1]
                    else:
                        day = dayrep[:2]
                        dil = dayrep.split("to")[-1]

                    new_dict["dilution"] = f"{dil}x"
                else:
                    day = string_to_parse [j+2]

                sample_date = f"{year}-{month}-{day}"
                new_dict["sample_date"] = sample_date
                j += 2
                pass
            elif (len(substr) == 1) & (substr.upper() in ("Q","W","Y","Z")):
                pcr_id = substr.upper()
                new_dict["PCR_ID"] = pcr_id

            elif substr in "spiked":
                new_dict["spiked"] = True
            elif (len(substr) == 1) & (substr.upper() in ("A","B","C","D","I")):
                new_dict["replicate"] = substr.upper()
            elif substr in ("ultraplex","up"):
                new_dict["pcr_mix"] = substr
            elif substr in ("qsv"):
                new_dict["pcr_mix"] = "qscript_virus"
            elif substr in ("BSA2ug"):
                new_dict["pcr_mix"] = new_dict["pcr_mix"] + "_2ugBSA"
            elif (substr.isalnum() & substr.endswith("x")):
                new_dict["dilution"] = substr
            elif substr in ("stericup","zymo"):
                new_dict["filt_protocol"] = substr
            elif substr in ("pmga","pmg2"):
                pmga_version = substr[-1]
                new_dict["volume_mL"] = 40
                new_dict["replicate"] = "C"
                if pmga_version == "a":
                    elution_offset_dil = int(new_dict["dilution"][:-1])*0.5
                    new_dict["dilution"] = f"{elution_offset_dil}x"
                # elif pmga_version == "2":
                #     new_dict["replicate"] = "C"
            elif substr in ("sendai","mhv","iav","ibv"):
                # new_dict["target"] = substr
                new_dict["sample_type"] = f"{substr}_stock"

            else:
                new_dict[f"unk_{j}"] = substr

            j += 1

        except IndexError:
            break

    return new_dict
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
        target_variant = droplet_cts.iloc[0] ##Assay name
        tot_ct = droplet_cts.iloc[1] #TotalDroplets
        dbl_ct = droplet_cts.iloc[2] #Double-positive droplets (either FAM+HEX or Cy5+HEX)
        fam_ct = droplet_cts.iloc[3] #FAM single-positive droplets
        hex_ct = droplet_cts.iloc[4] #HEX single-positive droplets
        cy5_ct = droplet_cts.iloc[5] #Cy5 single-positive droplets
    elif type(droplet_cts) == tuple:
        target_variant = droplet_cts.iloc[0] ##Assay name
        tot_ct = droplet_cts[1] #TotalDroplets
        dbl_ct = droplet_cts[2] #Double-positive droplets (either FAM+HEX or Cy5+HEX)
        fam_ct = droplet_cts[3] #FAM single-positive droplets
        hex_ct = droplet_cts[4] #HEX single-positive droplets
        cy5_ct = droplet_cts[5] #Cy5 single-positive droplets

    # if np.isnan(cy5_ct):##DELTA/OMICRON-BA2 Assay
    if (target_variant == "S:L452R (Delta)") | (target_variant == "S:L452R (Omicron-BA.4/5)"):##DELTA/OMICRON-BA4-5 Assay
        mutat_ct = dbl_ct #Mutation-positive droplets
        wldtp_ct = fam_ct #Mutation-negative, SARS-CoV-2 positive droplets ("wildtype")
        unkwn_ct = hex_ct #Unknown-yet-positive droplets; Right now unused
    # elif np.isnan(fam_ct):##OMICRON-BA1 Assay
    elif target_variant == "69-70del (Alpha, Omicron-BA.1)":##OMICRON-BA1 Assay
        mutat_ct = hex_ct #Mutation-positive droplets
        wldtp_ct = dbl_ct #Mutation-negative, SARS-CoV-2 positive droplets ("wildtype")
        unkwn_ct = cy5_ct #Unknown-yet-positive droplets; Right now unused
    elif target_variant == "ORF1a-d3675-3677":
        mutat_ct = cy5_ct #Mutation-positive droplets
        wldtp_ct = dbl_ct #Mutation-negative, SARS-CoV-2 positive droplets ("wildtype")
        unkwn_ct = hex_ct #Unknown-yet-positive droplets; Right now unused

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
# def calc_perc_mutation(droplet_cts):
#
#     tot_ct, mutat_ct, wldtp_ct, unkwn_ct = _parse_assay_counts(droplet_cts)
#
#     pos_ct = mutat_ct + wldtp_ct + unkwn_ct
#     perc_mut = mutat_ct / pos_ct*100
#
#     return perc_mut
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
def get_wwtp_df(wwtp_data_src):

    wwtp_df = pd.read_csv(wwtp_data_src,header=0)
    wwtp_df = wwtp_df[["ARA ID \n(Eawag/EPFL)", "ARANAME","Kanton","ARANR (GIS)",
                       "StadtNAME","Angeschlossene Einwohner\n(BAFU 2017)"]]
    wwtp_df.rename(columns = {"ARA ID \n(Eawag/EPFL)":"ARA_ID",
                              "Angeschlossene Einwohner\n(BAFU 2017)":"Population"}, inplace=True)
    wwtp_df.replace(to_replace="/",value="-",regex=True, inplace=True)
    wwtp_df["ARA_ID"] = wwtp_df["ARA_ID"].str.strip("_")
    wwtp_df["StadtNAME"] = wwtp_df["StadtNAME"].str.lower()

    wwtp_df.set_index("ARA_ID", inplace=True)

    return wwtp_df

##------------------------------------------------------------------------------------------------##
def calc_error(droplet_cts, mode, level=0.95):
    """
    Function derived from R, orig. author David Dreifuss:

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

    if mode == "lower":
        bound = 0
        a0 = 0.00001 ##Should be zero
        b0 = r_hat_val
    elif mode == "upper":
        bound = 1
        a0 = r_hat_val
        b0 = 0.99999 ##Should be one
    else:
        return ValueError

    if np.isnan(r_hat_val):
        return np.nan

    elif r_hat_val == bound:
        return bound

    else:
        lp_rhat = loglik_trinom_prof(droplet_cts, r_hat_val)
        chi_sq = sps.chi2.ppf(level, 1)

        try:
            root = spo.brentq(
                lambda r: 2*(loglik_trinom_prof(droplet_cts,r) - lp_rhat) + chi_sq,
                a=a0, b=b0
            )
            return root

        except ValueError:
            return 1

##------------------------------------------------------------------------------------------------##
def ddpcr_concentration(droplet_cts, droplet_volume_uL=0.000519):

    tot_ct, mutat_ct, wldtp_ct, unkwn_ct = _parse_assay_counts(droplet_cts)

    pos_ct = mutat_ct + wldtp_ct + unkwn_ct

    if pos_ct == 0:
        return 0
    else:
        return np.log(1 - (pos_ct / tot_ct)) * -(1/droplet_volume_uL)
##------------------------------------------------------------------------------------------------##


if __name__ == '__main__':

    gc = gspread.service_account(filename="cowwid-wave5-e31fa2771e32.json")
    wwtp_df = get_wwtp_df("data/wwtp_info.csv")

    os.chdir("/Volumes/PCR_Cowwid/01_dPCR_data/01_Stilla/NCX_CowwidVariants/04_Variant_Monitoring_Excels")
    excel_list = sorted(glob.glob("**/[!~]*.xlsx"))

    variant_df = pd.DataFrame()
    cols_to_use = ["ChamberName", "ChamberID", "Droplets"]
    for excel in excel_list:

        print(excel)

        results_df = pd.read_excel(excel,engine="openpyxl",sheet_name="Results")
        df_cols = [col for col in results_df.columns if re.search('|'.join(cols_to_use), col)]
        results_df = results_df[df_cols].copy()
        results_df.columns = results_df.columns.str.lower()

        chamber_df = results_df["chambername"].str.rstrip().str.split("_")
        fn_split = excel.split("_")
        run_date = fn_split[0]
        run_date = f"{run_date[:4]}-{run_date[4:6]}-{run_date[6:]}"
        sethr = fn_split[-1].split(".")[0]
        if (len(sethr) == 2) & sethr.isdigit():
            run_date = f"{run_date}-{sethr}"
        else:
            run_date = f"{run_date}-00"
        print(run_date)

        target = fn_split[2]
        if (target in ("Delta","L452R","157-158del")) & (run_date < "2022-06-01"):
            target = "S:L452R (Delta)"
            print("delta")
        elif (target in ("Delta","L452R","157-158del")) & (run_date > "2022-06-01"):
            target = "S:L452R (Omicron-BA.4/5)"
            print("omicron")
        elif target in ("69-70del","Omicron"):
            target = "69-70del (Alpha, Omicron-BA.1)"
        elif target in ("ORF1a-d3675-3677"):
            target = "ORF1a-d3675-3677"



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
            line_df.columns = line_df.columns.str.lower()


            new_df = new_df.append(line_df, ignore_index=True)

        results_df = pd.merge(new_df, results_df, on="chamberid")
        results_df.columns = results_df.columns.str.replace("numberof","")

        results_df["samplename"] = results_df["wwtp"].str.cat(
            results_df["sample_date"],
            sep="_",
            na_rep=results_df["run_date"]
        )

        channels = results_df.filter(regex="droplets$").columns.to_list()
        results_df.filter(regex="droplets$").columns.to_list()
        # print(channels)

        variant_df = pd.concat([variant_df,results_df], axis=0, ignore_index=True)
        # variant_df.drop(columns=variant_df.filter(regex="none_").columns, inplace=True)

    variant_df["wwtp"] = variant_df["wwtp"].apply(
        lambda x: wwtp_df.loc[x]["ARANAME"]
            if x in wwtp_df.index
            else pd.NA
    )

    variant_df.drop(columns=variant_df.filter(regex="_negativedroplets$").columns, inplace=True)

    # with pd.ExcelWriter(f"VariantResults_{datetime.now().date()}.xlsx", date_format="YYYY-MM-DD", engine="openpyxl") as writer:
    # variant_df.to_csv(f"VariantResults_{datetime.now().date()}.csv", index=True)
    # lock(f"VariantResults_{datetime.now().date()}.csv")
    os.chdir("/Users/dejavu/Data/cowwid/cowwid_website/variant_monitoring")
    ##Upload to Google_sheet

    droplet_cols = [
        "wwtp","sample_date","target_variant",
        "totaldroplets",
        "double_positivedroplets",
        "fam_single_positivedroplets",
        "hex_single_positivedroplets",
        "cy5_single_positivedroplets",
    ]

    final_variant_df = variant_df[droplet_cols].copy()
    final_variant_df.dropna(axis=0,subset=["sample_date"],inplace=True)
    final_variant_df = final_variant_df[final_variant_df["totaldroplets"]!=0]

    final_variant_df = final_variant_df.groupby([
        "wwtp",
        "sample_date",
        "target_variant"
    ]).sum(min_count=1).reset_index()

    calc_cols = [
        "target_variant",
        "totaldroplets",
        "double_positivedroplets",
        "fam_single_positivedroplets",
        "hex_single_positivedroplets",
        "cy5_single_positivedroplets",
        "double_positivedroplets",
    ]

    calc_df = final_variant_df[calc_cols]

    final_variant_df["negativedroplets"] = calc_df.apply(calc_neg_droplets,axis=1)

    final_variant_df["percmutation"] = calc_df.apply(r_hat,axis=1) * 100
    final_variant_df["lower_ci95"] = calc_df.apply(calc_error,axis=1,mode="lower",level=0.95)*100
    final_variant_df["upper_ci95"] = calc_df.apply(calc_error,axis=1,mode="upper",level=0.95)*100
    final_variant_df["hcov19concentration_(gc/rxn)"] = calc_df.apply(ddpcr_concentration,axis=1)
    final_variant_df.set_index("wwtp", inplace=True)
    final_variant_df.to_csv("data/variant_data.csv")

    gsheet = gc.open("COWWID LAB SHEET")
    worksheet = gsheet.worksheet("variant_monitoring")
    worksheet.clear()

    gsdf.set_with_dataframe(worksheet, final_variant_df, include_index=True)
