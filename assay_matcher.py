import pandas as pd

data_dict = pd.read_excel("/Users/dejavu/Data/cowwid/data/N1MHV.ResultsDPCR.2022-01-07.xlsx",
                    engine="openpyxl",
                    sheet_name=[0,1],
                    index_col = 0)
print("\nFinally finished reading Excel because Openpyxl sucks\n")

full_df = pd.concat([data_dict[0],data_dict[1]],axis=0)



variant_df = pd.read_csv("/Users/dejavu/Data/cowwid/cowwid_website/variant_monitoring/data/variant_data.csv"

)

v_df = variant_df[["wwtp","sample_date","target_variant","hcov19concentration_(gc/rxn)"]]
new_df = full_df[["wwtp","sample_date",'SARS-N1_(gc/rxn)']].reset_index(drop=True)

new_df =new_df.groupby(["wwtp","sample_date"]).mean()

v_df = v_df.groupby(["wwtp","sample_date","target_variant"]).mean()
v_df.reset_index(level=2,inplace=True)

for target in variant_df["target_variant"].unique():
    a= pd.concat([new_df, v_df[v_df["target_variant"]==target]],axis=1)
    a.dropna(subset=["hcov19concentration_(gc/rxn)"],axis=0,inplace=True)
