from bokeh.plotting import figure, show, output_file, ColumnDataSource
from bokeh.models import HoverTool, Whisker, Div,Legend,LegendItem,DatetimeTickFormatter,PrintfTickFormatter
from bokeh.models.widgets import PreText
from bokeh.layouts import row,column, gridplot
from bokeh.io import curdoc
from datetime import datetime

import pandas as pd
import numpy as np
from math import pi

data = "data/respv_data.csv"
respv_df = pd.read_csv(data)
respv_df.set_index("sample_date",inplace=True)
respv_df.index = pd.to_datetime(respv_df.index)
curdoc().theme = 'light_minimal'

wwtp_dict={
    "ARA Werdhölzli":"Zürich",
    "STEP Aire":"Geneva",
    "ARA Altenrhein":"Altenrhein",
    "ARA Chur":"Chur",
    "CDA Lugano":"Lugano",
    "ARA Sensetal":"Laupen"
}
channels = ("IAV-M","IBV-M","RSV")#,"SARS-N2")
colors = ("red","lightgreen","orange","lightblue")
labels = ("Influenza A", "Influenza B", "RSV","SARS-CoV-2")
TOOLS = "reset,pan,xwheel_zoom,box_zoom,save"

for chan in channels:
    respv_df[f"{chan}_low-range"] = respv_df[f"{chan}_normLoad-(gc/d/1e5)"] - respv_df[f"{chan}_normLoad-range"]*0.5
    respv_df[f"{chan}_hi-range"] = respv_df[f"{chan}_normLoad-(gc/d/1e5)"] + respv_df[f"{chan}_normLoad-range"]*0.5

plot_props = {
    "x":"sample_date",

    "line_width":3
}

whisker_props = {
    "base":"sample_date",
    "line_width":2,
    "line_cap":"round",
    "upper_head":None,
    "lower_head":None
}

today = pd.to_datetime(datetime.now().date())
prev_month = today - pd.Timedelta(days=60)
plot_list = []
for i,wwtp in enumerate(wwtp_dict.keys()):
    subplot_df = respv_df[respv_df["wwtp"]==wwtp].copy()

    if i > 0:
        x_range = plot_list[0].x_range
        y_range = plot_list[0].y_range
    else:
        x_range = (prev_month,today)
        y_range = (0, 1e13)

    p = figure(
        plot_height=350,
        plot_width=800,
        # sizing_mode='stretch_both',
        x_axis_type="datetime",
        title=wwtp_dict[wwtp],
        tools=TOOLS,
        toolbar_location="above",
        # y_axis_type="log",
        x_range=x_range,
        y_range=y_range
    )

    # p.toolbar.autohide = True
    p.title.text_font_size = '20pt'

    p.xaxis.axis_label = 'Sample Date'
    p.xaxis.axis_label_text_font_size = '12pt'
    p.xaxis.formatter=DatetimeTickFormatter(months = ['%d-%m-%Y'])
    p.xaxis.major_label_orientation = pi/4

    p.yaxis.axis_label = 'Viral Load (gc/day/100,000 people)'
    p.yaxis.axis_label_text_font_size = '12pt'
    p.yaxis.formatter= PrintfTickFormatter(format="%5e")

    respv_src = ColumnDataSource(subplot_df)

    for i,chan in enumerate(channels):

        p.line(source=respv_src,color=colors[i],legend_label=labels[i],y=f"{chan}_normLoad-7dMed",**plot_props)
        p.circle(source=respv_src,color=colors[i], y=f"{chan}_normLoad-(gc/d/1e5)",size=2,**plot_props)

        p.add_layout(
            Whisker(
                source=respv_src,
                upper = f"{chan}_hi-range",
                lower = f"{chan}_low-range",
                line_color=colors[i],
                **whisker_props
            )
        )

    p.legend.location = "top_left"
    p.legend.title = "Respiratory Virus"
    p.legend.click_policy="hide"

    plot_list.append(p)


output_file('plots.html',title="Respiratory Virus Monitoring Dashboard")


# div1 = Div(text=
#     """
#     <style>
#
#     </style>
#     <h2 style="color:black;">Influenza and Respiratory Syncytial Virus Prevalence in Swiss Wastewater</h2>
#     """
# )
# logo1 = Div(text=
#     """
#     <img src='images/ealogo-black.png' style="width:200px;height:42.5;border:10px #4587c1;">
#     """
# )
# logo2 = Div(text=
#     """
#     <img src='images/logo_epfl_black.png' style="width:169px;height:49;border:10px #4587c1;">
#     """
# )
# note1 = Div(text=
#     # """
#     # <p style="color:pink;">A more inclusive Assay for the Omicron variant (ORF1a-Δ3675-3677), which detects both BA.1 and BA.2 subclades
#     # has been test and will now be used in place of the original HV69-70 Assay</p>
#     # <br>
#     # <p style="color:white;">Note to the Viewer:</p>
#     # <p style="color:white;">These data represent percent values for the Delta and Omicron (subclade
#     # BA.1) variants, and should not be used to infer absolute numbers of SARS-CoV-2/hCoV-2019 in Wastewater
#     # </p>
#     # """
#     """\n\n\n"""
# )

# footer = Div(text=
#     """
#
#     """
# )
# header = row(div1,logo1,logo2,note1,background="#4587c1",sizing_mode="stretch_width")#, sizing_mode="stretch_both")
grid = gridplot(plot_list,ncols=1,toolbar_location="above",merge_tools=True)#,sizing_mode="stretch_width")
show(grid)
