from bokeh.plotting import figure, show, output_file, ColumnDataSource
from bokeh.models import HoverTool, Whisker, Div
from bokeh.layouts import row,column, gridplot
from bokeh.io import curdoc
from datetime import datetime

import pandas as pd
import numpy as np

data = "data/variant_data.csv"
variant_df = pd.read_csv(data)
variant_df.set_index("sample_date",inplace=True)
variant_df.index = pd.to_datetime(variant_df.index)

wwtp_list=["ARA Werdhölzli","STEP Aire","ARA Altenrhein","ARA Chur","CDA Lugano","ARA Sensetal"]

TOOLS = "reset,xpan,xwheel_zoom,xbox_zoom,save"

plot_props = {
    "x":"sample_date",
    "y":"percmutation",
    "line_width":3
}

whisker_props = {
    "base":"sample_date",
    "upper":"upper_ci95",
    "lower":"lower_ci95",
    "line_width":2,
    "line_cap":"round",
    "upper_head":None,
    "lower_head":None
}

plot_list = []
for i,wwtp in enumerate(wwtp_list):
    subplot_df = variant_df[variant_df["wwtp"]==wwtp].copy()

    delta_df = subplot_df[subplot_df["target_variant"] ==  "S:L452R (Delta)"]
    omicron1_df = subplot_df[subplot_df["target_variant"] ==  "69-70del (Alpha, Omicron-BA.1)"]
    omicron2_df = subplot_df[subplot_df["target_variant"] ==  "ORF1a-d3675-3677"]

    if i > 0:
        x_range = plot_list[0].x_range
        y_range = plot_list[0].y_range
    else:
        x_range = (pd.to_datetime("11-01-2021"), pd.to_datetime(datetime.now().date()))
        y_range = (-5, 105)

    p = figure(
        plot_height=400,
        plot_width=500,
        x_axis_type="datetime",
        title=f"{wwtp}",
        tools=TOOLS,
        x_range=x_range,
        y_range=y_range

    )
    p.title.text_font_size = '16pt'

    p.xaxis.axis_label = 'Sample Date'
    p.xaxis.axis_label_text_font_size = '10pt'

    p.yaxis.axis_label = 'Percent hCoV-2019 Mutation'
    p.yaxis.axis_label_text_font_size = '10pt'

    delta_src = ColumnDataSource(delta_df)
    omi1_src = ColumnDataSource(omicron1_df)
    omi2_src = ColumnDataSource(omicron2_df)

    p.line(source=delta_src,color="lightgreen",legend_label="S:L452R (Delta)",**plot_props)
    p.circle(source=delta_src,color="green", size=6,**plot_props)

    p.line(source=omi1_src,color="orange",legend_label="HV69-70 (Omicron-BA.1)",**plot_props)
    p.circle(source=omi1_src,color="darkorange", size=6,**plot_props)

    p.line(source=omi2_src,color="magenta",legend_label="ORF1a-Δ3675-3677",**plot_props)
    p.circle(source=omi2_src,color="red", size=6,**plot_props)

    p.add_layout(
        Whisker(
            source=delta_src,
            line_color="green",
            **whisker_props
        )
    )
    p.add_layout(
        Whisker(
            source=omi1_src,
            line_color="darkorange",
            **whisker_props
        )
    )
    p.add_layout(
        Whisker(
            source=omi2_src,
            line_color="red",
            **whisker_props
        )
    )
    p.legend.location = "center_right"
    # p.legend.click_policy="hide"

    plot_list.append(p)

curdoc().theme = 'dark_minimal'
output_file('variant_monitoring.html',title="CoWWID Variant Monitoring")

div1 = Div(text=
    """
    <style>
        body { background: #0F0F0F; }
    </style>
    <h1 style="color:white;">Coronavirus hCoV-2019 Variant Prevalence in Swiss Wastewater</h1>
    """
)
logo1 = Div(text=
    """
    <img src='images/ealogo-white.png' style="width:200px;height:42.5;border:10px solid black;">

    <img src='images/logo_epfl_white.png' style="width:169px;height:49;border:10px solid black;">
    """
)
note1 = Div(text=
    """
    <p style="color:pink;">A more inclusive Assay for the Omicron variant (ORF1a-Δ3675-3677), which detects both BA.1 and BA.2 subclades
    has been test and will now be used in place of the original HV69-70 Assay</p>
    <br>
    <p style="color:white;">Note to the Viewer:</p>
    <p style="color:white;">These data represent percent values for the Delta and Omicron (subclade
    BA.1) variants, and should not be used to infer absolute numbers of SARS-CoV-2/hCoV-2019 in Wastewater
    </p>
    """
)
footer = Div(text=
    """
    <p style="color:white;">These data made possible by the tireless work of:</p>
    <ul style="color:white;">
        <li>Lea Caduff</li>
        <li>Pravin Ganesanandamoorthy</li>
        <li>Charlie Gan</li>
        <li>Johannes Rusch</li>
        <li>Franziska Böni</li>
        <li>David Dreifuss</li>
        <li>A.J. Devaux (<a href = "mailto:alexander.devaux@eawag.ch">Webmaster</a>)</li>
        <li>And the rest of the Eawag & EPFL Team</li>
    </ul>
    """
)
header = row(div1,logo1, sizing_mode="scale_width")
grid = gridplot(plot_list,ncols=2,toolbar_location="left")
show(column(header,note1,grid,footer),background="black")
