from bokeh.plotting import figure, show, output_file, ColumnDataSource
from bokeh.models import HoverTool, Whisker, Div,Legend,LegendItem,DatetimeTickFormatter,PrintfTickFormatter
from bokeh.layouts import row,column, gridplot
from bokeh.io import curdoc
from datetime import datetime

import pandas as pd
import numpy as np

data = "data/respv_data.csv"
respv_df = pd.read_csv(data)
respv_df.set_index("sample_date",inplace=True)
respv_df.index = pd.to_datetime(respv_df.index)



wwtp_dict={
    "ARA Werdhölzli":"Zürich",
    "STEP Aire":"Geneva",
    "ARA Altenrhein":"Altenrhein",
    "ARA Chur":"Chur",
    "CDA Lugano":"Lugano",
    "ARA Sensetal":"Laupen"
}
channels = ("IAV-M","IBV-M","RSV")
colors = ("red","lightgreen","orange")
labels = ("Influenza A", "Influenza B", "RSV")
TOOLS = "reset,xpan,xwheel_zoom,box_zoom,save"

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

plot_list = []
for i,wwtp in enumerate(wwtp_dict.keys()):
    subplot_df = respv_df[respv_df["wwtp"]==wwtp].copy()

    if i > 0:
        x_range = plot_list[0].x_range
        y_range = plot_list[0].y_range
    else:
        x_range = (pd.to_datetime("09-01-2022"), pd.to_datetime(datetime.now().date()))
        y_range = (0, 9e12)

    p = figure(
        plot_height=400,
        plot_width=500,
        # sizing_mode='stretch_both',
        x_axis_type="datetime",
        title=wwtp_dict[wwtp],
        tools=TOOLS,
        x_range=x_range,
        y_range=y_range
    )

    # figure(, tools='pan', id="blue_fig")
    p.title.text_font_size = '16pt'

    p.xaxis.axis_label = 'Sample Date'
    p.xaxis.axis_label_text_font_size = '10pt'
    p.xaxis.formatter=DatetimeTickFormatter(months = ['%m/%Y'])

    p.yaxis.axis_label = 'Viral Load (gc/day/100,000 people)'
    p.yaxis.axis_label_text_font_size = '10pt'
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
    p.legend.click_policy="hide"

    plot_list.append(p)

# legend1 = Legend(items=[("point1", plot_list[0])],
#                  orientation="horizontal")
#
# p.add_layout(legend1, 'above')
# renderer_list = []
# color_list = []
# for plot in plot_list:
#     for i in range(5):
#         color = choice(palette)
#         renderer = plot.line(range(10),random(10),line_width=2,color=color)
#         renderer_list += [renderer]
#         color_list += [color]
# legend_items = [LegendItem(label=color,renderers=[renderer for renderer in renderer_list if renderer.glyph.line_color==color]) for color in colors]
#
# ## Use a dummy figure for the LEGEND
# dum_fig = figure(plot_width=300,plot_height=600,outline_line_alpha=0,toolbar_location=None)
# # set the components of the figure invisible
# for fig_component in [dum_fig.grid[0],dum_fig.ygrid[0],dum_fig.xaxis[0],dum_fig.yaxis[0]]:
#     fig_component.visible = False
# # The glyphs referred by the legend need to be present in the figure that holds the legend, so we must add them to the figure renderers
# dum_fig.renderers += renderer_list
# # set the figure range outside of the range of all glyphs
# dum_fig.x_range.end = 1005
# dum_fig.x_range.start = 1000
# add the legend
# dum_fig.add_layout( Legend(click_policy='hide',location='top_left',border_line_alpha=0,items=legend_items) )

curdoc().theme = 'light_minimal'
output_file('../respv_dashboard/respv_dashboard.html',title="Respiratory Virus Monitoring Dashboard")

div1 = Div(text=
    """
    <style>
        body { background: #FFFFFF; }
    </style>
    <h1 style="color:black;">Influenza and Respiratory Syncytial Virus Prevalence in Swiss Wastewater</h1>
    """
)
logo1 = Div(text=
    """
    <img src='images/ealogo-black.png' style="width:200px;height:42.5;border:10px solid white;">
    """
)
logo2 = Div(text=
    """
    <img src='images/logo_epfl_black.png' style="width:169px;height:49;border:10px solid white;">
    """
)
note1 = Div(text=
    # """
    # <p style="color:pink;">A more inclusive Assay for the Omicron variant (ORF1a-Δ3675-3677), which detects both BA.1 and BA.2 subclades
    # has been test and will now be used in place of the original HV69-70 Assay</p>
    # <br>
    # <p style="color:white;">Note to the Viewer:</p>
    # <p style="color:white;">These data represent percent values for the Delta and Omicron (subclade
    # BA.1) variants, and should not be used to infer absolute numbers of SARS-CoV-2/hCoV-2019 in Wastewater
    # </p>
    # """
    """\n\n\n"""
)
footer = Div(text=
    """

    <!doctype html>
    <html lang="en">
        <head>
            <!-- Required meta tags -->
            <meta charset="utf-8">
            <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
            <meta name="Description" content="">
            <meta name="author" content="eawag">
            <meta name="Copyright" content="Copyright (c) 2020 eawag, all rights reserved.">
            <title>SARS-CoV-2 in Wastewater</title>

        </head>
        <body style="background-color:#efefef">

            <nav class="navbar " style="background-color: #118EC6;">
            <a class="navbar-brand" href="#" style="padding-left:30px">
                <font style="color:white;font-size:30px;padding-left:10px">SARS-CoV-2 in Wastewater</font>
            </a>
            <div class="row">
                <a class="navbar-brand" href="https://www.eawag.ch/en/" style="padding-left:30px">
                    <img src="static/img/eawag.svg" width="180" height="50" >
                </a>
                <a class="navbar-brand" href="https://actu.epfl.ch/" style="padding-right:30px">
                    <img src="static/img/epfl.svg" width="160" height="40" >
                </a>
            </div>
    """
)
header = row(footer, sizing_mode="stretch_both")
grid = gridplot(plot_list,ncols=3,toolbar_location="left",merge_tools=True,sizing_mode="scale_both")
show(column(footer,note1,grid),background="black")
