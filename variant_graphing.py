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
    "y":"PercMutation",
    "line_width":3
}

whisker_props = {
    "base":"sample_date",
    "upper":"upper_ci",
    "lower":"lower_ci",
    "line_width":2,
    "line_cap":"round",
    "upper_head":None,
    "lower_head":None
}

plot_list = []
for i,wwtp in enumerate(wwtp_list):
    subplot_df = variant_df[variant_df["wwtp"]==wwtp].copy()

    delta_df = subplot_df[subplot_df["target_variant"] ==  "S:L452R (Delta)"]
    omicron_df = subplot_df[subplot_df["target_variant"] ==  '69-70del (Alpha, Omicron)']
    # colormap = {'S:L452R (Delta)': 'lightgreen', '69-70del (Alpha, Omicron)': 'cyan'}
    # subplot_df['colors'] = [colormap[x] for x in zh_df['target_variant']]
    if i > 0:
        x_range = plot_list[0].x_range
        y_range = plot_list[0].y_range
    else:
        x_range = (pd.to_datetime("08-01-2021"), pd.to_datetime(datetime.now().date()))
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
    omi_src = ColumnDataSource(omicron_df)


    p.line(source=delta_src,color="lightgreen",legend_label="S:L452R (Delta/B.1.617.2)",**plot_props)
    p.circle(source=delta_src,color="green", size=6,**plot_props)

    p.line(source=omi_src,color="orange",legend_label="HV69-70 (Omicron/B.1.1.529)",**plot_props)
    p.circle(source=omi_src,color="darkorange", size=6,**plot_props)


    p.add_layout(
        Whisker(
            source=delta_src,
            line_color="green",
            **whisker_props
        )
    )
    p.add_layout(
        Whisker(
            source=omi_src,
            line_color="darkorange",
            **whisker_props
        )
    )

    p.legend.location = "center_left"
    # p.legend.click_policy="hide"

    plot_list.append(p)
# p.circle('petal_length', 'petal_width', color='colors',
#          fill_alpha=0.2, size=10, source=ColumnDataSource(flowers))
# data['y_fixed'] = [1200 - val for val in data.y]
# x = zh_df.sample_date
# y = zh_df.PercMutation
# pc = particle_df.pc
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
show(column(header,grid,footer),background="black")
