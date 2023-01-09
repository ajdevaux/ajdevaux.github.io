import io
import pandas as pd
from bokeh.embed import components
from bokeh.models import HoverTool
from bokeh.models import LinearAxis, Range1d
from bokeh.plotting import figure
from bokeh.resources import CDN
from jinja2 import Template


template = Template(
    '''<!DOCTYPE html>
        <html lang="en">
            <head>
                <meta charset="utf-8">
                <title>Overview</title>
                {{ resources }}
                {{ script }}
                <style>
                    .embed-wrapper {
                        display: flex;
                        justify-content: space-evenly;
                    }
                </style>
            </head>
            <body>
                <div>
                    {{ table }}
                </div>
                <div class="embed-wrapper">
                    {{ div }}
                </div>
            </body>
        </html>
        ''')

# df: pd.DataFrame = get_data_frame()
# table_html = df.to_html()

plot = figure(x_axis_label='time', y_axis_label='value', x_axis_type='datetime',
                plot_width=1600, plot_height=800,
                tools='pan,wheel_zoom,zoom_in,zoom_out,box_zoom,reset,save,hover,tap')
plot.sizing_mode = 'scale_width'
# now continue setup your plot
# ...
#

# get bokeh parts
script_bokeh, div_bokeh = components(plot)
resources_bokeh = CDN.render()

# render everything together
html = template.render(resources=resources_bokeh,
                       script=script_bokeh,
                       # table=table_html,
                       div=div_bokeh)

# save to file
out_file_path = "test.html"
with io.open(out_file_path, mode='w') as f:
    f.write(html)
