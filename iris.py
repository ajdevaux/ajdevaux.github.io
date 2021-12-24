import pandas as pd
from bokeh.sampledata.iris import flowers

from bokeh.plotting import figure, show, output_file
from bokeh.models import ColumnDataSource, HoverTool

colormap = {'setosa': 'red', 'versicolor': 'green', 'virginica': 'blue'}
flowers['colors'] = [colormap[x] for x in flowers['species']]

hover = HoverTool(tooltips=[
    ("Sepal length", "@sepal_length"),
    ("Sepal width", "@sepal_width"),
    ("Petal length", "@petal_length"),
    ("Species", "@species")
    ])

p = figure(title = "Iris Morphology", plot_height=500, plot_width=500, tools=[hover, "pan,reset,wheel_zoom"])

p.xaxis.axis_label = 'Petal Length'
p.yaxis.axis_label = 'Petal Width'

p.circle('petal_length', 'petal_width', color='colors',
         fill_alpha=0.2, size=10, source=ColumnDataSource(flowers))

output_file('index.html')

show(p)
