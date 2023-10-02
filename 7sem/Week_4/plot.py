#!/usr/bin/python3

import numpy as np
import pandas as pd
import sys
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.interpolate import griddata


def main():
    
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <time_data_filename>")
        return

    df = pd.read_csv(sys.argv[1], sep=" ")
    x = np.array(df.Depth)
    y = np.array(df.Threads)
    z = np.array(df.Time)

    xi = np.linspace(x.min(), x.max(), 100)
    yi = np.linspace(y.min(), y.max(), 100)

    X, Y = np.meshgrid(xi, yi)

    Z = griddata((x, y), z, (X, Y), method='linear')


    figSurface = go.Figure(go.Surface(x=xi, y=yi, z=Z, opacity=0.75))
    figScatter = px.scatter_3d(df, x='Depth', y='Threads', z='Time')
    figScatter.update_traces(marker_size=9)

    figAll = go.Figure(data=figSurface.data + figScatter.data, layout=figScatter.layout)
    figAll.show()




if __name__ == '__main__':
    main()