import dash
from dash import dcc, html
import plotly.graph_objs as go

# Initialize the Dash app
app = dash.Dash(__name__)

# Create data for the 3D bar chart
trace = go.Bar3d(
    x=[1, 2, 3, 4],  # Data for MPI Workers
    y=[1, 2, 4, 8],  # Data for OpenMP Threads
    z=[0, 0, 0, 0],  # Starting point
    dx=1,
    dy=1,
    dz=[2, 3, 5, 7],  # Height of the bars
)

layout = go.Layout(
    scene=dict(
        xaxis=dict(title="MPI Workers"),
        yaxis=dict(title="OpenMP Threads"),
        zaxis=dict(title="Seconds"),
    )
)

# Create the figure object
fig = go.Figure(data=[trace], layout=layout)

# Set the app layout
app.layout = html.Div([
    dcc.Graph(figure=fig)
])

# Run the app
if __name__ == "__main__":
    app.run_server(debug=True)
