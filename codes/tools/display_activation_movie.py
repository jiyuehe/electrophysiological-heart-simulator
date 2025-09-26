import plotly.graph_objects as go # pip install plotly, pip install --upgrade nbformat. For 3D interactive plot: triangular mesh, and activation movie
import plotly.io as pio
pio.renderers.default = "browser" # simulation result mesh display in internet browser

def execute_on_mesh(vertex, face, map_color):
    # display movie (will open in internet browser)
    fig = go.Figure(
        data = [
            go.Mesh3d(
                x = vertex[:, 0], y = vertex[:, 1], z = vertex[:, 2],
                i = face[:, 0], j = face[:, 1], k = face[:, 2],
                vertexcolor = map_color[0])
        ],

        layout = go.Layout(
            updatemenus = [{
                "type": "buttons",
                "buttons": [
                    {"label": "Play", "method": "animate", "args": [None, {"frame": {"duration": 5}, "fromcurrent": True, "mode": "immediate"}]},
                    {"label": "Pause", "method": "animate", "args": [[None], {"frame": {"duration": 0}, "mode": "immediate", "transition": {"duration": 0}}]}
                ]
            }],

            sliders = [{
                'active': 0,
                'currentvalue': {"prefix": "Time: "},
                'pad': {"t": 50},
                'steps': [{'label': str(t), 'method': 'animate', 'args': [[str(t)], {'frame': {'duration': 0}, 'mode': 'immediate'}]}
                for t in range(len(map_color))]
            }]
        ),

        frames = [
            go.Frame(data = [go.Mesh3d(vertexcolor = map_color[t])], name = str(t))
            for t in range(len(map_color))
        ]
    )
    fig.show()