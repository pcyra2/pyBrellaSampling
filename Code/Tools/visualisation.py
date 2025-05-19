try:
    import nglview
except:
    pass
try:
    import ase
    from ase.io import read
except:
    pass
from pyBrellaSampling.Code.Tools.classes import Colour
import pyBrellaSampling.Code.Tools.io as io


import plotly

NottinghamBlue = Colour(16,38,59)
JubileeRed = Colour(185,28,46)
MandarinOrange = Colour(249,129,9)
RebelsGold = Colour(222,180,6)
PioneeringPink = Colour(215,51,108)
CivicPurple = Colour(121,45,133)
ForestGreen = Colour(0,95,54)
BramleyApple = Colour(147,213,0)
TrentTurquoise = Colour(55,180,176)
MalaysiaSkyBlue = Colour(0,155,193)
PortlandStone = Colour(250,246,239)

PlotColourList=[NottinghamBlue, JubileeRed, MandarinOrange, ForestGreen, TrentTurquoise, BramleyApple, RebelsGold, MalaysiaSkyBlue, PioneeringPink, CivicPurple,]

def VisXYZ(path: str):
    """Reads in a .xyz file and returns an nglview viewer object. View by running the object in jupyter notebook

    Args:
        path (str): _description_

    Returns:
        _type_: _description_
    """
    mol = read(path)
    view = nglview.show_ase(mol)
    return view

def Plot2D(x: list, y:list, ylabels:list, title:str, xlabel:str, ylabel:str, template="simple_white")->plotly.graph_objects.Figure():
    """
    Creates a 2D plot using Plotly library in Python.
    
    Args:
        x (list): A list of x-values for the plot.
        y (list): A list of lists where each sublist corresponds to y-values for one series.
        ylabels (list): A list of labels for each series, used as legend labels.
        title (str): The title of the plot.
        xlabel (str): Label for the x-axis.
        ylabel (str): Label for the y-axis.
    
    Returns:
        plotly.graph_objects.Figure: A Plotly Figure object containing the plotted data.
    """
    fig = plotly.graph_objects.Figure()
    j=0
    for i in range(len(y)):
        fig.add_scatter(x=x, y=y[i], 
                        name=ylabels[i], 
                        marker = {"color":PlotColourList[j].to_string(),})
        if j == len(PlotColourList)-1:
            j=0
        else:
            j +=1
    if title != None or title != "":
        fig.update_layout(title=title)
    fig.update_layout(xaxis_title = xlabel,
                    yaxis_title = ylabel,
                    template=template)
    return fig

def PlotTrajData(Data:dict, PlotKeys:list, Title:str, xlabel:str, ylabel:str):
    """
    Plots trajectory data from a dictionary. The dictionary format should be Data[XValue][YKey] = YValue where PlotKeys are a list of dictionary keys for Y plots. 

    Args:
        Data (dict): A dictionary containing the trajectory data with numeric keys and associated values.
        PlotKeys (list or str): The key(s) to plot on the y-axis. Can be a single string or a list of strings.
        Title (str): The title of the plot.
        xlabel (str): The label for the x-axis.
        ylabel (str): The label for the y-axis.

    Returns:
        Plot (plotly.graph_objects.Figure): The generated 2D plot with the trajectory data.
        XData (list): The X axis points
        YData (list): List of lists of Y data values corresponding to the Plot Keys
    """
    XData = [key for key in Data.keys() if type(key) == int or type(key) == float]
    if type(PlotKeys) == str:
        PlotKeys = [PlotKeys]
    YData = []
    for key in PlotKeys:
        Y = [Data[point][key] for point in XData]
        YData.append(Y)
    Plot = Plot2D(XData, YData, PlotKeys, Title, xlabel, ylabel)
    return Plot, XData, YData


def LatexPlottable(XData, YData, PlotKeys, path):
    lines = [str] * (len(XData)+1)
    lines[0] = "x"
    for key in PlotKeys:
        lines[0] = lines[0]+ "\t"+key 
    for i in range(len(XData)):
        lines[i+1] = str(XData[i])
        for j in range(len(YData)):
            lines[i+1] = lines[i+1] +"\t"+ str(YData[j][i])
    io.textDump(lines, path)