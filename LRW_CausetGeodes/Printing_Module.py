"""
Printing_Module.py module

This module contains all the graphical representation functions from the Causal_Set_Theory program. This visualization tools allow the user to
see all the results as a grafical representation of points and diagrams in the spacetime region of its choosing.

Functions:
    PrintCauset: This function represents the spacetime position of the points of the Causet.
    PrintHaseDiagram: This function draws the Hase Diagram of the poset.
    PrintComparison: This function draws the spacetime positions of the caset, and both continuum & discrete geodesics.
    PrintContinuumGeodesic: This function draws the continuum geodesic in a spacetime diagram.
    PrintFuture: This function draws positions of the points of the causet, in green the causal points related with a source, in red, the others.

Author: Cano Jones, Alejandro
linkedin: www.linkedin.com/in/alejandro-cano-jones-5b20a7136
github: https://github.com/Cano-Jones
"""

#Libraries used
from __future__ import annotations #Dedicated to function typing
from typing import TYPE_CHECKING #Dedicated to function typing
if TYPE_CHECKING:
    from .Class_Objects import *

import matplotlib.pyplot as plt #Grafical library
plt.style.use("FLRW_CausetGeodes/matplotlib_style.mplstyle") #Style of ploting (LaTex)




###################################



def PrintCauset(Causet: set=set(), directory: str = None) -> None: 
    """
    PrintCauset function

    Given a causet, dis function draws the position of said points on a spacetime diagram and then shows it on screen (saves it if a directory is given).

    Parameters:
        Causet (set): set of points withing a spacetime region.
        directory (str): Name of the directory for the image to be saved on (if none is given, the image is not saved, only showed).
    """

    #We split the data into two lists describig the spacial and temporal positions of the causet points
    X=[p[1] for p in Causet]
    Y=[p[0] for p in Causet]
    #The points are drawn in a scatterplot
    plt.scatter(X,Y, s=0.5, c='black')

    #Axis labeling and format
    plt.axis('equal')
    plt.xlabel(r"Space $(r)$")
    plt.ylabel(r"Time $(ct)$")

    #if a directory name is given, the image is saved there
    if directory != None: plt.savefig(directory, dpi=300, bbox_inches='tight')

    #The image is shown
    plt.show()


###################################################


def PrintHaseDiagram(links: dict={}, directory: str = None):
    """
    PrintHaseDiagram function
        This functions draws a grafical representation (Hase Diagram) of the causal relation between causet points. If an arrow connects two
        points, this represents that those two points are in direct (aka, there is no intermediate points with this property) causal relation.
        The causet points are drawn in their respective spacetime coordinates.
    
        Parameters:
            links (dict): dictionary describing causal links. For a given point (the key) the value of the dictionary is a set of points causaly connected to.
            directory (str): Name of the directory for the image to be saved on (if none is given, the image is not saved, only showed).

    """
    
    #In this loop we draw the arros connecting points
    for p in links.keys(): #For each point in the causet
        for link in links[p]: #For each point causally conected to point p
            plt.arrow(p[1], p[0], link[1]-p[1], link[0]-p[0], length_includes_head=True, alpha=0.1) #Arrow drawn

    #Axis labeling and format
    plt.axis('equal')
    plt.xlabel(r"Space $(r)$")
    plt.ylabel(r"Time $(ct)$")
    
    #if a directory name is given, the image is saved there
    if directory != None: plt.savefig(directory, dpi=300, bbox_inches='tight')

    #The image is shown
    plt.show()





###################################################





def PrintCausetGeodesic(Causet: set=set(), Geodesics: list=[[0],[0]], directory: str = None) -> None:
    """
    PrintCausetGeodesic function
        This functions draws all given geodesics between two points on a causet. All causet points are drawn in a spacetime diagram, where
        red arrows connect points of the geodesic. On top of that, initial and final points of the geodesic are marked with an 'X'.
    
    Parameters:
        Causet (set): set of points withing a spacetime region.
        Geodesics (list of 2D tuples): List of geodesics, each element of the list is a different geodesic (list of ordered causet points).
        directory (str): Name of the directory for the image to be saved on (if none is given, the image is not saved, only showed).
    """

    
    
    #We split the data into two lists describig the spacial and temporal positions of the causet points
    X=[p[1] for p in Causet]
    Y=[p[0] for p in Causet]
    #The points are drawn in a scatterplot
    plt.scatter(X,Y, s=0.5, c='black')

    #We plot each geodesic separately
    for geo in Geodesics: #Loop for each geodesic
        #We split the points of the geodesic into two lists
        X=[geo[i][1] for i in range(len(geo))]
        Y=[geo[i][0] for i in range(len(geo))]

        for i in range(len(geo)-1): #For each pair of connected points along the geodesic we plot an arrow
            plt.arrow(X[i], Y[i], X[i+1]-X[i], Y[i+1]-Y[i], color='red')
    
    #The initial and final points along the geodesics are defined as source and target respectively
    source=Geodesics[0][0]
    target=Geodesics[0][-1]
    #And then are marked with a red 'X'
    plt.scatter([source[1], target[1]], [source[0], target[0]], c='red', marker='X')


    #Axis labeling and format
    plt.axis('equal')
    plt.xlabel(r"Space $(r)$")
    plt.ylabel(r"Time $(ct)$")

    #if a directory name is given, the image is saved there
    if directory != None: plt.savefig(directory, dpi=300, bbox_inches='tight')

    #The image is shown
    plt.show()



###################################################




def PrintComparison(Causet: set=set(), CausGeodesics:list=[[[0],[0]]], ManiGeodesic:list=[[0],[0]], directory: str = None) -> None:
    """
    PrintComparison function
        Given a causet and two different kinds of geodesics (continuum and discrete), this function generates an grafic representation of
        how those types of geodesics compare from eachother (ie, are they close enough?). All other points of the causet are also drawn for
        perspective reasons.
    
    Parameters:
        Causet (set): set of points withing a spacetime region.
        CausGeodesics (list of 2D tuples): List of geodesics, each element of the list is a different geodesic (list of ordered causet points).
        ManiGeodesic (list of 2 lists of floats): List of points of a continuum geodesic.
        directory (str): Name of the directory for the image to be saved on (if none is given, the image is not saved, only showed).
        
    """

    
    #We split the data into two lists describig the spacial and temporal positions of the causet points
    X=[p[1] for p in Causet]
    Y=[p[0] for p in Causet]
    #The points are drawn in a scatterplot
    plt.scatter(X,Y, s=0.5, c='black')

    #We plot each geodesic separately
    for geo in CausGeodesics: #Loop for each geodesic
        #We split the points of the geodesic into two lists
        X=[geo[i][1] for i in range(len(geo))]
        Y=[geo[i][0] for i in range(len(geo))]

        for i in range(len(geo)-1): #For each pair of connected points along the geodesic we plot an arrow
            plt.arrow(X[i], Y[i], X[i+1]-X[i], Y[i+1]-Y[i], color='red')

    #We plot the continuum geodesic
    plt.plot(ManiGeodesic[1],ManiGeodesic[0], '--')

    #Axis labeling and format
    plt.axis('equal')
    plt.xlabel(r"Space $(r)$")
    plt.ylabel(r"Time $(ct)$")

    #if a directory name is given, the image is saved there
    if directory != None: plt.savefig(directory, dpi=300, bbox_inches='tight')

    #The image is shown
    plt.show()



###################################################



def PrintContinuumGeodesic(Geodesic: list[tuple[float]], directory: str = None) -> None:
    """
    PrintContinuumGeodesic function
        Given a list of points along a geodesic, this function draws said points conected through lines in a spacetime diagram.
    
    Parameters:
        Geodesic (list of 2 lists of floats): List of points of a continuum geodesic.
        directory (str): Name of the directory for the image to be saved on (if none is given, the image is not saved, only showed).

    """
    #We split the points along the geodesic into two lists for better handeling
    T, X = Geodesic
    plt.plot(X,T) #We plot those points

    #Axis labeling and format
    plt.axis('equal')
    plt.xlabel(r"Space $(r)$")
    plt.ylabel(r"Time $(ct)$")

    #if a directory name is given, the image is saved there
    if directory != None: plt.savefig(directory, dpi=300, bbox_inches='tight')
    
    #The image is shown
    plt.show()


###################################################



def PrintFuture(links: dict={}, future: set=set(), directory: str = None) -> None:
    """
    PrintFuture function
        This function allows the user to visualize the causal structure of a causet. It draws a Hase Diagram in a spacetime region; longside
        a color-coded causal relation of a particular point; in which all green points are (future) causally connected to it; and the rest are 
        represented in red.
    
    Parameters:
        links (dict): dictionary describing causal links. For a given point (the key) the value of the dictionary is a set of points causaly connected to.
        future (set): set of points that are in the causal future of a particular point.
        directory (str): Name of the directory for the image to be saved on (if none is given, the image is not saved, only showed).

    """

    #We get the complete Causet from the link dictionary
    causet = links.keys()

    #In this loop, the Hase Diagram is drawn
    for p in causet: #Loop for each point in the causet
        for link in links[p]: #If there is a causal (direct) conection between two points, an arrow is drawn
            plt.arrow(p[1], p[0], link[1]-p[1], link[0]-p[0], length_includes_head=True, alpha=0.1)
    
    #All points are ploted, following the color-coded rule (green if future, red if not)
    for p in causet:
        if p in future:
            plt.scatter(p[1], p[0], c="green")
        else:
            plt.scatter(p[1], p[0], c="red")

    #Axis labeling and format
    plt.axis('equal')
    plt.xlabel(r"Space $(r)$")
    plt.ylabel(r"Time $(ct)$")

    #if a directory name is given, the image is saved there
    if directory != None: plt.savefig(directory, dpi=300, bbox_inches='tight')

    #The image is shown
    plt.show()
