"""
CausalSetTheory_Geodesics.py module

This module compiles the main functions used for the CausetSimulation classes in the Causal_Set_Theory library. The main focus of this functions are 
the creation of the causet, and allows to work with them aftewards, analycing its causal estructure and determining geodesics within this framework.

Functions:
    SetCauset: Given a set of parameters for the simulation through a CausetSimulation class, this function creates a Causet.
    IsCausal: Given a manifold metric and two points, this function determines if they are causally connected.
    ChronologicalFuture: Given a causet, a manifold metric and a point, this functions computes all points in the causall future of the original.
    GetLinks: Given a CausetSimulation class, this functions creates a dictionary specifying the causal relation between points.
    GetGraph: Given a link dictionary describing the causal relation between points, this function creates a directed graph describing the same info.
    GetGeodesic: Given a causal dictionary of links, and two points, this function computes all posible geodesic between said points.

Author: Cano Jones, Alejandro
linkedin: www.linkedin.com/in/alejandro-cano-jones-5b20a7136
github: https://github.com/Cano-Jones
"""


#Libraries used
from __future__ import annotations #Dedicated to function typing
from typing import TYPE_CHECKING #Dedicated to function typing
if TYPE_CHECKING:
    from .Class_Objects import *

from numpy.random import poisson, uniform #Random generator for Causet sprinkling
from math import sin, asin, sinh, asinh, cos #Trigonometric functions
from scipy.optimize import newton #Root solver f(x)=0
from scipy.integrate import quad #Numerical integration
import networkx as nx #Graph analysis


###################################




def SetCauset(sim: CausetSimulation) -> set:
    """
    SetCauset function
        This function uses the "sprinkling algorithm to create a (1+1) causet from a set of parameters describing the model (through the 
        CausetSimulation class).
    
    Parameters:
        sim (CausetSimulation class): set of parameters describing the model to be considered.
    
    Returns:
        Causet (set): set of points in a (1+1) spacetime.
    """    
    
    #This rho parameter is the (const) point number density to be used throughout the manifold
    rho=sim.PointNumber/sim.Metric.ComputeVolume(sim.TimeRange, sim.SpaceRange) #Density = # of points / Volume 
    
    Causet=set() #Empty set to be added

    #The algorithm divides the spacetime manifold into a set of submanifolds (divisions)
    #We compute the temporal and spacial length of said divisions
    Delta_T = (sim.TimeRange[1]-sim.TimeRange[0])/sim.Divisions[0]
    Delta_S = (sim.SpaceRange[1]-sim.SpaceRange[0])/sim.Divisions[1]

    #Now we will cicle through each submanifold, computing its volume and determining how many (and where) points will have
    for i in range(sim.Divisions[0]):
        for j in range(sim.Divisions[1]):

            #This are the upper and lower bounds of the spacetime submanifold to be used
            TRange=(sim.TimeRange[0]+i*Delta_T, sim.TimeRange[0]+(i+1)*Delta_T)
            SRange=(sim.SpaceRange[0]+j*Delta_S, sim.SpaceRange[0]+(j+1)*Delta_S)

            Vol=sim.Metric.ComputeVolume(TRange, SRange) #Volume of submanifold
            N=poisson(rho*Vol) #The number of points inside this region of spacetime is given by a poisson random distribution

            #We "sprinkle" those N points inside the region, in a uniform random distribution, then add them to the causet
            for _ in range(N):
                Causet.add((uniform(TRange[0],TRange[1]), (uniform(SRange[0],SRange[1]))))

    return Causet


###################################
    

def IsCausal(Metric: MetricTensor, point1: tuple[float], point2: tuple[float]) -> bool:
    """
    IsCausal function:
        Given a metric manifold (1+1 FLRW) and two points, this function determines if point2 is in the causal future of point1 (i.e. if there is
        a timelike geodesic between them).

    Parameters:
        Metric (MetricTensor class): class description of the (1+1) FLRW metric manifold.
        point1 (2D float tuple): first spacetime point.
        point2 (2D float tuple): second spacetime point.
    
    Returns:
        _ (Bool): If point2 is in the causal future of point2, returns True (False otherwise).
    """

    #If the time coordinate of point1 is larger than the one of point2, the causal future connection is not correct
    if point2[0]<point1[0]: return False


    #To determine if there is a timelike geodesic between the two points, knowing if the scale factor is a constant will save computational power.
    def a_IsConstant(Metric) -> bool: 
        """
        a_IsConstant function
            Given a MetricTensor class, this function determines if the scale factor of the (1+1) FLRW metric is constant in time.
        Parameters:
            Metric (MetricTensor class): class describing the metric manifold.
        """

        from sympy import simplify, symbols #Symbolic computation library functions

        is_it = None #Boolean variable, if true, the scale factor will be constant

        #To determine if the scale factor is a constant, we will check if its derivative (saved into the MetricTensor class) is zero.

        try:
            is_it = (simplify(Metric.da(symbols('t')))==0) #Convert the derivative of scale factor to symbolic function and check if its zero
        except:
            is_it = False
            
        return is_it
    
    inv_a = lambda t: 1/Metric.a(t) #This function will become handy afterwards
    r"""
    To be able to determine if there is a timelike geodesic between the two points, we will trace a null geodesic starting at point1 and
    ending at the same temporal coordinate as point2; and then checking whether point2 is in the lightcone of point1. To be able to do that, we
    solve for the differential equation 
    $\frac{\td r}{\td t}=\pm\frac{c}{a(t)}\sqrt{1-\kappa r^2}\frac{\td r}{\td t}=\pm\frac{c}{a(t)}\sqrt{1-\kappa r^2}$
    which has different forms for each value of kappa (-1,0,1); described bellow.
    """
    
    if Metric.kappa == 0: #Flat space

        if a_IsConstant(Metric=Metric): #Minkowskian spacetime
            #Geodesics are linear

            DeltaR = point2[1]-point1[1] #Spacial distance between points
            DeltaTxinv_a = (point2[0]-point1[0])*inv_a(point1[0]) #Maximum distance travelled by light

            #If the distance between the two points is less than the maximum distance travelled by light, they must be causally conected.
            if DeltaR <= DeltaTxinv_a and DeltaR >= -DeltaTxinv_a:  return True
            else: return False

        else:
            #Similar method as previous analysis, this time, trajectories are not linear, integration is needed to account for scale factor
            DeltaR = point2[1]-point1[1]
            Int_inv_a = quad(inv_a, point1[0], point2[0])[0] #Maximum distance travelled by light

             #If the distance between the two points is less than the maximum distance travelled by light, they must be causally conected.
            if DeltaR <= Int_inv_a and DeltaR >= -Int_inv_a:  return True
            else: return False
    
    elif Metric.kappa == -1: #Hyperbolic space
        #solutions of trajectory of the form $\asin(r(t))=\pm \int \frac{\text{d}t}{a(t)}$
        if a_IsConstant(Metric=Metric):
            #If scale factor constant, easy integration, we compute the maximum distance for a point to have travelled (each direction)
            Right = sinh(asinh(point1[1])+(point2[0]-point1[0])*inv_a(point1[0]))
            Left = sinh(asinh(point1[1])-(point2[0]-point1[0])*inv_a(point1[0]))

            #If the second point is between those maximum distances, the two points are causally connected
            if point2[1] <= Right and point2[1] >= Left: return True
            else: return False

        else:
            #Non constant, numerical integration needed to compute maximum distance for light
            Right = sinh(asinh(point1[1])+quad(inv_a, point1[0], point2[0])[0])
            Left = sinh(asinh(point1[1])-quad(inv_a, point1[0], point2[0])[0])

            #If the second point is between those maximum distances, the two points are causally connected
            if point2[1] <= Right and point2[1] >= Left: return True
            else: return False
    
    elif Metric.kappa == 1: #Spherical space, can be circumnavigated
        r"""
        Solutions for positive curvature are of the form $\sin(r(t))=\pm\int\frac{\text{d}t}{a(t)}$. 
        The periodic nature of this space adds an extra problematic with previous methods. Here we will try to find the time for witch a light
        ray will arrive at the r=\pm 1 mark, then check wether point2 time coordinate is passed that time (then any timelike trajectory would be
        able to arrive there); if this is not the case, we will need to check case by case the relative position of point2 in relation to this 
        null geodesic.
        To be able to find the solution of r(t)=\pm 1 we will use the Newton-Raphson method, which uses the derivative of r(t).
        """

        if a_IsConstant(Metric=Metric):
            #Constant scale factor, easy integration

            #Functions f(t)=0 used for the Newton-Raphson method, each for r(t)=1 and r(t)=-1
            def F_Rifght(t):
                if t<point1[0]: return 1e6 #Penalization function to ensure that the solution found is at a future time
                return sin(asin(point1[1])+(t-point1[0])*Metric.a(t))-1
            def F_Left(t):
                if t<point1[0]: return 1e6 #Penalization function to ensure that the solution found is at a future time
                return sin(asin(point1[1])-(t-point1[0])*Metric.a(t))+1

            #Derivatives of f(t) used for the Newton-Raphson method, each for r(t)=1 and r(t)=-1
            def DerF_Right(t):
                return cos(asin(point1[1])+(t-point1[0])*Metric.a(t))*inv_a(t)
            def DerF_Left(t):
                return -cos(asin(point1[1])-(t-point1[0])*Metric.a(t))*inv_a(t)

            #We usethe Newton-Raphson method to find the solution of r(t)=\pm 1
            t_right=newton(F_Rifght, point1[0], DerF_Right)
            t_left=newton(F_Left, point1[0], DerF_Left)

            #If point2 is passed the time at which r(t)=\pm 1, there must be a timelike geodesic connecting the two points
            if point2[0]>t_right and point2[0]>t_left: return True    

            #If point2 is not passed that time, we must check the relative position
            s_Right=sin(asin(point1[1])+(point2[0]-point1[0])*Metric.a(point2[0]))
            if point2[0]>t_left and point2[1]<s_Right: return True
            
            s_Left=sin(asin(point1[1])-(point2[0]-point1[0])*Metric.a(point2[0]))
            if point2[0]>t_right and point2[1]>s_Left: return True

        else:
            #Non constant scale factor, integration needed 

            #Functions f(t)=0 used for the Newton-Raphson method, each for r(t)=1 and r(t)=-1
            def F_Rifght(t):
                if t<point1[0]: return 1e6 #Penalization function to ensure that the solution found is at a future time
                return sin(asin(point1[1])+quad(inv_a, point1[0], t)[0])-1
            def F_Left(t):
                if t<point1[0]: return 1e6 #Penalization function to ensure that the solution found is at a future time
                return sin(asin(point1[1])-quad(inv_a, point1[0], t)[0])+1
            
            #Derivatives of f(t) used for the Newton-Raphson method, each for r(t)=1 and r(t)=-1
            def DerF_Right(t):
                return cos(asin(point1[1])+quad(inv_a, point1[0], t)[0])*inv_a(t)
            def DerF_Left(t):
                return -cos(asin(point1[1])-quad(inv_a, point1[0], t)[0])*inv_a(t)
            

            #We usethe Newton-Raphson method to find the solution of r(t)=\pm 1
            t_right=newton(F_Rifght, point1[0], DerF_Right)
            t_left=newton(F_Left, point1[0], DerF_Left)

            #If point2 is passed the time at which r(t)=\pm 1, there must be a timelike geodesic connecting the two points
            if point2[0]>t_right and point2[0]>t_left: return True    

            #If point2 is not passed that time, we must check the relative position
            s_Right=sin(asin(point1[1])+quad(inv_a, point1[0], point2[0])[0])
            if point2[0]>t_left and point2[1]<s_Right: return True
            
            s_Left=sin(asin(point1[1])-quad(inv_a, point1[0], point2[0])[0])
            if point2[0]>t_right and point2[1]>s_Left: return True
    
    return False #If no True value has been returned, then the two points cannot be causally conected


###################################



def ChronologicalFuture(point: tuple[float], Causet: set[tuple[float]], metric: MetricTensor):
    """
    ChronologicalFuture function:
        This function computes all points on a causet that are part of the causal future of another.
          (A point is not considered to be in its future).
    
    Parameters:
        point (2D float tuple): point for which its future wants to be computed.
        Causet (set of 2d float tuples): set of points that can be on the previous future.
        metric (MetricTensor class): manifold metric (1+1 FLRW) determining causal structure.
    """

    future = set() #Empty set of points

    for p in (Causet-{point}): #Loop wor each point in causet except for the point to be studied
        if IsCausal(metric, point, p): future.add(p) #if p and point are causally conected, we add p to the future set

    return future #Return the set of all points in the causal future of point.


###################################



def GetLinks(sim: CausetSimulation) -> dict:
    """
    GetLinks function:
        This function will describe the causal structure of a causet in a particular way: it will create a dictionary, in which the keys are
        all points 'p' in a causet, the corresponding value is a subset {q1,q2,q3...} of the causet such that there is no point g in the
        causet such that p<g<qi.
    Parameters:
        sim (CausetSimulation class): set of parameters describing the model to be considered.
    Returns:
        link_dic (dic): dictiornary encoding causal structure of a causet.
    """

    dic={} #Empty (auxiliary) dictionary to be completed
    
    for p in sim.Causet: #Loop for each point in causet
        dic[p]=ChronologicalFuture(p,sim.Causet,sim.Metric) #We add a new entry in the dictionary with key as a point, and value as a set of its causal future
    
    link_dic={} #Empty dictionary to be completed

    #We want to arrive at the proper subset where there are no intermediary causal points
    for p in sim.Causet:
        links=dic[p]
        #If a point is both in both future of another and in the future of a point in its future, then we discart it
        for caus in dic[p]:
            links=links.difference(dic[caus])
        link_dic[p]=links
    
    return link_dic #Final causal structure dictionary


##################################



def GetGraph(Links: dict[tuple, set]) -> nx.digraph:
    """
    GetGraph function:
        Obtains a directed graph describing the (direct) causal structure of a causet from a causal Links dictionary for better handling in
        future analysis.
    Parameters:
        Links (dict): dictionary describing (direct) causal structure of the causet.
    Returns:
        G (nx.digraph): Graph nx data structure describing direct causal structure of a causet. 
    """

    G=nx.DiGraph() #Empty graph data structure to be completed

    keys = Links.keys() #Causet

    for node in keys: #We convert each point of the causet into a graph node
        G.add_node(node)

    #For each pair of points in the causet, we add an edge between them in the graph if they are in direct causal connection.
    for node in keys: 
        for edge in Links[node]:
            G.add_edge(node, edge)

    return G


##################################



def GetGeodesic(links: dict[tuple, set], source: tuple[float], target: tuple[float]) -> list:
    """
    GetGeodesic function:
        In the theory of causal sets is postulated that geodesics correspond to the longest paths in a causal graph, this function computes
        said longest path (in general not unique).
    Parameters:
        links (dict): dictionary describing (direct) causal structure of the causet.
        source (2D float tuple): starting point in the geodesic.
        target (2D float tuple): ending point of the geodesic.
    Returns:
        longest_paths (list of lists of 2D tuples): list containing all possible geodesics between source and target.
    """

    G=GetGraph(links) #nx graph data structure encoding causal structure of the causet

    longest_paths = [] #Empty list of paths in the geodesic from  soruce to target

    
    paths = nx.all_simple_paths(G, source=source, target=target) #List of all posible paths between source and target

    longest_path_length = 0 #Auxiliary variable (encounter with longest path will upgrade this value)
    #In the following loop we will select only the paths with maximal length
    for path in paths:
        if len(path) > longest_path_length: #If the path has a larger length than previously found
            longest_path_length = len(path)
            longest_paths.clear()
            longest_paths.append(path)
        elif len(path) == longest_path_length: #If it has the same length as the longest
            longest_paths.append(path)

    return longest_paths
