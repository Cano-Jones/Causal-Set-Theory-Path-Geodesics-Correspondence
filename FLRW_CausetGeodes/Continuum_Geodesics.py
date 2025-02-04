"""
Continuum_Geodesics.py module

This module contains the main functionalities needed to compute (numerically) geodesics in a continuous (1+1) FLRW metric manifold, charactericed by
a scale factor and a space curvature constant (-1,0,1).

Functions:
    Equations_Motion: This function serves the numerical integration method, encoding the differential equations of the geodesic.
    ComputeVt: Since temporal and spacial velocity are correlated, only one is truly needed, with this function we compute Vt from Vr.
    ComputeGeodesic: Using numerical integration methods, we solve the equations of motion of the geodesic.
    CutGeodesic: This function extracts all points of the geodesic inside a specific region of spacetime.

Author: Cano Jones, Alejandro
linkedin: www.linkedin.com/in/alejandro-cano-jones-5b20a7136
github: https://github.com/Cano-Jones
"""

#Libraries used
from __future__ import annotations #Dedicated to function typing
from typing import TYPE_CHECKING #Dedicated to function typing
if TYPE_CHECKING:
    from .Class_Objects import *
    
from scipy.integrate import odeint #Numerical integration methods
from math import sqrt #Usual square root


###################################


def Equations_Motion(state: list[float], tau: float, Metric: MetricTensor) -> list[float]:
    """
    Equations_Motion function:
        This function encodes the equations of motion for geodesics in a (1+1) FLRW metric manifold. TThe two second order equations of motion (
        one for time, one for space) are separated into four first order equations (including equations for the time derivatives); this last 
        set of equations being the result of this function, used in a later function for numerical integration.

    Parameters:
        state (4D list float): initial state of the sistem (temporal and spacial coordinates, alongside temporal and spacial time derivatives).
        tau (float): proper time parameter (not used).
        Metric (MetricTensor class): class object encoding parameters of the metric manifold.
    Returns:
        _ (4D list floats): list of values of derivatives of the initial state (equations of motion).
    """

    #Characteristics of the metric manifols used in the equations of motion
    a = Metric.a #Scale factor
    da = Metric.da #Time derivative of the scale factor
    k = Metric.kappa #Spacial curvature
    
    #We separate the list of the state into different variables
    x1, y1, x2, y2 = state #[x1,y1,x2,y2]=[t,r,dt,dr]

    #Equations of motion
    dx2 = -da(x1)*a(x1)*y2*y2/(1-k*y1*y1) #ddt
    dy2 = -k*y2*y2*y1/(1-k*y1*y1)-2*da(x1)*x2*y2/a(x1) #ddr

    return [x2, y2, dx2, dy2] #Values of the derivatives of the state



###################################



def ComputeVt(Metric: MetricTensor, source: tuple[float], Vr: float, g_type: str) -> float:
    """
    ComputeVt function:
        Given that temporal and spacial velocities are related (different relation for time,light or spacial -like geodesics), one can compute one
        from the other. The choice in this library is to specify the spacial velocity and to compute the temporal velocity using this function.
    Parameters:
        Metric (MetricTensor class): class object encoding the (1+1) FLRW metric manifold.
        source (2d float tuple): a point in spacetime alongside a geodesic.
        Vr (float): spacial velocity of said point alongside a geodesic.
        g_type (str): type of geodesic (timelike, lightlike or null, spacelike).
    """
    #Parameters of the metric manifol
    a = Metric.a #Scale factor
    kappa = Metric.kappa #Space curvature

    t0,x0 = source # Spacetime coordinates of the point

    #Depending on the type of geodesic, the relation between Vr and Vt is different, the relation are listed bellow

    if g_type == 'timelike':
        return sqrt(a(t0)*a(t0)*Vr*Vr/(1-kappa*x0*x0)+1)
            
    elif type == 'null' or g_type == 'lightlike':
        return a(t0)*Vr/sqrt(1-kappa*x0*x0)
            
    else:
        return sqrt(a(t0)*a(t0)*Vr*Vr/(1-kappa*x0*x0)-1)


###################################



def ComputeGeodesic(source: tuple[float], Vr: float, Vt: float, t_span: list, Metric: MetricTensor) -> list:
    """
    ComputeGeodesic function:
        Numerical integration method designed to compute all points alongside a specific geodesic.

    Parameters:
        source (2D float tuple): initial spacetime coordinates of the geodesic.
        Vr (float): initial spacial velocity of source.
        Vr (float): initial temporal velocity of source.
        t_span (list of float): list of proper time values to compute geodesic from.
        Metric (MetricTensor class): encoded information of the (1+1) FLRW metric manifold.
    Returns:
        _ (2D list of lists of points): List of spacetime coordinates alongside the geodesic.
    """
    #We will use the scipy.odeint integration method, which  needs the initial state of the system in the following structure
    S0 = [source[0], source[1], Vt, Vr]

    #The integration of the geodesic is done in the following line
    solution = odeint(Equations_Motion, S0, t_span, (Metric,))
    
    return solution[:, 0], solution[:, 1] #We are interested in the temporal and Spacial coordinates alongside the geodesic


###################################


def CutGeodesic(Geodesic: list[list[float]], TimeRange: tuple[float], SpaceRange: tuple[float]) -> list:
    """
    CutGeodesic function:
        Given a set of points representing a geodesic in a (1+1) FLRW metric manifold, this function aims to "clean" the data in such a way that
        only points of the geodesic inside a particular region of spacetime is considered, thus "cutting it".
    
    Parameters:
        Geodesic (2D list of float lists): List of spacetime coordinates alongside a geodesic.
        TimeRange (2D float tuple): Range of times for which to consider points.
        SpaceRange (2D float tuple): Range of space for which to consider points.
    
    Returns:
        _ (2D list of float lists): List of spacetime coordinates alongside a geodesic inside a region of spacetime.
    """

    T , X = Geodesic #We separate the temporal and spacial coordinates alongside the geodesic
    
    #Two empty lists that will be added points of the geodesic inside the region of spacetime
    T_new = []
    X_new = []

    #We add to the previous lists, the points inside the TimeRange section
    for t in range(len(T)):
        if T[t]>=TimeRange[0] and T[t]<=TimeRange[1]:
            T_new.append(T[t])
            X_new.append(X[t])
    
    #Original Points of the geodesic are rewriten so that only the ones inside TimeRange are considered
    T = T_new
    X = X_new

    #New lists to be considered for the spacerange
    T_new = []
    X_new = []

    #We extract again, only the points inside SpaceRange
    for x in range(len(X)):
        if X[x]>=SpaceRange[0] and X[x]<=SpaceRange[1]:
            T_new.append(T[x])
            X_new.append(X[x])

    return T_new, X_new #The ponts left must be inside TimeRange and SpaceRange










__all__ = ['ComputeGeodesic', 'CutGeodesic', 'ComputeVt']