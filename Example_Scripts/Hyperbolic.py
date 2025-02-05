"""
Hyperbolic.py

This example script shows the path-geodesic correspondence in a hyperbolic background for a timelike trajectory.
"""
from FLRW_CausetGeodes import *
from numpy import linspace

if __name__ == "__main__":

    #First we must specify the metric manifold
    g = MetricTensor(kappa=-1, a=lambda t: 1) # egative curvature, no expansion

    #Here we specify the simulation parameters 
    simCauset = CausetSimulation(Metric=g, TimeRange=(0,3), SpaceRange=(-3,3), PointNumber=350, Divisions=(100,250))

    #For the continuum geodesic we employ the following code
    simContinuum = ContinuumSimulation(Metric=g, source=(0,0), SpacialVelocity=-0.75, tau_span=linspace(0,100,1000), g_type='timelike')
    simContinuum.ComputeGeodesic(TimeRange=simCauset.TimeRange, SpaceRange=simCauset.SpaceRange) #We keep the data inside the spacetime region

    #Causet
    simCauset.CreateCauset() #Original causet sprinkling
    simCauset.Causet.add(simContinuum.source) #We add the starting point of the continuum geodesic
    simCauset.Causet.add(simContinuum.tarjet) #We add the ending point of the continuum geodesic
    simCauset.PrintCauset("Causet") #Visualization of the Causet

    #Causal Structure
    simCauset.GetLinks() #We compute the causal structure dictionary
    simCauset.PrintHasseDiagram("Hasse_Diagram") # Visualization of causal structure

    #Path-Geodesic correspondence
    simCauset.Geodesics(simContinuum.source, simContinuum.tarjet) #We compute the causet geodesic
    PrintComparison(simCauset.Causet, simCauset.Geodesic, simContinuum.Geodesic, "Correspondence") #Visualization of correspondence