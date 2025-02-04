"""
Class_Objects.py Module

This module compiles all class objects used in the Causal_Set_Theory library. This classes describe (in that particular order) the
metric tensor used in the manifold (1+1 FLRW), main description and utilities for causets & main description and utilities of continuum
geodesics.

Classes:
    MetricTensor: (1+1) FLRW metric descriptor.
    CausetSimulation: description & utilities for Causets.
    ContinuumSimulation: description & utilities for continuum geodesics.

Author: Cano Jones, Alejandro
linkedin: www.linkedin.com/in/alejandro-cano-jones-5b20a7136
github: https://github.com/Cano-Jones
"""


#Used libraries & Methods
from scipy.integrate import dblquad #Computation of volume
from sympy.utilities.lambdify import lambdify #Derivative definition
from sympy import symbols #Derivative definition
from .CausalSetTheory_Geodesics import SetCauset, GetLinks, ChronologicalFuture, GetGeodesic #Causet utilities
from .Printing_Module import PrintCauset, PrintHaseDiagram, PrintContinuumGeodesic, PrintFuture #Printing utilities
from .Continuum_Geodesics import ComputeGeodesic, ComputeVt, CutGeodesic #Continuum utilities
from math import sqrt #square root, volume computation


#################################################################

class MetricTensor():
    """
    Class MetricTensor

    This object completely describes the (1+1) Friedmann-Lemaître-Robertson-Walker (FLRW) metric through its constant space curvature 'k' and scale factor 'a(t)'. The chosen coordinates are those such that k=-1,0,1 , a(t) has lenth units and c=1.

    Class attributes:
        kappa (int): space curvature constant (-1,0,1)
        a (callable): scale factor (time-coordinate function).
        da (callable): time derivative of scale factor (time-coordinate function)
    
    Class methods:
        ComputeVolume: Given the spacetime boundary conditions of a given region, returns the volume of said region.
        Derivative: Given a one parameter symbolic function, returns its derivative as a numerical callable function.
    """

    def __init__(self, kappa: int = 0, a: 'callable' = lambda t: 1) -> None:
        """
        Constructor for MetricTensor class

        Parameters (Defaults to Minkowski k=0, a(t)=1)
            kappa (int): space curvature constant (-1,0,1)
            a (callable symbolic): scale factor (time-coordinate function) written with sympy functions.
        """
        self.kappa=kappa # Space curvature constant
        self.a=lambdify(symbols('t'), a(symbols('t')), 'math') # Scale factor (converted from Symbolic to numerical function)
        self.da = self.Derivative(a) # Time derivative of the scale factor (converted to numerical function) 
    
    
    def ComputeVolume(self, TimeRange: tuple[float],
                   SpaceRange: tuple[callable]) -> float:
        """
        ComputeVolume method

        This method computes the volumen of a (1+1) region of spacetime described by its space and time boundaries. To correctly describe any region, the space boundaries can be a function of the time coordinates.

        Parameters:
            TimeRange (2D float tuple): (2D) Tuple consisting in the lower and upper time bounds of region.
            SpaceRange (2D callable tuple): (2D) Tuple consisting in the lower and upper space bounds of region (can be callable, through a function of time).

        Returns:
            Volume (float): Volume of region of spacetime described by TimeRange & SpaceRange.
        """
        # To compute the volume of the region we must compute the integral considering the correct boundaries and volume element
    
        def det_g(r,t): #Volume element
            # The volume element of a metric is given by the square root of the tensor deteminant times -1
            k = self.kappa
            a = self.a
            return a(t)/sqrt(1-k*r*r) # For a FLRW metric, this is the correct volume element

        # We extract the boundaries for future integration

        T_inf, T_sup = TimeRange # Lower and upper time bounds
        #We will consider the space boundaries as function of time, next lines make sure that this boundaries are callables
        S_inf = SpaceRange[0] if callable(SpaceRange[0]) else lambda t: SpaceRange[0] # Lower space bound
        S_sup = SpaceRange[1] if callable(SpaceRange[1]) else lambda t: SpaceRange[1] # Upper space bound

        # The integration is done using the scipy.integrate.dblquad function
        return dblquad(det_g, T_inf, T_sup, S_inf, S_sup)[0] #We return only the zeroth element, the first would be estimated error
    


    def Derivative(self, f: 'callable') -> 'callable':
        """
        Derivative method

        This method computes the derivative of a symbolic callable function 'f' and converts it into a callable numerical function.

        Parameters:
            f (callable symbolic): function to derivate, written with sympy library functions.
        
        Returns:
            Derivative (callable): derivative of the 'f' parameter function, as numerical callable.
        
        """
        
        from sympy import diff

        f_expression = f(symbols('t')) # Symbolic expression of function 'f' with parameter 't'
        derivative = diff(f_expression, symbols('t')) # Symbolic derivative of function 'f' with respect to 't'
        return lambdify(symbols('t'), derivative, 'math') # Returns the result as a callable written with 'math' library



##################################################


class CausetSimulation():
    """
    Class Simulation

    This object encapsulates the needed attributes & methods to create and describe a Causet and its causal relation through its links.

    Class attributes:
        Metric (MetricTensor class): Describes the metric manifold to be considered in the causet sprinkling
        TimeRange (2D float tuple): Describes the upper and lower time limits of the simulation
        SpaceRange (2D float tuple): Describes the upper and lower space limits of the simulation
        PointNumber (int): (Average) number of points in the causet to be generated
        Divisions (2D int tuple): Number of divisions of the spacetime range to be considered in the sprinkling
        Causet (set): Causal set of spacetime points
        Links (Dict): Dictionary containing the (direct) future of a given point 
        Geodesic([T,X] list, where X & Y are float lists): List describig spacetime coordinates along the geodesic 
    
    Class methods:
        CreateCauset: Creates a causet from the given class attributes and saves it into the Causet attribuite of the class
        GetLinks: Creates de Links dictionary and saves it into the Links attribute of the class
        PrintCauset: Prints a 2D scatter plot of the causet spacetime diagram and saves it on a given directory.
        PrintHaseDiagram: Prints the Hase Diagram of the Causet and saves it on a given directory
        PrintFuture: Prints the Hase Diagram of the Causet signaling the chronological future of a point, saves it on a given directory
        Geodesics: Computes all posible geodesics between two points.

    """
    def __init__(self, Metric: MetricTensor, TimeRange: tuple[float], SpaceRange: tuple[float],
                 PointNumber: int, Divisions: tuple[int]) -> None:
        """
        Constructor for CausetSimulation class

        Parameters:
            Metric (MetricTensor class): Describes the metric manifold to be considered in the causet sprinkling
            TimeRange (2D float tuple): Describes the upper and lower time limits of the simulation
            SpaceRange (2D float tuple): Describes the upper and lower space limits of the simulation
            PointNumber (int): (Average) number of points in the causet to be generated
            Divisions (2D int tuple): Number of divisions of the spacetime range to be considered in the sprinkling
        """
        
        self.Metric = Metric #MetricTensor class atribute
        self.TimeRange = TimeRange # 2D tuple describing the upper and lower time limits of the simulation
        self.SpaceRange = SpaceRange #2D tuple describing the upper and lower space limits of the simulation
        self.PointNumber = PointNumber #(Average) number of points in the causet to be generated
        self.Divisions = Divisions #Number of divisions of the spacetime range to be considered

        self.Causet = set() #Causal set of spacetime points
        self.Links = {} #Dictionary containing the (direct) future of a given point
        self.Geodesic = None #Points of the geodesic (to be computed)

    
    def CreateCauset(self) -> None:
        """
        CreateCauset method

        This method generates a causet withing a given spacetime region by the "Sprinkling" Poisson distribution method
       
        """
        self.Causet = SetCauset(self)
    
    def GetLinks(self) -> None:
        """
        GetLinks method

        This method analizes the causality of the causet to create the Links dictionary, in wich each point in the causet
        is a key, with a value consisting of a set containing every direct future causal point to the key.
        """
        self.Links = GetLinks(self)
    
    def PrintCauset(self, directory: str = None) -> None: #Add save image in directory
        """
        PrintCauset method

        This method prints a 2D scatter plot of the causet spacetime diagram and saves it on a given directory.
        """
        PrintCauset(self.Causet, directory=directory)
    
    def PrintHasseDiagram(self, directory: str = None)  -> None: #Add save image in directory
        """
        PrintHasseDiagram method

        This method prints the Hasse Diagram of the Causet and saves it on a given directory
        """
        PrintHaseDiagram(self.Links, directory=directory)
    
    def PrintFuture(self, source: tuple[float], directory: str = None) -> None:
        """
        PrintFuture method

        This method prints a spacetime diagram of positions of the causet, alongside a color-coded causal relation with the 
        source point: if green, the point is causally related with the sorce; red otherwise.
        """
        future = ChronologicalFuture(source, self.Causet, self.Metric)
        PrintFuture(self.Links, future, directory=directory)
    
    def Geodesics(self, source: tuple[float], tarjet: tuple[float]) -> list:
        self.Geodesic = GetGeodesic(self.Links, source=source, target=tarjet)
    



#################################################################




class ContinuumSimulation():
    """
    Class ContinuumSimulation

    This object gathers all methods concerining with simulation and computing of geodesics in the continuum.

    Class attributes:
        Metric (MetricTensor Class): Describes the metric manifold to be considered in the causet sprinkling
        Source (2D float tuple): spacetime coordinates of the worldline starting point
        Tarjet (2D float tuple): spacetime coordinates of the worldline ending point
        Vr (float): spacial initial velocity
        Vt (float): temporal initial velocity
        tau_span (float list): proper time range for the integration method simulation
        type (str): Type of geodesic (timelike, null or lightlike, or spacelike)
        Geodesic([T,X] list, where X & Y are float lists): List describig spacetime coordinates along the geodesic 
    
    Class methods:
        ComputeGeodesic: Given a metric manifold and initial conditions, this function computes the corresponding geodesic.
        PrintGeodesic: Given a list of points describing a geodesic, this method prints the spacetime diagram of the trajectory.
    """
    def __init__(self, Metric: MetricTensor, source: tuple[float],
                  SpacialVelocity: float, tau_span: list[float], g_type: str = 'timelike') -> None:
        """
        Constructor for ContinuumSimulation class

        Parameters:
            Metric (MetricTensor Class): Describes the metric manifold to be considered in the causet sprinkling
            Source (2D float tuple): spacetime coordinates of the worldline starting point
            Tarjet (2D float tuple): spacetime coordinates of the worldline ending point
            Vr (float): spacial initial velocity
            Vt (float): temporal initial velocity
            tau_span (float list): proper time range for the integration method simulation
            type (str): Type of geodesic (timelike, null or lightlike, or spacelike)
        """
        self.Metric = Metric #MetricTensor class atribute
        self.source = source #Initial point for the geodesic
        self.tarjet = None #Final point of the geodesic (to be computed)
        self.Vr = SpacialVelocity
        self.Vt = ComputeVt(Metric=self.Metric, source=self.source, Vr=self.Vr, g_type=g_type) #Uses ComputeVt self method to obtain
        self.tau_span = tau_span #Proper time values to be considered
        self.type = g_type #Type of geodesic
        self.Geodesic = None #Points of the geodesic (to be computed)


        
    def ComputeGeodesic(self, TimeRange: tuple[float]=None, SpaceRange: tuple[float]=None):
        """
        ComputeGeodesic method

            This function computes all the points allong the geodesic for the given proper time parameters and initial conditions.
            The resulting set of points are cut to be withing a given temporal and spacial range to ensure consistency with rest of functions.

        Parameters: 
            TimeRange (2D float tuple): Lower and upper temporal bounds to be computed. (If non given this is not considered)
            SpaceRange (2D float tuple): Lower and upper spacial bounds to be computed. (If non given this is not considered)

        Returns:
            Does not return anything; results are automatically saved on self.Geodesic attribute.
        """

        #The points alongside the proper time range are computed using the ComputeGeodesic function from Continuum_Geodesics.py module.
        geodesic = ComputeGeodesic(Metric=self.Metric, source=self.source,
                               Vr=self.Vr, Vt=self.Vt, t_span=self.tau_span)
        
        #The final point of the geodesic is then considered to be the point tarjet.

        #If there are no temporal of spacial bounds, the final point is the tarjet
        if TimeRange is None and SpaceRange is None:
            self.Geodesic = geodesic
            self.tarjet = (geodesic[0][-1], geodesic[1][-1])
        
        # If there are temporal or spacial bounds, only the last point within said bounds is considered
        else:
            #To obtain the proper tarjet within bounds, the geodesic is cut using the CutGeodesic funcion from Continuum_Geodesics.py module.
            self.Geodesic = CutGeodesic(Geodesic=geodesic, TimeRange=TimeRange, SpaceRange=SpaceRange)
            self.tarjet = (self.Geodesic[0][-1], self.Geodesic[1][-1])
        
    def PrintGeodesic(self, directory: str = None): #Add save image in directory
        """
        PrinGeodesic method
            This function creates and saves a spacetime diagram of the geodesic.
        Parameters:
            Only the directory for saving the figure is needed, if none is given the figure will not be saved (only showed)
        """
        PrintContinuumGeodesic(Geodesic=self.Geodesic, directory=directory) 
