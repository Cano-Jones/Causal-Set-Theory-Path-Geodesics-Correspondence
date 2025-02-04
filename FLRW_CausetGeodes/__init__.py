"""
FLRW_CausetGeodes python Library

This python library, aims to study the path-geodesic corresponce in (1+1)FLRW metric manifolds.
A conjectured by Myrheim states that the most natural path in a causal set as correspondence to
a continuum geodesic, must be maximal in length. This conjecture has only been proved by 
Brightwell ans co-workers for a Minkowskian space. This library explores said correspondence for
de standard cosmological solution.

A script using this library must describe the metric tensor to be used (determining the scale
factor and space curvature), simulation parameters (incluiding point number, spacetime range,
spacetime subdivisions) & a continuum geodesic (can be computed using the provided functions).

A proper expansion of this project should include a Xi^2 study prooving the correspondence
between the continuum and the causet simulation; this is not provided in this library.

Author: Cano Jones, Alejandro
linkedin: www.linkedin.com/in/alejandro-cano-jones-5b20a7136
github: https://github.com/Cano-Jones
"""

from FLRW_CausetGeodes.Class_Objects import *
from FLRW_CausetGeodes.Printing_Module import PrintComparison