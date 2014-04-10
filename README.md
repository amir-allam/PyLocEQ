INTRODUCTION
-----------------------------------------------------------------
The eqloc3d package is a set of python codes to perform automated real-time earthquake location using complex 3D velocity models and no other assumptions. The software is optimized for use on real-time Antelope systems, although it can be used on non-Antelope systems with plain ascii input files. The process occurs in 3 steps:
1) Computation of traveltimes from each node to each station. This occurs only once for each station unless the velocity model is updated.
2) A 2-stage grid search to locate the grid node with the lowest residual arrival time
3) Computation of sub-grid location using least squares and location derivatives calculated by trilinear interpolation

Because nearly all of the significant computation occurs in step 1, steps 2 and 3 can occur in real time. Regular-grid binary files also facilitate the grid search in step 2.

1) Traveltime Computation
-----------------------------------------------------------------
Traveltimes are computed using the Fast Marching Method (FMM) software written by N. Rawlinson (de Kool et al., 2006; GJI). The FMM solves the eikonal traveltime equation using finite differences on a regular grid in spherical coordinates cut by interfaces (e.g., topography, moho). In principle, any traveltime solver can be used for this step, but the FMM has the following advantages:
- computes traveltimes at every grid point in a single simulation
- is stable even for extremely heterogeneous media
- can include known (or assumed) interfaces such as the moho
- includes topography for accurate station elevations and topographic effects on the wavefield
- spherical coordinates allow use on local, regional, or global scales

After parameterizing a velocity model (or more than one) on a regular structured grid, traveltimes are computed for each station by placing a source at the station location. Due to seismic reciprocity, this is equivalent to computing a traveltime from an event at each grid node in the volume to the receiver. These traveltimes are saved in simple binary float files with known byte order, which is the same for every station. Thus, the traveltime at grid point m starts at the 4*m byte in the traveltime files for each station. Because this number is the same for every file, traveltime lookup is rapid during step 2 below, and does not significantly increase with file size (or number of grid points).

After the traveltime files are computed, no new computations need to occur unless a new station is added or the velocity model is updated.

2) Two-stage Grid Search
-----------------------------------------------------------------
First, a coarse search is performed for every 8th node (by default) on the grid. For each arrival and for each potential grid node, a candidate origin time is computed by subtracting the traveltime from the observed arrival time. Misfit at each node is then calculated as the standard deviation of the computed origin time for all stations. The rationale is that the best-fit grid location will have a low standard deviation in origin time. By characterizing the misfit in this way, neither origin time nor P-S separation time need be assumed.

After the coarse minimum is located, a finer grid search is performed on every node within 10 nodes (by default; slightly larger than the coarse search spacing). Misfit is characterized in the same way using origin times. Because misfit formulated in this way is relatively smooth, the minimum found during the finer grid search is the global minimum.

3) Sub-grid Least Squares Estimate
-----------------------------------------------------------------
To solve for location at the sub-grid scale, an inverse problem is formulated:

d = G * m

Where d is a vector of observed arrival times, G are the derivatives with respect to x,y,z, and origin time with respect to the stations corresponding to d, and m are x,y,z, and origin time which are the unknowns. Using the best-fit grid node, the derivatives in G are computed by trilinear interpolation taking the differences with the surrounding nodes. The system is then solved using the SciPy linalg python module.


















