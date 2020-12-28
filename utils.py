import numpy as np
import skfmm
from scipy.signal import convolve2d


Ky = np.array([[1/2, 0., -1/2]])
Kx = Ky.T

def computeEikonal(xi, ti, grid, vel_field):
    """Computes the eikonal solution given the EASs xi and their timings ti.

    Parameters
    ----------
    xi : np.ndarray (int or float)
        [N, 2] array of float positions, or [N] array of DOF index, marking the EAS location
    ti : np.ndarray (float)
        [N] array, prescribing the initiation time of each EAS
    grid : np.ndarray (float)
        [K, L, 2] array describing the grid coordinates
    vel_field : np.ndarray (float)
        [K, L] array of velocities for each grid point

    Returns
    -------
    np.ndarray (float)
        [K, L] array of the computed activation times for each grid point
    """

    dx = grid[1, 0, 0] - grid[0, 0, 0]

    if np.issubdtype(xi.dtype, float):
        xi, ti = upsampleEAS(xi, ti, grid, vel_field)

    #The chosen eikonal solver can only handle one EAS at a time, so we have to compute the eikonal solution for each EAS separately
    #and then use the geodesic property to combine all solutions.
    #Note that the custom eikonal solver from the paper only needs to solve eikonal here once
    phi_total = np.ones(grid.shape[:-1]) * 1e5
    for i in range(ti.size):
        phi_current = np.ones_like(phi_total)
        phi_current[xi[i, 0], xi[i, 1]] = 0.
        phi_current = skfmm.travel_time(phi_current, speed=vel_field, order=1, dx=float(dx))
        phi_total = np.minimum(phi_total, phi_current + ti[i])

    return phi_total

def computeEikonalGradient(phi, grid, vel_field):
    """Computes the gradient -D \nabla \phi for a given phi

    Parameters
    ----------
    phi : np.ndarray (float)
        [K, L] array of activation times
    grid : np.ndarray (float)
        [K, L, 2] array describing the grid coordinates
    vel_field : np.ndarray (float)
        [K, L] array of velocities for each grid point

    Returns
    -------
    np.ndarray (float)
        [K, L, 2] array. Holds the gradient at each grid point
    """
    dx = grid[1, 0, 0] - grid[0, 0, 0]
    phi = np.pad(phi, [[1, 1], [1, 1]], mode='symmetric')
    nabla_phi = np.stack([convolve2d(phi[:, 1:-1], Kx/dx, mode='valid'),
             convolve2d(phi[1:-1], Ky/dx, mode='valid')], axis=-1)
    nabla_phi = -vel_field[..., np.newaxis]**2 * nabla_phi
    
    return nabla_phi

def getBoundingBoxes(pos, grid):
    """Computes the bounding box of each point in pos on the grid

    Parameters
    ----------
    pos : np.ndarray (float)
        [N, 2] array of position for which to compute the bounding boxes
    grid : np.ndarray (float)
        [K, L, 2] array describing the grid coordinates

    Returns
    -------
    tuple (np.ndarray (int), np.ndarray (float))
        Returns a tuple of two arrays:
        * [4, N, 2] array holding the 4 indices of the bounding box for each position in both axes
        * [4, N, 2] array holding the coordinates of the bounding box for each position
    """
    dx = grid[1, 0, 0] - grid[0, 0, 0]
    dims = grid.shape[0]
    lower_left_corners = (pos // dx).astype(np.int32) #Indices

    #Special case for elements sitting exactly on the top or right boundary
    lower_left_corners = np.where(lower_left_corners == dims - 1, lower_left_corners - 1, lower_left_corners)
    lr_corner = lower_left_corners.copy()
    lr_corner[..., 0] += 1
    ul_corner = lower_left_corners.copy()
    ul_corner[..., 1] += 1
    ur_corner = ul_corner.copy()
    ur_corner[..., 0] += 1
    bbox_i = np.stack([lower_left_corners, 
                    lr_corner, ur_corner, ul_corner], axis=0)

    return bbox_i, grid[bbox_i[..., 0], bbox_i[..., 1]]


def upsampleEAS(xi, ti, grid, vel_field):
    """Computes the discrete indices of the EAS for continuously given xi. Note that this requirement comes from
    the implementation of many solvers, necessitating the Dirichlet boundary conditions to coincide with a DOF.

    Parameters
    ----------
    xi : np.ndarray (float)
        [N, 2] array holding the EAS positions
    ti : np.ndarray (float)
        [N] array holding the EAS timings
    grid : np.ndarray (float)
        [K, L, 2] array describing the grid coordinates
    vel_field : np.ndarray (float)
        [K, L] array of velocities for each grid point

    Returns
    -------
    tuple (np.ndarray (int), np.ndarray (float))
        Returns a tuple of two arrays
        * The bounding box indices of each EAS
        * The timings on the bounding boxes. The timings account for the time it takes to travel from xi to bbox(xi) with linear assumptions.
    """
    eas_indices, bbox = getBoundingBoxes(xi, grid)
    vel_xi = vel_field[eas_indices[..., 0], eas_indices[..., 1]] #Each grid point of the bounding box takes its own velocity
    eas_ti = ti + 1/vel_xi * np.linalg.norm(xi[np.newaxis] - bbox, axis=-1)

    eas_indices = eas_indices.reshape([-1, 2])
    eas_ti = eas_ti.reshape([-1])

    #TODO: Duplicates
    flat_indices = np.ravel_multi_index(eas_indices.T, grid.shape[:-1])
    assert(flat_indices.size == np.unique(flat_indices.size))

    return eas_indices, eas_ti


def quadUpsampling(x, grid, vals):
    """Computes the values at each position x on the grid, using a bi-linear upsampling method, or equivalently,
    a quadrilateral element assumption.

    Parameters
    ----------
    x : np.ndarray (float)
        [N, 2] array of positions where to evaluate the function defined by vals
    grid : np.ndarray (float)
        [K, L, 2] array describing the grid coordinates
    vals : np.ndarray (float)
        [K, L] array of fixed values at the DOFs/grid coordinates

    Returns
    -------
    np.ndarray
        [N] array of vals evaluated at x, using bi-linear upsampling.
    """
    dx = grid[1, 0, 0] - grid[0, 0, 0]
    bbox_i, bbox = getBoundingBoxes(x, grid)
    vals_bbox = vals[bbox_i[..., 0], bbox_i[..., 1]]
    ll_corners = bbox[0]
    x_local = (x - ll_corners) / dx

    #Compute the basis functions at the current points
    phi1 = (1 - x_local[..., 0]) * (1 - x_local[..., 1])
    phi2 = x_local[..., 0] * (1 - x_local[..., 1])
    phi3 = x_local[..., 0] * x_local[..., 1]
    phi4 = (1 - x_local[..., 0]) * x_local[..., 1]

    phis = np.stack([phi1, phi2, phi3, phi4], axis=0)

    return np.sum(phis[..., np.newaxis] * vals_bbox, axis=0)

def computeGeodesics(nabla_phi, grid, start_points, proj_op, alpha=1e-2, max_iterations=1000):
    """Computes the geodesics by solving a simple geodesic flow ODE using the midpoint method. For more information on the actual
    ODE that is solved, we refer to the arXiv paper (see README.md). Note that the ODE will stop if all geodesics have converged,
    which is defined here by a change smaller than 1e-10 between two ODE iterations.

    Parameters
    ----------
    nabla_phi : np.ndarray (float)
        [K, L, 2] array holding -D \nabla \phi (see computeEikonalGradient)
    grid : np.ndarray (float)
        [K, L, 2] array describing the grid coordinates
    start_points : np.ndarray (float)
        [N, 2] array of start points/initial values (\gamma(0)) where the ODE will start
    proj_op : function/lambda
        A function that can project all ODE solutions in each step back onto the grid
    alpha : float, optional
        Step size for the ODE solver, by default 1e-2
    max_iterations : int, optional
        Maximum amount of iterations for each ODE before quitting, by default 1000

    Returns
    -------
    np.ndarray (float)
        [max_iterations, N, 2] array holding all geodesics. If the algorithm converges before max_iterations
        were reached, the remaining iterations will be assigned the last, converged solution.
    """
    geodesics = np.zeros(shape=[max_iterations] + list(start_points.shape))    

    current_pos = start_points
    #In contrast to the paper, we do not shuffle the converged geodesics here
    for i in range(max_iterations):
        geodesics[i] = current_pos

        nabla_phi_x = quadUpsampling(current_pos, grid, nabla_phi)
        current_pos_h_2 = current_pos + alpha/2 * nabla_phi_x
        nabla_phi_x_h_2 = quadUpsampling(current_pos_h_2, grid, nabla_phi)
        current_pos = proj_op(current_pos + alpha * nabla_phi_x_h_2)

        if i > 1 and np.all(np.linalg.norm(geodesics[i] - geodesics[i-1], axis=-1) < 1e-10):
            geodesics[i+1:] = geodesics[i:i+1]
            break

    return geodesics

def compSourceGradients(geodesics, xi, grid, nabla_phi, vel_field, neighborhood_eps):
    """Computes the 'source gradients', which are the gradients of the activation times $\phi$ at the points
    start_points in the function computeGeodesics. With this, the Jacobian matrix is built to approximate the influence
    of changing xi/ti on start_points. For more information on the mathematical background of this function, 
    we refer to the arXiv paper (see README.md). Note that if a geodesic does not reach any xi within neighborhood_eps,
    it is marked as invalid and its corresponding jacobian entry is set to zero. This can be a result of numerical
    approximation errors (e.g. neighborhood_eps chosen too small), or start_points and xi coinciding.

    Parameters
    ----------
    geodesics : np.ndarray (float)
        [max_iterations, N, 2] array holding all geodesics.
    xi : np.ndarray (float)
        [M, 2] array holding the EAS positions
    grid : np.ndarray (float)
        [K, L, 2] array describing the grid coordinates
    nabla_phi : np.ndarray (float)
        [K, L, 2] array holding -D \nabla \phi (see computeEikonalGradient)
    vel_field : np.ndarray (float)
        [K, L] array of velocities for each grid point
    neighborhood_eps : float
        Tells the function in which neighborhood the gradient should be calculated, since the actual EAS is not always hit.
        Should be optimally smaller than the grid spacing.

    Returns
    -------
    tuple (np.ndarray (float), np.ndarray (float))
        Returns a tuple of
        * \dot{\gamma}(t) at the end of the optimization
        * [N, 3*M] array holding the Jacobian matrix. The columns hold the flattened parameter influences. [:, 0:2*M] are reserved for xi,
          whereas [:, 2*M:] are for the corresponding ti.
    """
    nr_xi = xi.shape[0]
    nr_geodesics = geodesics.shape[1]

    #Compute all distances
    xi_dists = np.linalg.norm(geodesics[..., np.newaxis, :] - xi[np.newaxis, np.newaxis], axis=-1)

    #Neighborhood reached?
    xi_reached = np.any(xi_dists < neighborhood_eps, axis=0)
    source_geodesic_valid = np.any(xi_reached, axis=1)

    #Compute to which xi the geodesic converges to
    xi_ind = np.argmax(xi_reached, axis=1)

    #Find the earliest ODE iteration in which the geodesic reaches the neighborhood
    converg_iter = np.argmax(xi_dists[:, np.arange(nr_geodesics), xi_ind] < neighborhood_eps, axis=0)
    converg_pos = geodesics[converg_iter, np.arange(nr_geodesics)]

    #Evaluate the gradient at this position (already weighted with velocities)
    vel_xi = quadUpsampling(converg_pos, grid, vel_field[..., np.newaxis])
    gamma_dot_converg = 1/vel_xi**2 * quadUpsampling(converg_pos, grid, nabla_phi)
    gamma_dot_converg = gamma_dot_converg / (1/vel_xi * np.linalg.norm(gamma_dot_converg, axis=-1, keepdims=True))

    jacobian = np.zeros(shape=[nr_geodesics, nr_xi*3])
    jacobian[np.arange(nr_geodesics), 2*xi_ind] = gamma_dot_converg[..., 0]
    jacobian[np.arange(nr_geodesics), 2*xi_ind+1] = gamma_dot_converg[..., 1]

    #Influence of the EAS timing
    jacobian[np.arange(nr_geodesics), 2*nr_xi + xi_ind] = 1.

    #Remove invalid source geodesics
    jacobian[~source_geodesic_valid, :] = 0.

    return gamma_dot_converg, jacobian

