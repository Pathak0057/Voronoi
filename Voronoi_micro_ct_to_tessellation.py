"""
==============================================================================
 Micro_ct_abaqus.py
 -------------------
 PURPOSE
 -------
 This script converts Micro-CT fiber centroid data into a mesoscale Voronoi
 tessellation of a composite Representative Volume Element (RVE).

 The physical motivation is that short-fiber composites have locally
 non-uniform fiber distributions. Different sub-regions of the same
 specimen have different fiber volume fractions and orientation statistics.
 A single homogeneous material model cannot capture this scatter.

 Instead, we divide the RVE into discrete Voronoi "cells", each of which
 represents a local neighbourhood of fibers. Each cell will later receive
 its own set of material properties (computed from its local fiber statistics)
 before being assembled into the full mesoscale FE model in Abaqus.

 ALGORITHM SUMMARY
 -----------------
 1. Load fiber line-segment data from the Excel file
    (each fiber is defined by its two 3D endpoints: ep0 and ep1).
 2. Place a regular 3D grid of sphere centers across the domain.
 3. For each sphere center:
      a. Find all fiber segments that pass through (intersect) the sphere.
      b. Compute the mid-point of each intersecting fiber segment inside
         the sphere.
      c. Average those mid-points (plus the sphere center itself) to get
         a "physics-informed" Voronoi seed — the seed is biased toward
         where fibers actually are, not just the geometric center.
 4. Remove duplicate seeds.
 5. For each seed, compute its Voronoi cell by intersecting perpendicular
    bisector half-spaces from neighbouring seeds with the bounding cube.
    The intersection is found by exhaustive triple-plane vertex enumeration
    followed by a ConvexHull to recover the cell faces.
 6. Export each Voronoi cell as an STL triangular mesh.
 7. For each cell, clip all fiber segments to the cell volume and compute:
      - Cell volume (convex hull)
      - Total fiber volume inside
      - Local fiber volume fraction
    Export a summary CSV.
 8. Produce four diagnostic figures and an animated GIF of the sphere sliding.

 INPUTS
 ------
 final_fiber_data_CentrelinePoints_Abaqus_static.xlsx
   Columns used: ep0_x, ep0_y, ep0_z  (fiber start point)
                 ep1_x, ep1_y, ep1_z  (fiber end point)
   Each row is one fiber segment measured from Micro-CT data.

 OUTPUTS
 -------
 VoronoiSTL/region_001.stl ... region_NNN.stl   -- per-cell geometry
 voronoi_fiber_summary.csv                       -- per-cell fiber statistics
 Voronoi_model_sphere_sliding.gif                -- visualization of seeding
 4 Matplotlib figure windows                     -- diagnostic plots

 DEPENDENCIES
 ------------
 pip install numpy pandas scipy matplotlib numpy-stl pillow scikit-learn

 HOW TO RUN
 ----------
 python Micro_ct_abaqus.py
 (Must be run from the folder containing the Excel file, or update excel_file path)

==============================================================================
"""

import os
import math
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')          # change to 'Qt5Agg' / 'Agg' if needed
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist
from sklearn.neighbors import KDTree
from stl import mesh as stl_mesh
from PIL import Image
import io

# =============================================================================
# PARAMETERS
# =============================================================================

# --- Input file ---
# The Excel file exported from Micro-CT image analysis.
# Each row is one fiber. Columns ep0_x/y/z and ep1_x/y/z are the two
# 3D endpoints of the fiber centroidal axis segment.
excel_file = 'final_fiber_data_CentrelinePoints_Abaqus_static.xlsx'

# --- Sphere sliding grid ---
# A uniform 3D grid of sphere centers is placed between grid_min and grid_max
# with spacing grid_step. Each sphere defines the "local neighbourhood" used
# to compute one Voronoi seed.
# Rule of thumb: grid_step should be roughly equal to or slightly smaller than
# sphere_radius so that spheres overlap and every part of the RVE is sampled.
grid_min  = -70
grid_max  =  70
grid_step =  50   # spacing between sphere centers (same units as fiber coords)

# --- Sphere radius ---
# The sampling sphere radius. Fibers whose line segments intersect a sphere
# of this radius centered at the current grid point are included in that
# sphere's local fiber cluster.
# Larger radius → more fibers captured per sphere → smoother seeds but less
# sensitivity to local heterogeneity. Smaller → more local but risk of empty spheres.
sphere_radius = 50

# --- RVE bounding cube ---
# The cube defines the overall RVE extent. All Voronoi cells will be clipped
# to lie within [cube_min, cube_max] in each dimension.
# This should match the physical RVE size inferred from Micro-CT.
cube_min = -50
cube_max =  50

# --- Neighbour search radius for Voronoi construction ---
# When building the half-space intersection for cell i, we only consider
# seeds j whose distance to seed i is <= neighbor_radius.
# Setting this too small can miss far-away seeds and produce unbounded cells;
# setting too large increases computation time.
neighbor_radius = 80

# --- Output ---
output_folder = 'VoronoiSTL'     # folder where STL files are written
os.makedirs(output_folder, exist_ok=True)
output_csv = 'voronoi_fiber_summary.csv'   # per-region statistics CSV

# --- Numerical tolerances ---
tol_det  = 1e-10   # determinant threshold for near-singular 3x3 systems
                   # (avoids solving degenerate plane intersections)
tol_ineq = 1e-7    # tolerance for the half-space feasibility check
                   # (a vertex is "inside" if A*x <= b + tol_ineq for all planes)
snap_tol = 1e-6    # rounding tolerance for vertex deduplication
                   # (vertices closer than snap_tol are treated as identical)

# --- Figure appearance ---
face_alpha      = 0.5      # transparency of Voronoi cell faces in 3D plots
show_edges      = False    # draw cell edges (can clutter dense tessellations)
max_plot_cells  = np.inf   # limit how many cells are drawn (use int to subsample)

fiber_line_width    = 2.0
clip_fibers_to_cube = True   # clip fiber lines to the cube extents before plotting

fiber_color_mode    = "byCell"   # "byCell" = fibers colored by which Voronoi cell
                                  # they belong to; "uniform" = single color
fiber_uniform_color = np.array([0.85, 0.1, 0.1])

# --- 3D cylinder visualization of fibers ---
fiber_vis_radius = 1.2    # visual radius of fiber cylinders in Figure 1 (not physical)
fiber_vis_n_circ = 12     # number of circumferential facets on each cylinder

# --- Animated GIF (sphere sliding visualization) ---
enable_gif_generation = True
gif_filename          = 'Voronoi_model_sphere_sliding.gif'
gif_delay_time        = 0.15    # seconds per frame
gif_max_frames        = np.inf  # limit total frames (use int to shorten GIF)

sphere_rgba = np.array([0.0, 0.0, 1.0, 0.25])   # blue semi-transparent sphere
fiber_rgba  = np.array([1.0, 0.0, 0.0, 1.0])    # solid red fibers in GIF

# --- Fiber volume calculation ---
# r_fiber is the physical cross-sectional radius of a fiber (same units as coords).
# Used to compute fiber volume = pi * r_fiber^2 * length_inside_cell,
# and hence the local fiber volume fraction for each Voronoi cell.
r_fiber = 4.0

# Small tolerance for half-space plane clipping (fiber-in-cell detection)
inside_tol_plane = 1e-9


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def sphere_segment_intersection(sc, R, A, B):
    """
    Test whether the line segment A→B intersects a sphere centered at sc
    with radius R.

    Uses the standard ray-sphere quadratic formula. The segment parameter
    t is clamped to [0, 1] so only the portion between A and B is tested.

    Returns
    -------
    hit : bool   – True if the segment crosses the sphere
    t1  : float  – parameter of entry point along A→B (clamped to [0,1])
    t2  : float  – parameter of exit point along A→B (clamped to [0,1])

    Why we need this: each sphere center represents a "local window" in the
    composite. We want to find all fibers that pass through that window so we
    can compute a local fiber centroid and thus a physics-informed Voronoi seed.
    """
    AB = B - A
    a  = np.dot(AB, AB)
    if a < 1e-14:
        return False, 0.0, 0.0        # degenerate zero-length segment
    AS = A - sc
    b  = 2.0 * np.dot(AB, AS)
    c  = np.dot(AS, AS) - R**2
    D  = b*b - 4*a*c
    if D < 0:
        return False, 0.0, 0.0        # no real intersection (misses sphere)
    sd = math.sqrt(D)
    u1 = (-b - sd) / (2*a)
    u2 = (-b + sd) / (2*a)
    t1 = max(0.0, min(1.0, min(u1, u2)))
    t2 = max(0.0, min(1.0, max(u1, u2)))
    hit = t1 < t2
    return hit, t1, t2


def compute_seed(sc, R, P1s, P2s):
    """
    Compute a Voronoi seed for a sphere centered at sc with radius R.

    For each fiber (P1→P2) that intersects the sphere, we find the mid-point
    of the segment portion inside the sphere. All these mid-points, together
    with the sphere center itself, are averaged to produce the seed.

    Why this instead of just using the grid center?
    -----------------------------------------------
    If we used the raw grid centers as seeds, the resulting Voronoi cells would
    be purely geometric and would not reflect the actual fiber distribution.
    By biasing each seed toward where fibers are concentrated, the resulting
    Voronoi cells naturally align with local fiber clusters. Cells with more
    fibers become smaller (higher local fiber density → seed pulled inward),
    which is physically meaningful.

    If no fibers are found inside the sphere, the sphere center is returned
    unchanged (fall-back to geometric seed).
    """
    fiber_centroids = []
    for j in range(len(P1s)):
        hit, t1, t2 = sphere_segment_intersection(sc, R, P1s[j], P2s[j])
        if hit:
            A = P1s[j] + t1 * (P2s[j] - P1s[j])
            B = P1s[j] + t2 * (P2s[j] - P1s[j])
            fiber_centroids.append((A + B) / 2.0)
    if len(fiber_centroids) == 0:
        return sc.copy()
    pts = np.vstack(fiber_centroids + [sc])
    return pts.mean(axis=0)


def clip_segment_to_aabb(A, B, bmin, bmax):
    """
    Clip a 3D line segment A→B to the axis-aligned bounding box
    [bmin, bmax]^3 using the parametric Liang-Barsky algorithm.

    Returns (ok, A_clipped, B_clipped).
    ok is False if the segment lies entirely outside the box.

    Used to: (1) avoid drawing fibers outside the RVE cube in figures,
             (2) clip fiber segments to the cube before computing their
                 length inside each Voronoi cell.
    """
    d  = B - A
    t0 = 0.0
    t1 = 1.0
    for axis in range(3):
        if abs(d[axis]) < 1e-14:
            if A[axis] < bmin or A[axis] > bmax:
                return False, A.copy(), B.copy()
        else:
            inv_d = 1.0 / d[axis]
            t_near = (bmin - A[axis]) * inv_d
            t_far  = (bmax - A[axis]) * inv_d
            if t_near > t_far:
                t_near, t_far = t_far, t_near
            t0 = max(t0, t_near)
            t1 = min(t1, t_far)
            if t0 > t1:
                return False, A.copy(), B.copy()
    A2 = A + t0 * d
    B2 = A + t1 * d
    return True, A2, B2


def build_convex_halfspaces_from_faces(V, F, centroid):
    """
    Given a convex polyhedron defined by vertices V and triangular faces F,
    build a set of inward-facing half-spaces [n | d] where n·x <= d defines
    the interior side of each face.

    The inward direction is enforced by checking that the centroid satisfies
    n·centroid <= d; if not, the normal is flipped.

    Why we need this: to determine whether a fiber segment passes through a
    specific Voronoi cell, we represent the cell as a set of linear inequalities
    and clip the segment against them. This is more general and exact than
    bounding-box clipping.
    """
    planes = []
    for tri in F:
        p1 = V[tri[0]]
        p2 = V[tri[1]]
        p3 = V[tri[2]]
        n  = np.cross(p2 - p1, p3 - p1)
        nn = np.linalg.norm(n)
        if nn < 1e-14:
            continue                   # degenerate face, skip
        n = n / nn
        d = np.dot(n, p1)
        if np.dot(n, centroid) > d:    # centroid is on the wrong side — flip
            n = -n
            d = -d
        planes.append(np.append(n, d))
    return np.array(planes) if planes else np.zeros((0, 4))


def clip_segment_with_convex_halfspaces(A, B, planes, tol):
    """
    Clip the segment A→B against a convex polytope defined by half-spaces.
    Each row of `planes` is [nx, ny, nz, d] where n·x <= d is the interior.

    Uses the Cyrus-Beck algorithm: for each plane, compute the parametric
    intersection of the segment with the plane and update the [t0, t1] interval.

    Returns (hit, A_clip, B_clip, length_inside).
    hit = True if any part of the segment lies inside the polytope.

    Used to determine what fraction of each fiber lies inside each Voronoi cell,
    which is needed for computing local fiber volume fractions.
    """
    dvec = B - A
    t0   = 0.0
    t1   = 1.0
    for row in planes:
        n = row[:3]
        d = row[3]
        num = d - np.dot(n, A)
        den = np.dot(n, dvec)
        if abs(den) < 1e-14:
            if np.dot(n, A) > d + tol:
                return False, A.copy(), B.copy(), 0.0   # parallel and outside
            continue
        t = num / den
        if den > 0:
            t1 = min(t1, t)   # entering constraint
        else:
            t0 = max(t0, t)   # exiting constraint
        if t0 > t1:
            return False, A.copy(), B.copy(), 0.0   # clipped to empty
    t0 = max(0.0, min(1.0, t0))
    t1 = max(0.0, min(1.0, t1))
    if t0 >= t1:
        return False, A.copy(), B.copy(), 0.0
    Aclip      = A + t0 * dvec
    Bclip      = A + t1 * dvec
    len_inside = np.linalg.norm(Bclip - Aclip)
    hit        = len_inside > 0
    return hit, Aclip, Bclip, len_inside


# ---------- Drawing helpers --------------------------------------------------
# These functions produce 3D visualizations for diagnostic figures.
# They do not affect any computed results.

def draw_cube_wireframe(ax, cmin, cmax):
    """Draw the 12 edges of the RVE bounding cube as black lines."""
    P = np.array([
        [cmin, cmin, cmin], [cmax, cmin, cmin],
        [cmax, cmax, cmin], [cmin, cmax, cmin],
        [cmin, cmin, cmax], [cmax, cmin, cmax],
        [cmax, cmax, cmax], [cmin, cmax, cmax]])
    edges = [(0,1),(1,2),(2,3),(3,0),
             (4,5),(5,6),(6,7),(7,4),
             (0,4),(1,5),(2,6),(3,7)]
    for e in edges:
        xs = [P[e[0],0], P[e[1],0]]
        ys = [P[e[0],1], P[e[1],1]]
        zs = [P[e[0],2], P[e[1],2]]
        ax.plot(xs, ys, zs, 'k-', linewidth=1.2)


def make_cylinder_mesh(p1, p2, radius, n_circ):
    """
    Build a parametric cylinder surface mesh between points p1 and p2.
    Returns (X, Y, Z) grid arrays suitable for ax.plot_surface().

    We use a local coordinate frame aligned with the fiber axis (p2-p1).
    Two orthogonal radial directions (u, w) are computed via cross-product
    with a reference vector. Circumferential points are rotated around these
    two directions to form the cylinder surface.

    This is for visualization only (Figure 1 and the GIF).
    """
    v    = p2 - p1
    L    = np.linalg.norm(v)
    if L < 1e-12:
        return None
    vhat = v / L
    ref  = np.array([0,0,1]) if abs(vhat[2]) < 0.9 else np.array([0,1,0])
    u    = np.cross(vhat, ref); u /= np.linalg.norm(u)
    w    = np.cross(vhat, u)

    theta  = np.linspace(0, 2*np.pi, n_circ)
    z_vals = np.array([0.0, L])
    tt, zz = np.meshgrid(theta, z_vals)

    X_loc = radius * np.cos(tt)
    Y_loc = radius * np.sin(tt)

    Px = p1[0] + u[0]*X_loc + w[0]*Y_loc + vhat[0]*zz
    Py = p1[1] + u[1]*X_loc + w[1]*Y_loc + vhat[1]*zz
    Pz = p1[2] + u[2]*X_loc + w[2]*Y_loc + vhat[2]*zz
    return Px, Py, Pz


def draw_cylinder_segment(ax, p1, p2, radius, n_circ, rgba):
    """Render a fiber as a 3D cylinder surface on axes ax."""
    result = make_cylinder_mesh(p1, p2, radius, n_circ)
    if result is None:
        return
    Px, Py, Pz = result
    color = rgba[:3]
    alpha = rgba[3] if len(rgba) > 3 else 1.0
    ax.plot_surface(Px, Py, Pz,
                color=color,
                alpha=alpha,
                linewidth=0,
                antialiased=True,
                shade=False)


def draw_simple_sphere(ax, center, radius, rgba):
    """
    Render a translucent sphere at a given center.
    Used in the GIF animation to show the sliding sampling sphere.
    Returns the surface handle so it can be removed between frames.
    """
    u_ang = np.linspace(0, 2*np.pi, 25)
    v_ang = np.linspace(0, np.pi, 25)
    xs = radius * np.outer(np.cos(u_ang), np.sin(v_ang)) + center[0]
    ys = radius * np.outer(np.sin(u_ang), np.sin(v_ang)) + center[1]
    zs = radius * np.outer(np.ones(25),   np.cos(v_ang)) + center[2]
    color = rgba[:3]; alpha = rgba[3] if len(rgba) > 3 else 0.25
    return ax.plot_surface(xs, ys, zs, color=color, alpha=alpha,
                           linewidth=0, antialiased=True)


def draw_cube_faces(ax, cmin, cmax, face_color, face_alpha_val):
    """Draw the 6 faces of the RVE cube as semi-transparent polygons."""
    verts = np.array([
        [cmin,cmin,cmin],[cmax,cmin,cmin],[cmax,cmax,cmin],[cmin,cmax,cmin],
        [cmin,cmin,cmax],[cmax,cmin,cmax],[cmax,cmax,cmax],[cmin,cmax,cmax]])
    faces_idx = [
        [0,1,2,3],[4,5,6,7],[0,1,5,4],
        [1,2,6,5],[2,3,7,6],[3,0,4,7]]
    polys = [[verts[fi] for fi in face] for face in faces_idx]
    pc = Poly3DCollection(polys,
        facecolor=list(face_color)+[face_alpha_val],
        edgecolor='k', linewidth=1.0)
    ax.add_collection3d(pc)


# =============================================================================
# LOAD FIBER DATA
# =============================================================================
# Read the Micro-CT centroidal axis data.
# Each row is a fiber represented as a line segment between two 3D endpoints.
# P1s[i] = start point of fiber i, P2s[i] = end point of fiber i.
tbl      = pd.read_excel(excel_file)
splines  = tbl['Spline_name'].unique()
P1s      = tbl[['ep0_x','ep0_y','ep0_z']].values.astype(float)
P2s      = tbl[['ep1_x','ep1_y','ep1_z']].values.astype(float)
n_fibers = len(P1s)
print(f'Loaded {n_fibers} fibers.')

# =============================================================================
# BUILD SEED GRID
# =============================================================================
# Create a regular 3D grid of sphere center positions.
# np.meshgrid with indexing='ij' gives a grid in (x,y,z) order consistent
# with our coordinate convention.
steps = np.arange(grid_min, grid_max + grid_step, grid_step)
gx, gy, gz = np.meshgrid(steps, steps, steps, indexing='ij')
sphere_centers_all = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])
n_seeds_grid = len(sphere_centers_all)

# For each grid point, compute the physics-informed Voronoi seed.
# compute_seed() returns the average position of all fiber midpoints inside
# the sphere, biased toward regions of high fiber density.
raw_seeds = np.array([
    compute_seed(sphere_centers_all[i], sphere_radius, P1s, P2s)
    for i in range(n_seeds_grid)
])

# Remove duplicate seeds (multiple grid points may produce identical seeds
# if they sample the same fiber cluster).
# np.unique on axis=0 returns sorted unique rows; we restore original order
# via argsort of the index array to keep seeds in a stable spatial order.
_, ia = np.unique(raw_seeds, axis=0, return_index=True)
ia = np.sort(ia)
seeds          = raw_seeds[ia]
sphere_centers = sphere_centers_all[ia]
n_seeds        = len(seeds)
print(f'Seeds after deduplication: {n_seeds}')

# =============================================================================
# PRECOMPUTE SEED DISTANCES
# =============================================================================
# D[i,j] = Euclidean distance between seeds i and j.
# Precomputing the full distance matrix avoids recomputing distances inside
# the Voronoi cell construction loop, which would be O(n^3).
D = cdist(seeds, seeds)

# =============================================================================
# HALF-SPACES FOR CUBE
# =============================================================================
# The RVE bounding cube is represented as 6 half-space inequalities:
#   x <= cube_max,  -x <= -cube_min  (i.e. x >= cube_min)
#   y <= cube_max,  -y <= -cube_min
#   z <= cube_max,  -z <= -cube_min
# These are included in every cell's half-space set to clip cells to the cube.
A_cube = np.array([
    [ 1, 0, 0], [-1, 0, 0],
    [ 0, 1, 0], [ 0,-1, 0],
    [ 0, 0, 1], [ 0, 0,-1]], dtype=float)
b_cube = np.array([cube_max, -cube_min,
                   cube_max, -cube_min,
                   cube_max, -cube_min], dtype=float)

# =============================================================================
# BUILD VORONOI CELLS
# =============================================================================
# For each seed i:
#   1. Collect all neighbouring seeds j within neighbor_radius.
#   2. For each neighbour j, add the perpendicular bisector half-space:
#        (sj - si) · x <= 0.5*(|sj|^2 - |si|^2)
#      This is the set of points closer to si than to sj — the definition
#      of the Voronoi cell boundary between i and j.
#   3. Combine these with the 6 cube half-spaces.
#   4. Find all vertices of the resulting polytope by solving every possible
#      combination of 3 planes simultaneously (triple-plane intersection).
#      Each solution is a candidate vertex; it is valid only if it satisfies
#      ALL other half-space inequalities.
#   5. Round and deduplicate vertices (floating-point noise can produce
#      near-identical vertices from different plane triples).
#   6. Compute the ConvexHull of valid vertices to get the cell faces.
#   7. Export the cell as an STL triangular mesh.
#
# WHY triple-plane enumeration instead of scipy.spatial.Voronoi?
# scipy.Voronoi only works for unbounded Voronoi diagrams and does not
# support clipping to an arbitrary convex domain. The half-space intersection
# approach directly produces the bounded, cube-clipped Voronoi cells we need.

all_vertices = [None] * n_seeds
all_faces    = [None] * n_seeds

for i in range(n_seeds):
    si        = seeds[i]
    neigh_idx = np.where((D[i] > 0) & (D[i] <= neighbor_radius))[0]

    if len(neigh_idx) == 0 and np.any((si < cube_min) | (si > cube_max)):
        continue   # seed is outside the cube and has no neighbours → skip

    # Start with the 6 cube planes
    A_planes = A_cube.copy()
    b_planes = b_cube.copy()

    # Add one perpendicular bisector half-space per neighbour
    for j in neigh_idx:
        sj   = seeds[j]
        n_ij = sj - si                          # normal points from i to j
        c_ij = 0.5 * (np.dot(sj, sj) - np.dot(si, si))  # midplane offset
        A_planes = np.vstack([A_planes, n_ij])
        b_planes = np.append(b_planes, c_ij)

    m_planes   = len(b_planes)
    verts_list = []

    # Enumerate all triples of planes and find their intersection point
    for p in range(m_planes - 2):
        for q in range(p+1, m_planes - 1):
            for r in range(q+1, m_planes):
                M = np.array([A_planes[p], A_planes[q], A_planes[r]])
                if abs(np.linalg.det(M)) < tol_det:
                    continue   # nearly coplanar planes → no unique intersection
                rhs = np.array([b_planes[p], b_planes[q], b_planes[r]])
                try:
                    x = np.linalg.solve(M, rhs)
                except np.linalg.LinAlgError:
                    continue
                # Accept only vertices that satisfy all constraints
                if np.all(A_planes @ x <= b_planes + tol_ineq):
                    verts_list.append(x)

    if not verts_list:
        continue

    verts = np.array(verts_list)
    verts = np.round(verts / snap_tol) * snap_tol   # snap to grid to merge near-duplicates
    verts = np.unique(verts, axis=0)                # remove exact duplicates
    if len(verts) < 4:
        continue   # need at least 4 points to form a 3D convex hull

    try:
        hull = ConvexHull(verts)
    except Exception:
        continue   # degenerate (co-planar) vertex set → skip

    all_vertices[i] = verts
    all_faces[i]    = hull.simplices

    # --- Export this Voronoi cell as an STL ---
    # We write a binary STL with one triangle per ConvexHull simplex.
    # These STLs are the input to Step 2 (FreeCAD shape healing).
    n_tri     = len(hull.simplices)
    cell_mesh = stl_mesh.Mesh(np.zeros(n_tri, dtype=stl_mesh.Mesh.dtype))
    for k, tri in enumerate(hull.simplices):
        for dim in range(3):
            cell_mesh.vectors[k][dim] = verts[tri[dim]]
    stl_name = os.path.join(output_folder, f'region_{i+1:03d}.stl')
    cell_mesh.save(stl_name)

print('Done. Saved geometry and STL files.')

# =============================================================================
# CSV SUMMARY
# =============================================================================
# For each Voronoi cell, compute:
#   - Cell volume via ConvexHull.volume
#   - For each fiber: clip to cell half-spaces, get length inside cell
#   - Total fiber volume = sum(pi * r_fiber^2 * length_inside)
#   - Fiber volume fraction = total_fiber_volume / cell_volume
# Write all statistics to a CSV that can be used to:
#   (a) Validate the tessellation quality
#   (b) Compute per-cell material properties for Abaqus (micromechanics)

region_ids          = list(range(1, n_seeds+1))
voronoi_volumes     = np.zeros(n_seeds)
total_fiber_volumes = np.zeros(n_seeds)
fiber_volume_fracs  = np.zeros(n_seeds)
fiber_indices_col   = [''] * n_seeds
fiber_lengths_col   = [''] * n_seeds
fiber_p1_col        = [''] * n_seeds
fiber_p2_col        = [''] * n_seeds
min_corner_col      = [''] * n_seeds
cuboid_size_col     = [''] * n_seeds

print('\nComputing per-region fiber summary for CSV...')

for i in range(n_seeds):
    V = all_vertices[i]
    F = all_faces[i]
    if V is None or F is None:
        continue

    # Axis-aligned bounding box of the cell (for reference)
    mn = V.min(axis=0)
    mx = V.max(axis=0)
    min_corner_col[i]  = f'[{mn[0]:.3f} {mn[1]:.3f} {mn[2]:.3f}]'
    cuboid_size_col[i] = f'[{mx[0]-mn[0]:.3f} {mx[1]-mn[1]:.3f} {mx[2]-mn[2]:.3f}]'

    try:
        hull = ConvexHull(V)
        vol  = hull.volume
    except Exception:
        vol = 0.0
    voronoi_volumes[i] = vol

    # Build half-space representation of this cell for fiber clipping
    centroid = V.mean(axis=0)
    planes   = build_convex_halfspaces_from_faces(V, F, centroid)

    f_ids     = []
    f_lengths = []
    fp1s      = []
    fp2s      = []
    total_vol = 0.0

    for j in range(n_fibers):
        hit, Aclip, Bclip, len_in = clip_segment_with_convex_halfspaces(
            P1s[j], P2s[j], planes, inside_tol_plane)
        if hit and len_in > 0:
            f_ids.append(j+1)
            f_lengths.append(round(len_in, 4))
            fp1s.append(f'[{Aclip[0]:.4f} {Aclip[1]:.4f} {Aclip[2]:.4f}]')
            fp2s.append(f'[{Bclip[0]:.4f} {Bclip[1]:.4f} {Bclip[2]:.4f}]')
            # Fiber volume = cross-section area × length (cylinder model)
            total_vol += math.pi * r_fiber**2 * len_in

    total_fiber_volumes[i] = total_vol
    fiber_volume_fracs[i]  = total_vol / max(vol, 1e-12)  # avoid /0

    fiber_indices_col[i] = ', '.join(map(str, f_ids))
    fiber_lengths_col[i] = ', '.join(map(str, f_lengths))
    fiber_p1_col[i]      = ', '.join(fp1s)
    fiber_p2_col[i]      = ', '.join(fp2s)

df_out = pd.DataFrame({
    'RegionID':           region_ids,
    'VoronoiVolume':      voronoi_volumes,
    'TotalFiberVolume':   total_fiber_volumes,
    'FiberVolumeFrac':    fiber_volume_fracs,
    'FiberIndices':       fiber_indices_col,
    'FiberLengthsInside': fiber_lengths_col,
    'FiberP1_inside':     fiber_p1_col,
    'FiberP2_inside':     fiber_p2_col,
    'MinCorner':          min_corner_col,
    'CuboidSize':         cuboid_size_col,
})
df_out.to_csv(output_csv, index=False)
print(f'\nSaved CSV summary: {output_csv}')

# =============================================================================
# FIGURE 1: ALL FIBERS
# =============================================================================
# Shows all fiber segments rendered as 3D cylinders inside the RVE cube.
# Used to visually verify the Micro-CT data was loaded correctly and that
# fibers span a reasonable distribution within the domain.
fig1 = plt.figure(figsize=(9,7), facecolor='white')
ax1  = fig1.add_subplot(111, projection='3d')
ax1.set_title('Figure (1): All Fibers')
ax1.set_xlabel('X'); ax1.set_ylabel('Y'); ax1.set_zlabel('Z')
draw_cube_wireframe(ax1, cube_min, cube_max)

for k in range(n_fibers):
    A, B = P1s[k].copy(), P2s[k].copy()
    if clip_fibers_to_cube:
        ok, A, B = clip_segment_to_aabb(A, B, cube_min, cube_max)
        if not ok:
            continue
    draw_cylinder_segment(ax1, A, B, fiber_vis_radius, fiber_vis_n_circ, fiber_rgba)

ax1.set_box_aspect([1,1,1])
plt.tight_layout()

# =============================================================================
# GIF: SPHERE SLIDING THROUGH FIBERS
# =============================================================================
# Animates the sphere sliding through each seed position in turn.
# Each frame: the static fibers are drawn as red cylinders, the moving
# sampling sphere is shown as a translucent blue surface.
# This animation illustrates the seeding process to the reader/sponsor.
if enable_gif_generation:
    fig_gif = plt.figure(figsize=(9,7), facecolor='white')
    ax_gif  = fig_gif.add_subplot(111, projection='3d')
    ax_gif.set_title('GIF: Sphere sliding through fibers')
    ax_gif.set_xlabel('X'); ax_gif.set_ylabel('Y'); ax_gif.set_zlabel('Z')

    draw_cube_faces(ax_gif, cube_min, cube_max, [0.85,0.85,0.85], 0.08)

    pad = 20
    ax_gif.set_xlim(cube_min-pad, cube_max+pad)
    ax_gif.set_ylim(cube_min-pad, cube_max+pad)
    ax_gif.set_zlim(cube_min-pad, cube_max+pad)
    ax_gif.set_box_aspect([1,1,1])

    for k in range(n_fibers):
        A, B = P1s[k].copy(), P2s[k].copy()
        if clip_fibers_to_cube:
            ok, A, B = clip_segment_to_aabb(A, B, cube_min, cube_max)
            if not ok:
                continue
        draw_cylinder_segment(ax_gif, A, B, fiber_vis_radius, fiber_vis_n_circ, fiber_rgba)

    gif_frames = []
    n_frames = int(min(len(sphere_centers), gif_max_frames))

    for iF in range(n_frames):
        c0 = sphere_centers[iF]
        surf_handle = draw_simple_sphere(ax_gif, c0, sphere_radius, sphere_rgba)
        fig_gif.canvas.draw()

        buf = io.BytesIO()
        fig_gif.savefig(buf, format='png', dpi=72)
        buf.seek(0)
        gif_frames.append(Image.open(buf).copy())
        buf.close()
        surf_handle.remove()    # remove sphere before drawing next frame

    if gif_frames:
        gif_frames[0].save(
            gif_filename,
            save_all=True,
            append_images=gif_frames[1:],
            loop=0,
            duration=int(gif_delay_time * 1000))
        print(f'GIF saved: {gif_filename}')

    plt.close(fig_gif)

# =============================================================================
# FIGURE 2: VORONOI ONLY  (axes scaled 0→2 mm)
# =============================================================================
# Scale: internal coords [-50, 50]  →  display [0, 2] mm
# Transform:  mm = (coord - cube_min) / (cube_max - cube_min) * 2.0
_fig2_scale  = 2.0 / (cube_max - cube_min)   # 0.02
_fig2_offset = -cube_min                       # 50

# Distinct MATLAB-style pastel colors matching the reference figure
_voronoi_palette = np.array([
    [0.38, 0.75, 0.72],   # teal
    [0.20, 0.63, 0.17],   # green
    [0.58, 0.20, 0.58],   # purple
    [0.75, 0.55, 0.28],   # brown/tan
    [0.25, 0.45, 0.75],   # steel blue
    [0.85, 0.60, 0.18],   # amber
    [0.45, 0.72, 0.42],   # sage green
    [0.65, 0.30, 0.65],   # violet
    [0.30, 0.65, 0.65],   # cyan-teal
    [0.75, 0.38, 0.25],   # terracotta
    [0.35, 0.55, 0.80],   # periwinkle
    [0.60, 0.70, 0.35],   # olive green
])

fig2 = plt.figure(figsize=(9,7), facecolor='white')
ax2  = fig2.add_subplot(111, projection='3d')
ax2.set_title('Figure (2): Voronoi Regions Only')
ax2.set_xlabel('X (mm)'); ax2.set_ylabel('Y (mm)'); ax2.set_zlabel('Z (mm)')

# Draw cube wireframe in mm coords
_P_cube = np.array([
    [cube_min, cube_min, cube_min], [cube_max, cube_min, cube_min],
    [cube_max, cube_max, cube_min], [cube_min, cube_max, cube_min],
    [cube_min, cube_min, cube_max], [cube_max, cube_min, cube_max],
    [cube_max, cube_max, cube_max], [cube_min, cube_max, cube_max]])
_P_cube_mm = (_P_cube + _fig2_offset) * _fig2_scale
_edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
for e in _edges:
    ax2.plot([_P_cube_mm[e[0],0], _P_cube_mm[e[1],0]],
             [_P_cube_mm[e[0],1], _P_cube_mm[e[1],1]],
             [_P_cube_mm[e[0],2], _P_cube_mm[e[1],2]], 'k-', linewidth=1.2)

plotted     = 0
cell_colors = np.full((n_seeds, 3), np.nan)

for i in range(n_seeds):
    V, F = all_vertices[i], all_faces[i]
    if V is None or F is None:
        continue
    plotted += 1
    if plotted > max_plot_cells:
        break
    # Pick color from palette (cycles), consistent across figures
    c = _voronoi_palette[i % len(_voronoi_palette)]
    cell_colors[i] = c
    # Scale vertices to mm
    V_mm = (V + _fig2_offset) * _fig2_scale
    tris = [V_mm[tri] for tri in F]
    ec   = 'k' if show_edges else 'none'
    pc   = Poly3DCollection(tris, facecolor=list(c)+[face_alpha],
                            edgecolor=ec, linewidth=0.4 if show_edges else 0)
    ax2.add_collection3d(pc)

# Axis ticks every 0.5 mm
ax2.set_xlim(0, 2); ax2.set_ylim(0, 2); ax2.set_zlim(0, 2)
ax2.set_xticks([0, 0.5, 1.0, 1.5, 2.0])
ax2.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
ax2.set_zticks([0, 0.5, 1.0, 1.5, 2.0])
ax2.set_box_aspect([1,1,1])
plt.tight_layout()

# =============================================================================
# FIGURE 3: VORONOI + FIBER OVERLAY
# =============================================================================
# Overlays the fiber segments on the Voronoi tessellation.
# Fibers are colored by which Voronoi cell their midpoint falls into
# (same color as the cell they belong to), making the cell-fiber assignment
# immediately visible.
fig3 = plt.figure(figsize=(9,7), facecolor='white')
ax3  = fig3.add_subplot(111, projection='3d')
ax3.set_title('Figure (3): Voronoi Regions + Fiber Overlay')
ax3.set_xlabel('X (mm)'); ax3.set_ylabel('Y (mm)'); ax3.set_zlabel('Z (mm)')
# Draw cube in mm
for e in _edges:
    ax3.plot([_P_cube_mm[e[0],0], _P_cube_mm[e[1],0]],
             [_P_cube_mm[e[0],1], _P_cube_mm[e[1],1]],
             [_P_cube_mm[e[0],2], _P_cube_mm[e[1],2]], 'k-', linewidth=1.2)
seeds_mm = (seeds + _fig2_offset) * _fig2_scale
ax3.scatter(seeds_mm[:,0], seeds_mm[:,1], seeds_mm[:,2], s=18, c='k')

plotted = 0
for i in range(n_seeds):
    V, F = all_vertices[i], all_faces[i]
    if V is None or F is None:
        continue
    plotted += 1
    if plotted > max_plot_cells:
        break
    c  = cell_colors[i]
    if np.any(np.isnan(c)):
        c = _voronoi_palette[i % len(_voronoi_palette)]
    ec   = 'k' if show_edges else 'none'
    V_mm = (V + _fig2_offset) * _fig2_scale
    tris = [V_mm[tri] for tri in F]
    pc   = Poly3DCollection(tris, facecolor=list(c)+[face_alpha],
                            edgecolor=ec, linewidth=0.4 if show_edges else 0)
    ax3.add_collection3d(pc)

# Assign each fiber to its nearest seed using a KD-tree for efficiency
mid_pts  = 0.5 * (P1s + P2s)
kdt      = KDTree(seeds)
idx_cell = kdt.query(mid_pts, k=1)[0].ravel().astype(int)

for k in range(n_fibers):
    A, B = P1s[k].copy(), P2s[k].copy()
    if clip_fibers_to_cube:
        ok, A, B = clip_segment_to_aabb(A, B, cube_min, cube_max)
        if not ok:
            continue
    if fiber_color_mode == "byCell":
        ci = int(idx_cell[k])
        c  = cell_colors[ci]
        if np.any(np.isnan(c)):
            c = np.zeros(3)
    else:
        c = fiber_uniform_color
    A_mm = (A + _fig2_offset) * _fig2_scale
    B_mm = (B + _fig2_offset) * _fig2_scale
    ax3.plot([A_mm[0],B_mm[0]], [A_mm[1],B_mm[1]], [A_mm[2],B_mm[2]],
             '-', linewidth=fiber_line_width, color=c)

ax3.set_xlim(0, 2); ax3.set_ylim(0, 2); ax3.set_zlim(0, 2)
ax3.set_xticks([0, 0.5, 1.0, 1.5, 2.0])
ax3.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
ax3.set_zticks([0, 0.5, 1.0, 1.5, 2.0])
ax3.set_box_aspect([1,1,1])
plt.tight_layout()

# =============================================================================
# FIGURE 4: SPECIFIC REGION + SUMMARY TABLE
# =============================================================================
# Interactive: the user enters a region index (1-based) and the script plots
# that single Voronoi cell with all fiber segments clipped to its interior.
# A summary table on the right shows volume, fiber count, volume fraction,
# and the precise clipped fiber endpoint coordinates.
# This is used to validate individual cells before running Abaqus.
try:
    region_idx = int(input('\nEnter Voronoi region index (e.g., 1 for region_001): '))
except (ValueError, EOFError):
    region_idx = -1

if (region_idx < 1 or region_idx > n_seeds or all_vertices[region_idx-1] is None):
    warnings.warn('Invalid region index or region is empty. Skipping Figure (4).')
else:
    ri = region_idx - 1
    V  = all_vertices[ri]
    F  = all_faces[ri]

    vor_vol     = voronoi_volumes[ri]
    tot_fib_vol = total_fiber_volumes[ri]
    fvf         = fiber_volume_fracs[ri]

    centroid = V.mean(axis=0)
    planes   = build_convex_halfspaces_from_faces(V, F, centroid)

    f_ids      = []
    f_lengths  = []
    fp1_inside = []
    fp2_inside = []

    for j in range(n_fibers):
        hit, Aclip, Bclip, len_in = clip_segment_with_convex_halfspaces(
            P1s[j], P2s[j], planes, inside_tol_plane)
        if hit and len_in > 0:
            f_ids.append(j+1)
            f_lengths.append(len_in)
            fp1_inside.append(Aclip.copy())
            fp2_inside.append(Bclip.copy())

    fp1_inside = np.array(fp1_inside) if fp1_inside else np.zeros((0,3))
    fp2_inside = np.array(fp2_inside) if fp2_inside else np.zeros((0,3))

    fig4 = plt.figure(figsize=(16, 8), facecolor='white')
    fig4.suptitle(f'Figure (4): region{region_idx:03d} (Cell + Fibers Inside)', fontsize=13)

    ax4 = fig4.add_axes([0.04, 0.08, 0.55, 0.85], projection='3d')
    ax4.set_xlabel('X'); ax4.set_ylabel('Y'); ax4.set_zlabel('Z')

    # No grid, no background panes, no outer box
    ax4.grid(False)
    ax4.set_axis_off()

    tris = [V[tri] for tri in F]
    pc   = Poly3DCollection(tris, facecolor = [0.04, 0.27, 0.39, 0.55],
                            edgecolor='none', linewidth=0)   # no black edge lines
    ax4.add_collection3d(pc)
    ax4.scatter([seeds[ri,0]], [seeds[ri,1]], [seeds[ri,2]], s=70, c='k')

    for ii in range(len(f_ids)):
        Aclip = fp1_inside[ii]
        Bclip = fp2_inside[ii]
        if clip_fibers_to_cube:
            ok, Aplot, Bplot = clip_segment_to_aabb(Aclip, Bclip, cube_min, cube_max)
            if not ok:
                Aplot, Bplot = Aclip, Bclip
        else:
            Aplot, Bplot = Aclip, Bclip
        draw_cylinder_segment(
                                ax4,
                                Aplot,
                                Bplot,
                                fiber_vis_radius ,            # <-- this is already 4.0 in your code
                                fiber_vis_n_circ,   # smoothness (18 is good)
                                [0.85, 0.1, 0.1, 1.0]  # red solid fiber
                            )

    ax4.set_box_aspect([1,1,1])

    # Right side: summary tables (Region properties, fiber indices, lengths, endpoints)
    summary_data = [
        ['RegionID',        str(region_idx)],
        ['VoronoiVolume',   f'{vor_vol:.6g}'],
        ['TotalFiberVol',   f'{tot_fib_vol:.6g}'],
        ['FiberVolumeFrac', f'{fvf:.6g}'],
        ['NumFibersInside', str(len(f_ids))],
        ['FiberRadius',     str(r_fiber)],
        ['CSV_Output',      output_csv],
    ]
    ax_tbl1 = fig4.add_axes([0.62, 0.76, 0.36, 0.22])
    ax_tbl1.axis('off')
    tbl1 = ax_tbl1.table(cellText=summary_data, colLabels=['Property', 'Value'],
                         loc='center', cellLoc='left')
    tbl1.auto_set_font_size(False); tbl1.set_fontsize(8); tbl1.scale(1, 1.4)

    ax_tbl2 = fig4.add_axes([0.62, 0.52, 0.17, 0.22])
    ax_tbl2.axis('off')
    fi_data = [[str(x)] for x in f_ids]
    tbl2 = ax_tbl2.table(cellText=fi_data, colLabels=['FiberIdx'],
                         loc='center', cellLoc='center')
    tbl2.auto_set_font_size(False); tbl2.set_fontsize(8); tbl2.scale(1, 1.2)

    ax_tbl3 = fig4.add_axes([0.80, 0.52, 0.17, 0.22])
    ax_tbl3.axis('off')
    fl_data = [[f'{x:.4f}'] for x in f_lengths]
    tbl3 = ax_tbl3.table(cellText=fl_data, colLabels=['LenInside'],
                         loc='center', cellLoc='center')
    tbl3.auto_set_font_size(False); tbl3.set_fontsize(8); tbl3.scale(1, 1.2)

    ax_tbl4 = fig4.add_axes([0.62, 0.27, 0.36, 0.22])
    ax_tbl4.axis('off')
    p1_data = [[f'[{r[0]:.4f} {r[1]:.4f} {r[2]:.4f}]'] for r in fp1_inside] \
              if len(fp1_inside) else [['—']]
    tbl4 = ax_tbl4.table(cellText=p1_data, colLabels=['FiberP1_inside (Start)'],
                         loc='center', cellLoc='left')
    tbl4.auto_set_font_size(False); tbl4.set_fontsize(7); tbl4.scale(1, 1.2)

    ax_tbl5 = fig4.add_axes([0.62, 0.03, 0.36, 0.22])
    ax_tbl5.axis('off')
    p2_data = [[f'[{r[0]:.4f} {r[1]:.4f} {r[2]:.4f}]'] for r in fp2_inside] \
              if len(fp2_inside) else [['—']]
    tbl5 = ax_tbl5.table(cellText=p2_data, colLabels=['FiberP2_inside (End)'],
                         loc='center', cellLoc='left')
    tbl5.auto_set_font_size(False); tbl5.set_fontsize(7); tbl5.scale(1, 1.2)

    print(f'Region {region_idx:03d}:\n'
          f'  VoronoiVolume    = {vor_vol:.6g}\n'
          f'  TotalFiberVolume = {tot_fib_vol:.6g}\n'
          f'  FiberVolumeFrac  = {fvf:.6g}\n'
          f'  FibersInside     = {len(f_ids)}')

plt.show()
