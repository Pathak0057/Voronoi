# Micro-CT → Fiber Network → Voronoi Regions → Abaqus RVE Pipeline

This repository contains scripts and example data to process **micro-CT–derived fiber centrelines**, generate **bounded 3D Voronoi regions / subdomains**, compute **regional statistics** (e.g., local fiber volume fraction), and export geometry + inputs for **Abaqus-based mechanical simulations** of short-fiber composites.

> **Goal:** Provide a transparent, reusable pipeline that converts measured fiber microstructure into simulation-ready RVEs and region-wise descriptors for modeling and surrogate learning.

---

## Contents (What’s in this repo)

### Core scripts
- **`Micro_ct_to_Voronoi.m`**  
  Main MATLAB pipeline: reads fibre centrelines, prepares fiber geometry, creates seeding strategy, builds bounded Voronoi regions, and exports region summaries.

- **`clipSegmentWithPolyhedron.m`**  
  Utility function to clip fiber segments by a polyhedral domain (used to keep fibers inside the cuboid/RVE or within Voronoi cells).

- **`loadSTLsmart.m`**  
  Robust STL import helper for MATLAB (handles common STL issues).

- **`Voronoi_Abaqus.py`**  
  Python script for Abaqus CAE: imports generated geometry/region definitions and automates meshing + job setup (supports region-wise processing).

- **`STL_refinement_Freecad.py`**  
  Optional FreeCAD-based STL cleanup/refinement step prior to Abaqus import.

### Example data (replace with your own)
- **`final_fiber_data_CentrelinePoints.xlsx`**  
  Example fiber centreline points exported from micro-CT processing.

- **`final_fiber_data_CentrelinePoints_Abaqus_static_S1_AF.xlsx`**  
  Example dataset formatted for Abaqus static analysis workflow.

- **`materials_db.csv`**  
  Material database / parameter table used by scripts (edit for your material system).

---

## Workflow Overview

1. **Micro-CT centreline data** (fiber startpoints/ endpoints)  
2. **Geometry preparation**: fiber segmentation + bounding cuboid enforcement  
3. **Seeding strategy**: seed points generated from fiber midpoints/centroids  
4. **3D Voronoi tessellation** within bounded cuboid  
5. **Region statistics**: local volume fractions, region-wise fiber lists  
6. **Export** to Abaqus / CAD:
   - region-wise geometry (STL/STEP as configured)
   - region summary CSV files
7. **Abaqus automation**: mesh, assign sets/surfaces, apply BCs, run jobs

---

## Requirements

### MATLAB
- MATLAB R2020a+ recommended

### Python
- Python 3.9+ recommended
- Abaqus Python (comes with Abaqus/CAE) for `Voronoi_Abaqus.py`
- FreeCAD (optional) for `STL_refinement_Freecad.py`

### Software
- Abaqus/CAE 2020+ recommended  
- FreeCAD (optional) for STL refinement

---

## Quick Start (Typical Use)

### 1) Prepare input
Place your fibre centreline file in the repo root (or update paths inside scripts):
- `final_fiber_data_CentrelinePoints.xlsx`

Expected minimum fields typically include:
- Fiber ID
- (x,y,z) for start / end, or ordered centreline points
- Units must be consistent (e.g., µm or mm)

### 2) Run MATLAB pipeline
In MATLAB:
```matlab
run('Micro_ct_to_Voronoi.m');
