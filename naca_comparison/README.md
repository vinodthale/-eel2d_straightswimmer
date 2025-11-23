# NACA Airfoil Method Comparison

This directory contains a comprehensive comparison between the **original/standard NACA 4-digit airfoil method** and the **IBAMR perpendicular stacking method** currently implemented in the eel2d project.

## Overview

### Method 1: Original NACA 4-Digit Method (Standard Textbook)

The standard NACA 4-digit airfoil formulation as described in aerospace engineering textbooks:

**Key Equations:**

1. **Thickness Distribution** (for NACA 00xx symmetric profiles):
   ```
   y_t = 5t [a0*√x + a1*x + a2*x² + a3*x³ + a4*x⁴]
   ```
   where:
   - t = maximum thickness as fraction of chord
   - Coefficients: a0=0.2969, a1=-0.1260, a2=-0.3516, a3=0.2843, a4=-0.1015

2. **Camber Line** (for cambered airfoils):
   ```
   y_c = f(x)  [manuscript uses: y_c = Amax * ((x+0.03125)/1.03125) * sin(2πx)]
   ```

3. **Surface Construction** (perpendicular to camber line):
   ```
   θ = atan(dy_c/dx)
   x_upper = x - y_t * sin(θ)
   y_upper = y_c + y_t * cos(θ)
   x_lower = x + y_t * sin(θ)
   y_lower = y_c - y_t * cos(θ)
   ```

**Characteristics:**
- ✓ Smooth, continuous surface representation
- ✓ Direct analytical formulation
- ✓ Minimal point count for surface definition
- ✓ Standard aerospace industry method
- ✗ No internal structure (just surface)
- ✗ Not directly compatible with IBAMR volume discretization

---

### Method 2: IBAMR Perpendicular Stacking Method (Current Implementation)

The discretized layering approach for IBAMR immersed boundary simulations:

**Key Approach:**

1. **Determine Layer Count** at each chordwise station:
   ```
   k = floor(y_t / ds)
   ```
   where ds = dx_fine (finest mesh spacing)

2. **Stack Layers Perpendicular to Camber**:
   - Upper layers: j = 0, 1, 2, ..., k (total: k+1 points)
   - Lower layers: j = 1, 2, ..., k (total: k points)

3. **Vertex Positions**:
   ```
   For each layer j:
     offset = j * ds
     Upper: [x, y_c] + offset * [n_x, n_y]
     Lower: [x, y_c] - offset * [n_x, n_y]
   ```
   where n = [-sin(θ), cos(θ)] is the unit normal

**Characteristics:**
- ✓ Full 3D volumetric representation
- ✓ Internal structure for fluid-structure coupling
- ✓ Mesh-adaptive (ds tied to finest grid resolution)
- ✓ IBAMR-compatible vertex ordering
- ✗ Higher point count (~2-10× more points)
- ✗ Discretization-dependent (accuracy depends on ds)
- ✗ Outer envelope doesn't exactly match analytical surface

---

## Key Differences

| Aspect | Original Method | Stacking Method |
|--------|----------------|-----------------|
| **Purpose** | Analytical airfoil definition | CFD immersed boundary simulation |
| **Output** | Surface coordinates only | Full volume with internal layers |
| **Point Distribution** | Continuous surface (256-512 pts) | Discrete layers (1000-6000 pts) |
| **Accuracy** | Exact analytical formula | Discretization error ~O(ds) |
| **Mesh Dependency** | Independent | Coupled to CFD mesh (ds = dx_fine) |
| **IBAMR Compatibility** | No | Yes |
| **Typical Use** | Design, analysis, manufacturing | Fluid-structure interaction CFD |

---

## Volume/Area Comparison

Both methods compute airfoil cross-sectional area, but differently:

### Analytical Area (Exact):
```matlab
A_exact = ∫₀¹ 2*y_t dx = ∫₀¹ 10t [a0*√x + a1*x + a2*x² + a3*x³ + a4*x⁴] dx
```

### Original Method Area:
Uses polygon area of smooth upper/lower surfaces:
```matlab
A_orig = polyarea([x_u, flip(x_l)], [y_u, flip(y_l)])
```

### Stacking Method Areas:
1. **Surface Area**: Polygon of analytical surface points on stacking grid
2. **Envelope Area**: Polygon of outermost layer (j=k_max at each x)

**Expected Results:**
- Original method: Error < 1e-6 (very accurate)
- Stacking surface: Error ~ 1e-6 to 1e-5 (good)
- Stacking envelope: Error ~ 1e-4 to 1e-3 (discretization dependent)

---

## Files in This Directory

### 1. `naca_method_comparison.m`
**Main comparison script** - Runs both methods side-by-side and generates:
- Comprehensive comparison plots (6 subplots per thickness)
- Detailed quantitative comparison tables
- Error analysis and computational efficiency metrics

**Usage:**
```matlab
cd naca_comparison
naca_method_comparison
```

**Outputs:**
- `comparison_plots/comparison_NACA00XX.png` - Detailed visualizations
- Console tables showing volume, error, and point count comparisons

### 2. `README.md` (this file)
Complete documentation of both methods and their comparison.

---

## Understanding the Results

### Subplot Guide (Comparison Plots)

Each generated plot contains 6 subplots:

1. **Top-Left**: Original method full airfoil (chord frame)
2. **Top-Center**: Stacking method full airfoil (world frame)
3. **Top-Right**: Thickness distribution comparison
4. **Bottom-Left**: Leading edge detail (original, showing point distribution)
5. **Bottom-Center**: Leading edge detail (stacking, showing layering)
6. **Bottom-Right**: Area/volume bar chart comparison

### Key Metrics to Check

1. **Volume Accuracy**:
   - How close is Vf_surface to Vf_analytic?
   - Both methods should be within 0.1% for validation

2. **Point Count Ratio**:
   - Stacking typically uses 2-10× more points
   - Ratio increases with thickness (thicker profiles need more layers)

3. **Layer Count (k_max)**:
   - Maximum number of layers at thickest point
   - k_max = floor(y_t_max / ds)
   - For t=0.24, expect k_max ~ 15-25 layers

---

## References

### Standard NACA Documentation:
1. Abbott, I. H., & Von Doenhoff, A. E. (1959). *Theory of Wing Sections*. Dover Publications.
2. NACA Report 824: "Summary of Airfoil Data" (1945)
3. ERAU Aerospace Engineering Textbook: [Airfoil Geometries Chapter](https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/airfoil-geometries/)

### IBAMR Methodology:
4. Griffith, B. E., & Luo, X. (2017). "Hybrid finite difference/finite element immersed boundary method." *International Journal for Numerical Methods in Biomedical Engineering*, 33(12), e2888.
5. IBAMR Documentation: https://ibamr.github.io/

---

## Validation Checklist

When running the comparison, verify:

- [ ] Vf_analytic matches theoretical value for each thickness
- [ ] Vf_orig_surface within 1e-6 of Vf_analytic
- [ ] Vf_stack_surface within 1e-5 of Vf_analytic
- [ ] Leading edge points properly clustered (original method)
- [ ] Perpendicular layering visible in stacking plots
- [ ] k_max increases proportionally with thickness
- [ ] No self-intersecting polygons in either method
- [ ] Smooth airfoil shape in all visualizations

---

## When to Use Each Method

### Use Original Method When:
- Designing airfoils for aerodynamic performance
- Creating CAD models for manufacturing
- Performing 2D potential flow analysis
- Need exact analytical surface definition
- Minimizing point count is critical

### Use Stacking Method When:
- Running IBAMR fluid-structure interaction simulations
- Need internal structure for Lagrangian forcing
- Coupling to Eulerian fluid mesh
- Simulating flexible/deformable bodies
- Mesh-adaptive refinement required

---

## Contact & Issues

For questions about this comparison:
- Review the manuscript methodology section
- Check IBAMR documentation for vertex file format
- Verify grid parameters (MAX_LEVELS, REF_RATIO, ds)

---

**Last Updated:** 2025-11-23
**Compatibility:** MATLAB R2019b or newer
