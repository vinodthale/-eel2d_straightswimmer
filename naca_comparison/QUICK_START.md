# Quick Start Guide - NACA Comparison

## What's in This Directory?

This directory contains a comprehensive comparison between two NACA airfoil generation methods:

1. **Original/Standard Method** - Classic aerospace engineering approach
2. **IBAMR Stacking Method** - Current implementation for fluid-structure simulations

## How to Run

### Option 1: Full Comparison (Recommended First Run)

```matlab
cd naca_comparison
naca_method_comparison
```

**What it does:**
- Compares both methods for 6 thickness ratios (0.04, 0.08, 0.12, 0.16, 0.20, 0.24)
- Generates detailed comparison plots (6 subplots each)
- Outputs comprehensive tables (volume, error, computational cost)
- Runtime: ~30-60 seconds

**Output:**
- `comparison_plots/comparison_NACA00XX.png` (6 files)
- Console tables with quantitative metrics

---

### Option 2: Original Method Only

```matlab
cd naca_comparison
naca_original_method
```

**What it does:**
- Implements standard NACA 4-digit formulation
- Generates airfoils with cosine point clustering
- Saves coordinate files in standard format
- Runtime: ~10-20 seconds

**Output:**
- `original_method_plots/NACA00XX_original.png` (plots)
- `original_method_data/NACA00XX_coords.dat` (coordinate files)
- `original_method_data/NACA00XX_area.txt` (area analysis)

---

### Option 3: Validation Test (Single Airfoil)

```matlab
cd naca_comparison
validation_test
```

**What it does:**
- Detailed analysis of NACA 0012 (single airfoil)
- Comprehensive diagnostics and metrics
- Side-by-side visualization
- Runtime: ~5-10 seconds

**Output:**
- `validation_test_NACA0012.png`
- Detailed console diagnostics

---

## Understanding the Results

### Key Metrics to Check

1. **Volume/Area Accuracy**
   - Look for: `Error < 1e-5` for both methods
   - Original method typically: ~1e-7 to 1e-6
   - Stacking method typically: ~1e-6 to 1e-5

2. **Point Count**
   - Original: ~512 points (256 upper + 256 lower)
   - Stacking: ~1000-6000 points (depends on thickness)
   - Ratio: Stacking uses 2-10× more points

3. **Layer Count (k_max)**
   - Thicker airfoils → more layers
   - NACA 0012: ~10-15 layers
   - NACA 0024: ~20-30 layers

### What the Plots Show

**Full Comparison** (`naca_method_comparison.m`):
1. Top-left: Original method (smooth surface, cosine clustering)
2. Top-center: Stacking method (discrete layers, IBAMR format)
3. Top-right: Thickness distribution
4. Bottom-left: Leading edge detail (original)
5. Bottom-center: Leading edge detail (stacking with layers visible)
6. Bottom-right: Area comparison bar chart

**Validation Test** (`validation_test.m`):
- 6 diagnostic subplots
- Error analysis
- Layer distribution
- Point count comparison

---

## Customization

### Change Thickness Values

Edit the `thickness_list` variable:
```matlab
thickness_list = [0.08, 0.12, 0.16];  % Just 3 airfoils
```

### Change Camber Line

In `naca_original_method.m`, set:
```matlab
camber_type = 'none';        % Symmetric (flat camber)
camber_type = 'manuscript';  % Current manuscript midline
camber_type = 'naca2412';    % Standard NACA 4-digit camber
```

### Adjust Resolution

**Original method:**
```matlab
Nx = 512;  % More points (default: 256)
```

**Stacking method:**
```matlab
MAX_LEVELS = 5;  % Finer grid → smaller ds → more layers
```

---

## Interpreting Table Output

### Part 1: Volume/Area Analysis
```
  t     Vf_analytic      Vf_orig         Vf_stack_surf   Vf_stack_env
-----------------------------------------------------------------------------------
 0.12   1.297812345e-01  1.297812340e-01 1.297812355e-01 1.297815000e-01
```
- `Vf_analytic`: Exact integral of NACA formula
- `Vf_orig`: Polygon area (original method surface)
- `Vf_stack_surf`: Polygon area (stacking analytical surface)
- `Vf_stack_env`: Polygon area (outermost layer envelope)

**What to check:** All values should agree within ~0.01%

### Part 2: Error Analysis
```
  t     Error_orig       Error_stack     Error_env       Error_ratio
-------------------------------------------------------------------------
 0.12   1.234567e-07     4.567890e-06    1.234567e-04    37.0
```
- `Error_orig`: |Vf_orig - Vf_analytic|
- `Error_stack`: |Vf_stack_surf - Vf_analytic|
- `Error_env`: |Vf_stack_env - Vf_analytic| (discretization error)
- `Error_ratio`: Error_stack / Error_orig

**What to check:** Error_ratio should be 1-100 (stacking slightly less accurate)

### Part 3: Computational Efficiency
```
  t     Npts_orig   Npts_stack   Point_ratio   Max_layers
---------------------------------------------------------------
 0.12        512         2048          4.0            15
```
- Higher thickness → more layers → more points
- Point_ratio = Npts_stack / Npts_orig
- Max_layers = k_max (at thickest point)

**What to check:** Point_ratio increases with thickness

---

## Troubleshooting

### "Unable to open file for writing"
- Check that you have write permissions
- Make sure directories exist (run `mkdir comparison_plots` if needed)

### "Out of memory"
- Reduce `Nx_chord` in stacking method
- Process fewer thickness values
- Close other MATLAB figures

### Plots look strange
- Check that `axis equal` is set
- Verify domain placement (`Xcenter`, `Ycenter`)
- Ensure no NaN values in coordinates

### Different results than expected
- Verify NACA coefficients (a0, a1, a2, a3, a4)
- Check camber line formulation
- Confirm ds calculation (dx_fine)

---

## Next Steps

1. **Run full comparison** to see all 6 airfoils
2. **Review validation test** for detailed diagnostics
3. **Read README.md** for theoretical background
4. **Customize** for your specific needs

---

## Quick Reference

| File | Purpose | Runtime | Output |
|------|---------|---------|--------|
| `naca_method_comparison.m` | Full comparison (6 airfoils) | 30-60s | Plots + tables |
| `naca_original_method.m` | Original method only | 10-20s | Coords + plots |
| `validation_test.m` | Single airfoil diagnostics | 5-10s | Diagnostics + plot |
| `README.md` | Theoretical documentation | - | - |

---

**Questions?** Check the main README.md for detailed explanations of the methods and equations.
