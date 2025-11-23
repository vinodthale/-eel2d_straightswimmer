# 2D Eel Straight Swimmer

This MATLAB script generates geometry for a 2D eel simulation using the Immersed Boundary method. It is taken from the [IBAMR](https://github.com/IBAMR/IBAMR) (Immersed Boundary Adaptive Mesh Refinement) framework.

## Source

This code is from the IBAMR project:
- **Original Source**: [IBAMR/examples/ConstraintIB/eel2d/eel2d_straightswimmer.m](https://github.com/IBAMR/IBAMR/blob/master/examples/ConstraintIB/eel2d/eel2d_straightswimmer.m)
- **Created By**: Amneet Bhalla
- **License**: 3-clause BSD license (see IBAMR COPYRIGHT file)

## Description

This script generates the geometry and material points for a 2D eel body that swims in a straight line. The eel geometry is designed for use in fluid-structure interaction (FSI) simulations using the constraint immersed boundary (CIB) method.

### Key Features

- **2D Eel Geometry Generation**: Creates a streamlined eel body with proper head and tail tapering
- **Material Point Distribution**: Distributes Lagrangian material points along the eel body based on the background Eulerian mesh spacing
- **Swimming Kinematics**: Incorporates sinusoidal body deformation for swimming motion
- **Output Format**: Generates vertex file compatible with IBAMR input requirements

## Technical Details

### Geometry Parameters

- **Body Length (L)**: 1.0 unit
- **Head Width (Wh)**: 0.04 × L (4% of body length)
- **Head Length (Xh)**: 0.04 units
- **Background Mesh**: 8 × 4 domain with resolution controlled by refinement factors

### Body Sections

1. **Head Section** (0 to Xh):
   - Elliptical cross-section with width = √(2·Wh·x - x²)
   - Provides smooth nose geometry

2. **Tail Section** (Xh to L):
   - Linear taper from head width to tail
   - Width decreases as: Wh·(L - x)/(L - Xh)

### Swimming Motion

The centerline follows a sinusoidal pattern:
```
y = 0.125 × ((x + 0.03125)/1.03125) × sin(2πx)
```

This creates a traveling wave along the body for forward propulsion.

## Usage

### Prerequisites

- MATLAB (any recent version)

### Running the Script

1. Run the script in MATLAB:
   ```matlab
   eel2d_straightswimmer
   ```

2. The script will:
   - Generate the eel geometry
   - Display a plot of the eel body with material points
   - Create an output file: `eel2d.vertex`

### Output

The script generates `eel2d.vertex` containing:
- First line: Total number of material points
- Subsequent lines: x and y coordinates of each material point (tab-separated)

## Integration with IBAMR

This vertex file is designed to be used with IBAMR simulations:

1. Use the generated `eel2d.vertex` file as input to IBAMR
2. Configure the IBAMR input file to reference this vertex file
3. Set up appropriate fluid properties and simulation parameters
4. Run the IBAMR constraint IB simulation

## Background Mesh Configuration

The mesh resolution can be adjusted by modifying the refinement factors:
```matlab
Ny = 16*4*4*4; Ly = 4;  % 4096 cells in y-direction
Nx = 32*4*4*4; Lx = 8;  % 8192 cells in x-direction
```

Higher refinement provides better resolution but increases computational cost.

## License

Copyright (c) 2014 - 2019 by the IBAMR developers
All rights reserved.

This file is part of IBAMR.

IBAMR is free software and is distributed under the 3-clause BSD license. The full text of the license can be found in the file COPYRIGHT at the top level directory of IBAMR.

## References

For more information about IBAMR and the Immersed Boundary method:
- [IBAMR GitHub Repository](https://github.com/IBAMR/IBAMR)
- [IBAMR Documentation](https://ibamr.github.io/)

## Volume and Area Verification

The repository includes comprehensive tests for verifying geometric properties of the generated eel body.

### Running Volume and Area Tests

To run all volume and area verification tests:

```matlab
RUN_VOLUME_AREA_TESTS
```

This will execute:
1. Main eel geometry generation with area/volume calculations
2. Detailed verification tests with multiple calculation methods

### Test Coverage

The verification tests include:

#### 1. **Area Calculations**
   - **Polygon method**: Using MATLAB's `polyarea()` on body envelope
   - **Trapezoidal integration**: Integrating cross-sectional areas along body length
   - **Analytical approximation**: Comparing with theoretical formulas
   - **Material point method**: Area based on point count and grid spacing

#### 2. **Volume Calculations**
   - 2D area × unit depth (for 3D volume estimation)
   - Multiple methods for cross-validation

#### 3. **Section-Wise Analysis**
   - Head section area (elliptical cross-section)
   - Tail section area (linear taper)
   - Conservation verification (head + tail = total)

#### 4. **Geometric Properties**
   - Perimeter calculation
   - Compactness metric (isoperimetric ratio)
   - Aspect ratio
   - Cross-sectional height distribution

#### 5. **Discretization Verification**
   - Height discretization error analysis
   - Comparison of theoretical vs actual cross-sections
   - Grid resolution effects

### Output Files

Running the tests generates:

```
output/
├── vertices/
│   └── eel2d.vertex                      # Lagrangian material points
├── images/
│   ├── eel2d_straightswimmer.png         # Main geometry plot
│   └── volume_area_verification.png      # Verification plots
└── data/
    ├── eel_geometry_properties.txt       # Basic geometric properties
    └── volume_area_data.txt              # Detailed verification data
```

### Validation Metrics

The tests verify:
- Area calculation accuracy (< 0.1% error between methods)
- Volume conservation
- Height discretization error
- Section-wise area conservation

### Individual Test Files

- **`eel2d_straightswimmer.m`**: Main script with integrated area/volume calculations
- **`test_volume_area.m`**: Comprehensive verification test with multiple methods and visualizations

## Related Examples

This is part of the ConstraintIB examples in IBAMR. See other examples in the IBAMR repository for similar fluid-structure interaction simulations.
