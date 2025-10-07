# Multilayer VCT Implementation

This document describes the new multilayer capabilities added to the VCT (Voxel Cellular Tissue) simulation framework.

## Overview

The multilayer VCT implementation extends the original single-layer cellular Potts model to support multiple layers with cross-layer interactions. Each voxel can now interact with neighboring voxels in adjacent layers, enabling more realistic 3D tissue modeling.

## Key Features

### 1. **Cross-Layer Interactions**
- Voxels interact with neighbors in layers above and below (z±1)
- Separate coupling parameters for different cell type interactions:
  - `JZ_CMCM`: Cardiomyocyte-Cardiomyocyte cross-layer coupling
  - `JZ_FBFB`: Fibroblast-Fibroblast cross-layer coupling  
  - `JZ_CMFB`: Cardiomyocyte-Fibroblast cross-layer coupling

### 2. **Multiple Layer Types**
- **Basic layers**: Standard X×Y cell configuration with uniform probability
- **Probability matrix layers**: Custom probability distribution per voxel
- **Ready layers**: Predefined cell configurations (future enhancement)

### 3. **Backward Compatibility**
- Original `py_cpmfem()` function continues to work unchanged
- New `py_cpmfem_layers()` function for multilayer simulations

## API Reference

### New Function: `py_cpmfem_layers()`

```python
types_list, ctags, fibers, contacts = py_cpmfem_layers(
    layers,                    # List of layer configurations
    VOXSIZE,                   # Voxel size in mm
    sizeX, sizeY,              # Simulation domain size in mm
    scenario,                  # Simulation scenario
    NRINC,                     # Number of simulation steps
    **kwargs                   # Additional parameters including cross-layer coupling
)
```

#### Parameters

- **`layers`**: List of dictionaries, each defining one layer:
  ```python
  layers = [
      {
          "type": "basic",           # Layer type
          "NCX": 10,                 # Number of cells in X
          "NCY": 10,                 # Number of cells in Y  
          "PART": 0.3               # Fibroblast percentage
      },
      {
          "type": "prob",
          "NCX": 10,
          "NCY": 10,
          "PART_matrix": matrix      # Custom probability matrix
      }
  ]
  ```

- **Cross-layer coupling parameters** (optional):
  - `JZ_CMCM`: CM-CM cross-layer coupling strength (default: 0.0)
  - `JZ_FBFB`: FB-FB cross-layer coupling strength (default: 0.0)
  - `JZ_CMFB`: CM-FB cross-layer coupling strength (default: 0.0)

#### Returns

- **`types_list`**: List of arrays, one per layer, containing cell types
- **`ctags`**: 3D array with shape `(NL, NVY, NVX)` - cell tags per voxel
- **`fibers`**: 3D array with shape `(NL, NVY, NVX)` - fiber presence per voxel
- **`contacts`**: 3D array with shape `(NL, NVY, NVX)` - contact information per voxel

## Layer Types

### 1. Basic Layers
Standard configuration with uniform cell distribution:

```python
{
    "type": "basic",
    "NCX": 10,        # Number of cells in X direction
    "NCY": 10,        # Number of cells in Y direction
    "PART": 0.3       # Percentage of fibroblasts (0.0-1.0)
}
```

### 2. Probability Matrix Layers
Custom probability distribution for cell placement:

```python
{
    "type": "prob", 
    "NCX": 10,
    "NCY": 10,
    "PART_matrix": probability_matrix  # 2D numpy array (NVY x NVX)
}
```

The `PART_matrix` should have values between 0.0 and 1.0, where:
- 0.0 = 100% cardiomyocytes at that voxel
- 1.0 = 100% fibroblasts at that voxel
- 0.5 = 50% chance of fibroblast at that voxel

### 3. Ready Layers (Future)
Predefined cell configurations:

```python
{
    "type": "ready",
    "ready_ctag": ctag_array,    # Predefined cell tags (NVY x NVX)
    "ready_type": type_array     # Predefined cell types (NCX*NCY)
}
```

## Usage Examples

### Example 1: Basic Multilayer Simulation

```python
from libcpmfem import py_cpmfem_layers
import numpy as np

# Create 3 layers with different fibroblast percentages
layers = [
    {"type": "basic", "NCX": 8, "NCY": 8, "PART": 0.2},  # 20% FB
    {"type": "basic", "NCX": 8, "NCY": 8, "PART": 0.5},  # 50% FB  
    {"type": "basic", "NCX": 8, "NCY": 8, "PART": 0.8}   # 80% FB
]

# Run simulation with cross-layer coupling
types_list, ctags, fibers, contacts = py_cpmfem_layers(
    layers=layers,
    VOXSIZE=0.0025,
    sizeX=1.0,
    sizeY=1.0,
    scenario="monolayer_on_fiber",
    NRINC=1000,
    JZ_CMCM=0.3,    # Strong CM-CM coupling
    JZ_FBFB=0.2,    # Moderate FB-FB coupling
    JZ_CMFB=0.4     # Strong CM-FB coupling
)

print(f"Simulation completed with {len(types_list)} layers")
print(f"Output shapes: ctags={ctags.shape}, fibers={fibers.shape}")
```

### Example 2: Custom Probability Distributions

```python
# Create gradient probability matrix
nvx, nvy = 50, 50
gradient_matrix = np.zeros((nvy, nvx))
for i in range(nvy):
    for j in range(nvx):
        gradient_matrix[i, j] = j / (nvx - 1)  # Left-to-right gradient

layers = [
    {
        "type": "prob",
        "NCX": 6,
        "NCY": 6, 
        "PART_matrix": gradient_matrix
    }
]

types_list, ctags, fibers, contacts = py_cpmfem_layers(
    layers=layers,
    VOXSIZE=0.0025,
    sizeX=1.0,
    sizeY=1.0,
    scenario="monolayer_on_fiber",
    NRINC=500
)
```

### Example 3: Mixed Layer Types

```python
# Create center-focused probability matrix
center_matrix = np.zeros((50, 50))
center_x, center_y = 25, 25
for i in range(50):
    for j in range(50):
        dist = np.sqrt((i - center_y)**2 + (j - center_x)**2)
        center_matrix[i, j] = max(0, 1.0 - dist/20.0)

layers = [
    {"type": "basic", "NCX": 8, "NCY": 8, "PART": 0.3},
    {"type": "prob", "NCX": 8, "NCY": 8, "PART_matrix": center_matrix},
    {"type": "basic", "NCX": 8, "NCY": 8, "PART": 0.7}
]

types_list, ctags, fibers, contacts = py_cpmfem_layers(
    layers=layers,
    VOXSIZE=0.0025,
    sizeX=1.0, 
    sizeY=1.0,
    scenario="monolayer_on_fiber",
    NRINC=1000,
    JZ_CMCM=0.2,
    JZ_FBFB=0.1,
    JZ_CMFB=0.3
)
```

## Cross-Layer Coupling Details

### Energy Calculation
The total energy change for a voxel move now includes cross-layer interactions:

```
ΔH = ΔH_contact + ΔH_volume + ΔH_focals + ΔH_syncytium + ΔH_nuclei + ΔH_crosslayer
```

Where `ΔH_crosslayer` considers interactions with voxels in adjacent layers (z±1).

### Coupling Parameters
- **`JZ_CMCM`**: Energy cost/benefit for CM-CM interactions across layers
- **`JZ_FBFB`**: Energy cost/benefit for FB-FB interactions across layers  
- **`JZ_CMFB`**: Energy cost/benefit for CM-FB interactions across layers

### Interaction Rules
- Voxels at position (x,y) in layer L interact with voxels at (x,y) in layers L±1
- Medium-voxel interactions use reduced coupling (0.1 × base coupling)
- Same-type interactions use full coupling strength
- Mixed-type interactions use the specified CM-FB coupling

## Performance Considerations

### Memory Usage
- Memory scales linearly with number of layers
- Each layer requires: `NVX × NVY × sizeof(VOX)` bytes
- Total memory: `NL × NVX × NVY × sizeof(VOX)`

### Computational Complexity
- Cross-layer energy calculations add ~2× overhead per voxel move
- Simulation time scales approximately linearly with number of layers
- Recommended: 2-5 layers for optimal performance

## Migration Guide

### From Single-Layer to Multilayer

**Old code:**
```python
types, ctags, fibers, contacts = py_cpmfem(
    NCX=10, NCY=10, PART=0.3,
    VOXSIZE=0.0025, sizeX=1.0, sizeY=1.0,
    scenario="monolayer_on_fiber", NRINC=1000
)
```

**New code (equivalent):**
```python
layers = [{"type": "basic", "NCX": 10, "NCY": 10, "PART": 0.3}]
types_list, ctags, fibers, contacts = py_cpmfem_layers(
    layers=layers,
    VOXSIZE=0.0025, sizeX=1.0, sizeY=1.0,
    scenario="monolayer_on_fiber", NRINC=1000
)
# Note: ctags, fibers, contacts now have shape (1, NVY, NVX)
```

### Accessing Results
```python
# For single layer (old way)
layer_0_ctags = ctags[0]  # Shape: (NVY, NVX)

# For multiple layers
for i in range(len(types_list)):
    layer_i_ctags = ctags[i]  # Shape: (NVY, NVX)
    layer_i_types = types_list[i]  # Shape: (NCX*NCY,)
```

## Troubleshooting

### Common Issues

1. **Memory errors with many layers**
   - Reduce `sizeX`, `sizeY` or `VOXSIZE` to decrease voxel count
   - Use fewer layers or shorter simulations

2. **Slow performance**
   - Reduce `NRINC` for testing
   - Use smaller domain sizes
   - Consider fewer layers

3. **Unexpected cell distributions**
   - Check `PART_matrix` values are between 0.0 and 1.0
   - Verify layer configurations are correct
   - Adjust cross-layer coupling parameters

### Debugging Tips

```python
# Check layer configurations
for i, layer in enumerate(layers):
    print(f"Layer {i}: {layer}")

# Verify output shapes
print(f"ctags shape: {ctags.shape}")
print(f"Expected: ({len(layers)}, NVY, NVX)")

# Check cell type distributions
for i, types in enumerate(types_list):
    unique, counts = np.unique(types, return_counts=True)
    print(f"Layer {i} cell types: {dict(zip(unique, counts))}")
```

## Future Enhancements

- **Ready layer initialization**: Support for predefined cell configurations
- **3D visualization**: Tools for visualizing multilayer results
- **Advanced coupling**: Distance-dependent cross-layer interactions
- **Layer-specific parameters**: Different energy parameters per layer
- **Dynamic layer addition**: Add/remove layers during simulation

## Technical Implementation

### C API Changes
- Added `NL` parameter to `cpmfem()` function
- Extended energy calculations with cross-layer terms
- Modified output aggregation for multiple layers

### Python Interface
- New `py_cpmfem_layers()` function
- Backward-compatible `py_cpmfem()` function
- Enhanced error handling and validation

### Memory Management
- Automatic allocation/deallocation of 3D arrays
- Proper cleanup of temporary structures
- Efficient data copying between layers
