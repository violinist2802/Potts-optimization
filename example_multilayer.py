#!/usr/bin/env python3
"""
Simple example demonstrating multilayer VCT functionality
"""

import numpy as np
from libcpmfem import py_cpmfem_layers

def create_gradient_matrix(nvx, nvy):
    """Create a left-to-right gradient probability matrix"""
    matrix = np.zeros((nvy, nvx))
    for i in range(nvy):
        for j in range(nvx):
            matrix[i, j] = j / (nvx - 1)  # 0 to 1 gradient
    return matrix

def create_center_matrix(nvx, nvy):
    """Create a center-focused probability matrix"""
    matrix = np.zeros((nvy, nvx))
    center_x, center_y = nvx // 2, nvy // 2
    for i in range(nvy):
        for j in range(nvx):
            dist = np.sqrt((i - center_y)**2 + (j - center_x)**2)
            max_dist = np.sqrt(center_y**2 + center_x**2)
            matrix[i, j] = max(0, 1.0 - (dist / max_dist))
    return matrix

def main():
    print("Multilayer VCT Example")
    print("=" * 40)
    
    # Create 3 layers with different configurations
    layers = [
        {
            "type": "basic",
            "NCX": 6,
            "NCY": 6,
            "PART": 0.2  # 20% fibroblasts
        },
        {
            "type": "prob", 
            "NCX": 6,
            "NCY": 6,
            "PART_matrix": create_gradient_matrix(40, 40)
        },
        {
            "type": "prob",
            "NCX": 6, 
            "NCY": 6,
            "PART_matrix": create_center_matrix(40, 40)
        }
    ]
    
    print("Layer configurations:")
    for i, layer in enumerate(layers):
        print(f"  Layer {i+1}: {layer['type']} with {layer.get('PART', 'custom matrix')}")
    
    # Run multilayer simulation
    print("\nRunning simulation...")
    types_list, ctags, fibers, contacts = py_cpmfem_layers(
        layers=layers,
        VOXSIZE=0.0025,
        sizeX=1.0,
        sizeY=1.0,
        scenario="monolayer_on_fiber",
        NRINC=200,  # Short simulation for demo
        JZ_CMCM=0.3,  # Cross-layer CM-CM coupling
        JZ_FBFB=0.2,  # Cross-layer FB-FB coupling
        JZ_CMFB=0.4   # Cross-layer CM-FB coupling
    )
    
    print("Simulation completed!")
    print(f"Output shapes:")
    print(f"  ctags: {ctags.shape}")
    print(f"  fibers: {fibers.shape}")
    print(f"  contacts: {contacts.shape}")
    print(f"  types per layer: {[len(t) for t in types_list]}")
    
    # Analyze results
    print("\nResults analysis:")
    for i in range(len(layers)):
        layer_ctags = ctags[i]
        layer_types = types_list[i]
        
        # Count cell types
        unique_types, counts = np.unique(layer_types, return_counts=True)
        type_counts = dict(zip(unique_types, counts))
        
        # Count occupied voxels
        occupied_voxels = np.sum(layer_ctags > 0)
        total_voxels = layer_ctags.size
        
        print(f"  Layer {i+1}:")
        print(f"    Cell types: {type_counts}")
        print(f"    Occupied voxels: {occupied_voxels}/{total_voxels} ({100*occupied_voxels/total_voxels:.1f}%)")
    
    print("\nExample completed successfully!")
    return types_list, ctags, fibers, contacts

if __name__ == "__main__":
    main()
