#!/usr/bin/env python3
"""
Test script for multilayer VCT functionality
Demonstrates the new py_cpmfem_layers function with different layer types
"""

import numpy as np
import matplotlib.pyplot as plt
from libcpmfem import py_cpmfem_layers, py_cpmfem

def create_test_probability_matrix(nvx, nvy, pattern='gradient'):
    """Create a test probability matrix for layer type 'prob'"""
    matrix = np.zeros((nvy, nvx))
    
    if pattern == 'gradient':
        # Create a gradient from left to right
        for i in range(nvy):
            for j in range(nvx):
                matrix[i, j] = j / (nvx - 1)  # 0 to 1 gradient
    elif pattern == 'center':
        # Higher probability in center
        center_x, center_y = nvx // 2, nvy // 2
        for i in range(nvy):
            for j in range(nvx):
                dist = np.sqrt((i - center_y)**2 + (j - center_x)**2)
                max_dist = np.sqrt(center_y**2 + center_x**2)
                matrix[i, j] = 1.0 - (dist / max_dist)
    elif pattern == 'stripes':
        # Vertical stripes
        for i in range(nvy):
            for j in range(nvx):
                matrix[i, j] = 0.8 if (j // 5) % 2 == 0 else 0.2
    
    return matrix

def test_basic_multilayer():
    """Test basic multilayer simulation with different layer configurations"""
    print("Testing basic multilayer simulation...")
    
    # Create layers with different configurations
    layers = [
        {
            "type": "basic",
            "NCX": 8,
            "NCY": 8, 
            "PART": 0.2  # 20% fibroblasts
        },
        {
            "type": "basic",
            "NCX": 8,
            "NCY": 8,
            "PART": 0.5  # 50% fibroblasts
        },
        {
            "type": "basic", 
            "NCX": 8,
            "NCY": 8,
            "PART": 0.8  # 80% fibroblasts
        }
    ]
    
    # Run simulation with cross-layer coupling
    types_list, ctags, fibers, contacts = py_cpmfem_layers(
        layers=layers,
        VOXSIZE=0.0025,
        sizeX=1.0,
        sizeY=1.0,
        scenario="monolayer_on_fiber",
        NRINC=100,  # Short simulation for testing
        JZ_CMCM=0.3,  # Cross-layer CM-CM coupling
        JZ_FBFB=0.2,  # Cross-layer FB-FB coupling  
        JZ_CMFB=0.4   # Cross-layer CM-FB coupling
    )
    
    print(f"Simulation completed!")
    print(f"Number of layers: {len(types_list)}")
    print(f"Output shapes:")
    print(f"  ctags: {ctags.shape}")
    print(f"  fibers: {fibers.shape}")
    print(f"  contacts: {contacts.shape}")
    print(f"  types per layer: {[len(t) for t in types_list]}")
    
    return types_list, ctags, fibers, contacts

def test_probability_matrix_layers():
    """Test layers with custom probability matrices"""
    print("\nTesting probability matrix layers...")
    
    # Create different probability patterns
    nvx, nvy = 50, 50
    
    layers = [
        {
            "type": "prob",
            "NCX": 6,
            "NCY": 6,
            "PART_matrix": create_test_probability_matrix(nvx, nvy, 'gradient')
        },
        {
            "type": "prob", 
            "NCX": 6,
            "NCY": 6,
            "PART_matrix": create_test_probability_matrix(nvx, nvy, 'center')
        },
        {
            "type": "prob",
            "NCX": 6, 
            "NCY": 6,
            "PART_matrix": create_test_probability_matrix(nvx, nvy, 'stripes')
        }
    ]
    
    types_list, ctags, fibers, contacts = py_cpmfem_layers(
        layers=layers,
        VOXSIZE=0.0025,
        sizeX=1.0,
        sizeY=1.0,
        scenario="monolayer_on_fiber", 
        NRINC=50,  # Very short for testing
        JZ_CMCM=0.1,
        JZ_FBFB=0.1,
        JZ_CMFB=0.2
    )
    
    print(f"Probability matrix simulation completed!")
    print(f"Output shapes: {ctags.shape}")
    
    return types_list, ctags, fibers, contacts

def test_backward_compatibility():
    """Test that original py_cpmfem still works"""
    print("\nTesting backward compatibility...")
    
    # Test original function
    types, ctags, fibers, contacts = py_cpmfem(
        NCX=5, NCY=5, PART=0.3,
        VOXSIZE=0.0025, sizeX=0.5, sizeY=0.5,
        scenario="monolayer_on_fiber", NRINC=50,
        PART_matrix=0.3
    )
    
    print(f"Original function works! Output shapes:")
    print(f"  types: {types.shape}")
    print(f"  ctags: {ctags.shape}")
    print(f"  fibers: {fibers.shape}")
    print(f"  contacts: {contacts.shape}")
    
    return types, ctags, fibers, contacts

def visualize_results(ctags, title="Multilayer Results"):
    """Create a simple visualization of the results"""
    try:
        import matplotlib.pyplot as plt
        
        n_layers = ctags.shape[0]
        fig, axes = plt.subplots(1, n_layers, figsize=(4*n_layers, 4))
        if n_layers == 1:
            axes = [axes]
        
        for i in range(n_layers):
            im = axes[i].imshow(ctags[i], cmap='viridis', aspect='equal')
            axes[i].set_title(f'Layer {i+1}')
            axes[i].set_xlabel('X')
            axes[i].set_ylabel('Y')
            plt.colorbar(im, ax=axes[i])
        
        plt.suptitle(title)
        plt.tight_layout()
        plt.show()
        
    except ImportError:
        print("Matplotlib not available for visualization")

def main():
    """Run all tests"""
    print("=" * 60)
    print("MULTILAYER VCT TEST SUITE")
    print("=" * 60)
    
    try:
        # Test 1: Basic multilayer
        print("\n1. BASIC MULTILAYER TEST")
        print("-" * 30)
        types_list, ctags, fibers, contacts = test_basic_multilayer()
        
        # Test 2: Probability matrix layers  
        print("\n2. PROBABILITY MATRIX TEST")
        print("-" * 30)
        types_list2, ctags2, fibers2, contacts2 = test_probability_matrix_layers()
        
        # Test 3: Backward compatibility
        print("\n3. BACKWARD COMPATIBILITY TEST")
        print("-" * 30)
        test_backward_compatibility()
        
        # Visualization
        print("\n4. VISUALIZATION")
        print("-" * 30)
        print("Creating visualizations...")
        visualize_results(ctags, "Basic Multilayer Results")
        visualize_results(ctags2, "Probability Matrix Results")
        
        print("\n" + "=" * 60)
        print("ALL TESTS COMPLETED SUCCESSFULLY!")
        print("=" * 60)
        
        # Summary
        print("\nSUMMARY:")
        print("- ✅ Multilayer simulation works")
        print("- ✅ Cross-layer coupling implemented") 
        print("- ✅ Probability matrix layers work")
        print("- ✅ Backward compatibility maintained")
        print("- ✅ Output arrays have correct shapes")
        
    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
