#!/usr/bin/env python3
"""
Build script for multilayer VCT implementation
Compiles the updated C code and Cython extensions
"""

import os
import sys
import subprocess
import shutil
from pathlib import Path

def run_command(cmd, cwd=None):
    """Run a command and return success status"""
    print(f"Running: {cmd}")
    try:
        result = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error: {result.stderr}")
            return False
        print(f"Success: {result.stdout}")
        return True
    except Exception as e:
        print(f"Exception: {e}")
        return False

def build_c_library():
    """Build the C library with multilayer support"""
    print("Building C library...")
    
    vct_dir = Path("VCT")
    if not vct_dir.exists():
        print("Error: VCT directory not found")
        return False
    
    # Change to VCT directory
    os.chdir(vct_dir)
    
    # Clean previous builds
    print("Cleaning previous builds...")
    for ext in ['.o', '.a']:
        for file in Path('.').glob(f'*{ext}'):
            file.unlink()
    
    # Compile C sources
    c_sources = [
        'init.c', 'cellmoves.c', 'CPM_dH.c', 'CM.c', 
        'mylib.c', 'mt.c', 'write.c', 'cpmfem.c'
    ]
    
    compile_cmd = "gcc -c -O2 -fPIC"
    for source in c_sources:
        if Path(source).exists():
            cmd = f"{compile_cmd} {source}"
            if not run_command(cmd):
                return False
        else:
            print(f"Warning: {source} not found")
    
    # Create static library
    print("Creating static library...")
    cmd = "ar rcs libcpmfem.a *.o"
    if not run_command(cmd):
        return False
    
    print("C library built successfully!")
    return True

def build_cython_extension():
    """Build the Cython extension"""
    print("Building Cython extension...")
    
    # Go back to root directory
    os.chdir("..")
    
    # Check if setup.py exists
    if not Path("setup.py").exists():
        print("Error: setup.py not found")
        return False
    
    # Build Cython extension
    cmd = "python setup.py build_ext --inplace"
    if not run_command(cmd):
        return False
    
    print("Cython extension built successfully!")
    return True

def test_import():
    """Test if the module can be imported"""
    print("Testing module import...")
    
    try:
        import libcpmfem
        print("✅ libcpmfem imported successfully")
        
        # Test if new function exists
        if hasattr(libcpmfem, 'py_cpmfem_layers'):
            print("✅ py_cpmfem_layers function found")
        else:
            print("❌ py_cpmfem_layers function not found")
            return False
            
        if hasattr(libcpmfem, 'py_cpmfem'):
            print("✅ py_cpmfem function found (backward compatibility)")
        else:
            print("❌ py_cpmfem function not found")
            return False
            
        return True
        
    except ImportError as e:
        print(f"❌ Import failed: {e}")
        return False
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        return False

def run_tests():
    """Run basic tests"""
    print("Running basic tests...")
    
    try:
        from libcpmfem import py_cpmfem_layers
        import numpy as np
        
        # Simple test
        layers = [
            {"type": "basic", "NCX": 3, "NCY": 3, "PART": 0.3}
        ]
        
        types_list, ctags, fibers, contacts = py_cpmfem_layers(
            layers=layers,
            VOXSIZE=0.0025,
            sizeX=0.5,
            sizeY=0.5,
            scenario="monolayer_on_fiber",
            NRINC=10  # Very short test
        )
        
        print("✅ Basic multilayer test passed")
        print(f"   Output shapes: ctags={ctags.shape}, fibers={fibers.shape}")
        return True
        
    except Exception as e:
        print(f"❌ Test failed: {e}")
        return False

def main():
    """Main build process"""
    print("=" * 60)
    print("MULTILAYER VCT BUILD SCRIPT")
    print("=" * 60)
    
    # Check if we're in the right directory
    if not Path("VCT").exists():
        print("Error: Please run this script from the pyVCT_MRI root directory")
        return False
    
    success = True
    
    # Step 1: Build C library
    print("\n1. BUILDING C LIBRARY")
    print("-" * 30)
    if not build_c_library():
        print("❌ C library build failed")
        success = False
    
    # Step 2: Build Cython extension
    print("\n2. BUILDING CYTHON EXTENSION")
    print("-" * 30)
    if not build_cython_extension():
        print("❌ Cython extension build failed")
        success = False
    
    # Step 3: Test import
    print("\n3. TESTING IMPORT")
    print("-" * 30)
    if not test_import():
        print("❌ Import test failed")
        success = False
    
    # Step 4: Run tests
    print("\n4. RUNNING TESTS")
    print("-" * 30)
    if not run_tests():
        print("❌ Basic tests failed")
        success = False
    
    # Summary
    print("\n" + "=" * 60)
    if success:
        print("✅ BUILD COMPLETED SUCCESSFULLY!")
        print("✅ Multilayer VCT is ready to use")
        print("\nNext steps:")
        print("1. Run: python example_multilayer.py")
        print("2. Run: python test_multilayer.py")
        print("3. Read: MULTILAYER_README.md")
    else:
        print("❌ BUILD FAILED!")
        print("Check the error messages above")
    
    print("=" * 60)
    return success

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
