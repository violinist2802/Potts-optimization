# pyVCT
Virtual Cardiac Tissue Model – A Cellular Potts Model for cardiac monolayers that reproduces fibrotic patterns

## Build instructions
`git clone git@github.com:violinist2802/Potts-optimization.git`  
`cd VCT`  
`make`  
`cd ..`  
`python3 setup.py build_ext --inplace`    
## Usage
`import pyVCT`  
`types, c_tags, fibers, contacts = pyVCT.py_cpmfem(26, 26, 0.5, 0.0025, 1, 1, 'monolayer_without_fiber', 2001)`
