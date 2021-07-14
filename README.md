# Color Guided Depth Recovery from RGB-D Data Using an Adaptive Auto-regressive Model

This is the C implementation of the color guided AR model presented in "Color-guided Depth Recovery from RGB-D Data 
Using an Adaptive Auto-regressive Model" in TIP 2014. The original author provided MATLAB source code can be downloaded
here:http://cs.tju.edu.cn/faculty/likun/projects/depth_recovery/index.htm. Our implementation can produce the same results
as theirs. The difference is that we implement the method with C language and sparse matrix which is much faster and can handle
much larger depth maps.

*********************************************************************

**How to use**

run MexFile.m first to compile the cpp file, then run main.m to perform experiments.

************************************************************************
**Related Work:**
1. *"Robust Color Guided Depth Map Restoration.", Wei Liu, Xiaogang Chen, Jie Yang and Qiang Wu. In IEEE Transactions on Image Processing, 26(1), 315-327.* [Code](https://github.com/wliusjtu/Robust-Color-Guided-Depth-Map-Restoration)
2. *"A generalized framework for edge-preserving and structure-preserving image smoothing." Wei Liu, Pingping Zhang, Yinjie Lei, Xiaolin Huang, Jie Yang, and Michael Ng, IEEE Transactions on Transactions on Pattern Analysis and Machine Intelligence (2021).* [Code](https://github.com/wliusjtu/Generalized-Smoothing-Framework)

