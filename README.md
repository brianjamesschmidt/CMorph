CMorph
=============
Tracking algorithm for cell motion and morphology (MATLAB)
-------------

This is a multiparticle version of the CMorph algorithm previously developed. The algorithm takes the gradient of the image, and then applies decreasing thresholds and fill operations to attempt to identify objects. One it has found the objects, it attempts to pair them up with objects in the previous frame. It was designed to track cells and microbeads in flow cell assays, so it is designed to be biased to look for a flow direction.This is a basic cell tracker for use on uncompressed AVI files, ideal for use with videos where the image of the cell overlap much frame-to-frame. If you use it to generate data for publication please cite the conference paper in the readme file, which you might want to read anyway if you want to use. Sample video to track and paper available at lab website: bme.virginia.edu/lawrence.

Instructions:
Download the tracking algorithm and sample parameter file. 
The parameter file is set up for the sample movie, 
or use your own uncompressed avi and adjust the parameter file.
Place all 3 in your Matlab root directory.
Type: CMorph03('Parsampleflowgrav050901b1').

Reference available at IEEXplore: 

Schmidt, B.J., C.D. Paschall, W.H. Guilford, and M.B. Lawrence, "High-Resolution Optical Tracking to Identify Adhesive Events in Vitro." 
Conference Proceedings of the Asilomar Conference on Signals, Systems, and Computers, Pacific Grove, CA, November 4-7, 2007.  p. 1856-1860.

Requirements:
MATLAB release		developed with MATLAB 7.2 (R2006a)
Other requirements	Image Processing Toolbox

You might want to try MCShape first in case that works better.

Release history
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Release Date      Notes
03.0		5.15.2009	-Last update to package as distributed on MATLAB Central.
