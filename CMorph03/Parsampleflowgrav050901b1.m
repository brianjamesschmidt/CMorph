% Parameter file for Cmorph03.m
% Filename of avi stack to analyze
% The program will append the avi extension automatically.
inputavi='sampleflowgrav3um60X500fps050901b1';

% Enter any additional identifiers you want appended to the name of the
% output avi's and centroid matrix files
outputappend='track';

% Set the image feature desired to track based on.  By default, the program
% thresholds the gradient of the image to detect edges.  If track gradient
% is false, the program will invert and threshold the image itself.
trackgradient=true;

% Input a point near the center of the cell to be tracked in the first
% frame.  Pick this out of ImageJ, top/bottom matrix conventions accounted
% for in program.
horpointo=489;
verpointo=94;

% Specify whether the program should perform deconvolution while tracking
deconvolve=false;
psffile='measure_your_psf_and_specify_it_here.tif';

% Specify whether to correct for gain errors.  Each should be in a *.mat
% file saved as Slope and Intercept.  Program will append *.mat.
gaincorrect=false;
slopefile='measure_your_slope_and_specify_it_here';
interceptfile='measure_your_intercept_and_specify_it_here';

% Enter the beginning frame in the avi to use for tracking
% Input NaN, and the program will assume 1.
startframe=NaN;

% Enter the ending frame to use for image tracking
% Input NaN, and the program will stop at the end or when the object
% goes off the screen.
endframe=NaN;

% Enter the height of the image data (not including headers) in pixels
% Right now the program assumes there may be a header in the top of the
% vertical direction but no borders around the left and right edges
dheight=240;

% If there is header data, indicate whether it is appended to the top of
% the images or bottom. " headertop = true " or " headertop = false "
headertop=true;

% Enter whether you wish to reduce the program window when saving the 
% output images to save on file size.  If so, the tracking algorithm will
% limit the output to +/- 3 bead diameters. " reducewindow = true "
% or " reducewindow = false "
reducewindow=false;

% Enter the approximate diameter of the bead in micrometers.
beadsize=3;

% Image resolution in pixels/um:
xres=8.2;
yres=8.2;

% Enter record speed in fps
fps=500;

% Compression setting to use for output.  Pick 'indeo5' or 'none'.
compressionset='none';

% Indicate whether you would like encode shape data in the
% the generated movie file.  If set to true, the program will show
% the major and minor axis with lines and adjust object coloring according
% to deformation.  If false, it and coloring is on, it will highlight the 
% object with orange so it is easy to follow
showshapedata=false;

% Boolean variable indicating whether to apply a "blurring" filter while
% picking out edges.  Sometimes helpful for weak signal situations or with
% ragged edges.
% Blursize is the nXn mask used for blurring.  Choose an odd number if
% possible to keep the tracked area in the center of the object.  If an
% even number is chosen, however, the tracking will still be good since
% frame-to-frame differences are important.
blur=false;
filtersizeimage=1;
filtersizegradient=1;

% In case the center of the bead registers as a bright ring, allow for
% these pixels to be zeroed out.  Set the area limits here in pixels.  
% When set to "true," this option may slow the program a bit.
middlering=true;
minringarea=(beadsize/2*xres)^2*pi*.02;
maxringarea=(beadsize/2*xres)^2*pi*.8;

% Maximum number of times at each threshold value to look for bead in image.
% If set greater than 1, the program will look at additional data points in
% the vicinity of the first guess.
numberattemptsmax=1;

% Indicate whether separate binary image or overlay is desired.
binaryout=false;
