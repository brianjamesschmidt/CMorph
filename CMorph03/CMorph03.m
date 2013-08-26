function []=CMorph03(parameterfilename)
%  CMorph03(parameter filename as 'parameterfilename'.m)
%  Cell tracking program originally coded by Brian Schmidt at the 
%  University of Virginia, beginning in May 2005.
%  Two operation modes.  In gradient mode, This program uses the 
%  contrast in surrounding objects (Becke rings if defocused) to attempt 
%  to locate the edge.  It then uses the centroid of the generated shape to
%  track position and records morphologic properties.
%  Results are output into two files: one is a *.mat file with the tracking
%  results and the other is a video file that shows the tracked region.
%  *.mat file formating:
%  Column 1: Apparent object radius
%  Column 2: Time
%  Column 3: Centroid, X-Coordinate
%  Column 4: Centroid, Y-Coordinate
%  Column 5: Morphology data: bounding box x-width
%  Column 6: Morphology data: bounding box y-height
%  Column 7: Morphology data: orientation angle (degrees)
%  Column 8: Morphology data: object length along orientation angle axis
%  Column 9: Morphology data: Object width along line perpendicular to
%           orientation angle axis through centroid
%  The Shape Strain Index used to color the object is the ratio of
%  column 9 : column 8. 
%  The video should be watched to verify the algorithm correctly identified
%  the object.  This algorithm is ideal for low noise situations where 
%  there is substantial overlap of the object in subsequent images.
%  -----REVISION HISTORY-----
%  V03:    5/15/09:   Fixed the guessing strategy for when a an object was
%                     not detected with the previous object location.
%                     Thanks to Harrison Prentice-Mott for noticing this!
%  V02:    9/18/08:   Adjusted output movie coloring so the object would 
%                     be orange if showshapedata (showline in V01) is 
%                     set to false.
%  V01:   11/02/07:   Changed algorithm name to reflect recent changes in
%                     functionality.  Replaced gray 235 cmaps with gray
%                     256 so output appearance better matched input.  Added
%                     in the ability to automatically correct for linear
%                     gain errors in the camera pixel elements given the 
%                     slope and intercept of each pixel value.
%  V15:   6/05/07:    Added in additional morpholical operations to
%                     identify the major & minor object axis.  Also
%                     updated the search algorithm for a new point in
%                     the object if the attempt using the
%                     previous centroid fails.  Added in an intelligent
%                     selection of the ending frame so it can be
%                     optionally specified.
%  V14:   3/06/07:    Added an optional deconvolution step to improve
%                     accuracy of tracked object dimensions.  Corrected
%                     an error introduced in the avi write portion of V13
%                     that prevented ImageJ from properly recognizing the
%                     format of RGB tracked overlays.  Also sped the code
%                     up with the Matlab profiler.
%  V13:   2/12/07:    Added the object's length and width to the output.  
%                     Also incorporated the option to output a binary movie
%                     series depicting the shape.
%  V12:   12/10/06:   Limited the program to perform image calculations
%                     only in the vicinity of the last bead location to 
%                     speed up execution.  Also allow the user to save a
%                     small horizontal patch of the output file to save on
%                     disk space.
%  V11:   12/08/06:   Revamped the program to handle separate X and Y
%                     resolutions for handling video from analog sources
%                     that have been deinterlaced.  Also added an
%                     intensity-based tracking option.
%  V10:   4/14/06:    Updated the filtering scheme to add a blurring filter
%                     to the original image to help mask noise and switched
%                     the filter of the gradient image to a median to
%                     better preserve the edges.
%  V09:   9/01/05:    Added a backup of the center position every 10 frames
%                     in case program is terminated early.  Also switched
%                     the program to read in movies frame-by-frame rather
%                     than generate a big matrix to store all the frames
%                     in.
%  V06-08:8/19/05:    Improved the code to more accurately track position
%                     and eliminate edge switching phenomenon.
%  V05:   8/17/05:    Added an "intelligent" search for threshold 
%                     intensity.                  
%  V04:   8/01/05:    Ran into some bead data with areas of high intensity 
%                     gradients in the middle of the beads.  Added an
%                     algorithm to follow the outer gradient.
%  V03:   7/18/05:    Modified the parameter file to specify whether to
%                     apply a blurring filter and set the threshold value.
%                     Also modified the output file to highlight the bead
%                     in color and save an estimated bead radius.
%  V02:   7/11/05:    Made several revisions including writing distances in
%                     µm and times to the output file instead of pixels.
%  V01:   5/12/05:    Original Version.  Special thanks to Yinbo Li for
%                     brainstorming and help with some of the image 
%                     processing features.

tic
input_file=strcat(parameterfilename);

% Read the input parameters into the workspace
eval(input_file);

% Compute some data both the resolution that will be useful
maxres=max(xres,yres);
minres=min(xres,yres);
resratio=minres/maxres;

% AVI output file to write the overlay to.
avioutfile=strcat(inputavi,outputappend,'overlay.avi');

% The filename of the centroid location is automatically appended with .mat
% by matlab.
centertrackerfile=strcat(inputavi,outputappend);

inputavi=strcat(inputavi,'.avi');

% Assume the beginning frame for tracking is 1 if the user does not
% specify.  If the user does not specify an ending frame, the program
% will automatically stop when the object reaches the edge.
if isnan(startframe) == 1
    startframe = 1;
end

if isnan(endframe) == 1
    fileinfo = aviinfo(inputavi);
    endframe = fileinfo.NumFrames;
end

% First read in the first slice to get the image size.
nmovieframes=endframe-startframe+1;

Tempframe1 = aviread(inputavi,startframe);
[height,width]=size(getfield(Tempframe1,'cdata'));

% Update the estimate for the initial particle center to account for the
% data section on the file.
nhorpointo=horpointo;
if headertop ~= false
    nverpointo=verpointo-(height-dheight);
else
    nverpointo=verpointo;
end

% Initialize a matrix to track the centroids
Centertracker=zeros(1,9);

% Initialize the avi used to verify the centroid tracker.
aviobj = avifile(avioutfile);
aviobj.compression=compressionset;
aviobj.fps=30;
aviobj.keyframe=30;
aviobj.quality=100;

% Minimum number of pixels that must be detected in frame for pocessing
borderareamin=pi*beadsize*minres;

% Set a threshold size to make sure the entire background isn't being picked
% up.
dthresharea=min(dheight*width*.5,3*pi*beadsize^2*.25*xres*yres);

% Initialize the area check
lastarea=round((beadsize/2*maxres)^2*pi*resratio);

% Set up the filter for mean filtering.
if blur == 1
    h = ones(filtersizeimage,filtersizeimage) / filtersizeimage^2;
end

if deconvolve == true
    % Read in the image to deconvolve
    PSF=imread(psffile);
    PSF=double(PSF);
    % Scale the PSF so it has a maximum of 1
    PSF=PSF./(max(max(PSF)));
end

if gaincorrect==true
    % Read in the slope and intercept files
    Temp=open(strcat(interceptfile,'.mat'));
    Intercept=Temp.Intercept;
    clear Temp
    Temp=open(strcat(slopefile,'.mat'));
    Slope=Temp.Slope;
    clear Temp
end

counter=1;
inframe=true;
% Iterate until all of the frames have been processed or the object leaves
% the screen.
while counter <= nmovieframes && inframe == true
    
    Tempframe1 = aviread(inputavi,counter+startframe-1);
    Tempframe1 = (getfield(Tempframe1,'cdata'));
    % Pick the image data out of the matrix without the header data
    if headertop ~= false
        Movmatrix=Tempframe1(height-dheight+1:height,1:width);
    else
        Movmatrix=Tempframe1(1:dheight,1:width);
    end
    
    % Apply the gain error correction if instructed to.
    if gaincorrect==true
        Movmatrix=double(Movmatrix);
        Movmatrix=Slope.*Movmatrix+Intercept;
        if headertop ~= false
            Tempframe1(height-dheight+1:height,1:width)=Movmatrix;
        else
            Tempframe1(1:dheight,1:width)=Movmatrix;
        end   
    end
                    
    % The area is used as a flag to check that a cell has been detected.
    % Initialize it to zero for the current frame. 
    STATS.Area=0;

    % Cell size flag indicates whether the area of the detected cell in the
    % current frame has passed muster yet.
    cellsizegood=false;

    % Use the sum of the absolute value of the gradient in the x and
    % y-directions as an indication where image intensity is changing to 
    % detect edges.  Note that this value is arbitrarily assigned to the 
    % upper-left hand corner of the 3 pixels involved in the subtraction, 
    % but this is consistently applied so the changes in position should 
    % still be accurate.
    Tempmatrix=Movmatrix;
    Tempmatrix=double(Tempmatrix);

    % For tracking, focus on the area near where the user specified the
    % bead location to speed the computations.  Just read in pixels within 
    % the nearest +/- 2 bead sizes of the last centroid.  These are
    % relative to image without header.
    leftpixel=max(round(nhorpointo-2*beadsize*xres),1);
    rightpixel=min(round(nhorpointo+2*beadsize*xres),width);
    toppixel=max(round(nverpointo-2*beadsize*yres),1);
    bottompixel=min(round(nverpointo+2*beadsize*yres),dheight);      
        
    Tempmatrix=Tempmatrix(toppixel:bottompixel,leftpixel:rightpixel);
    
    % Use the Lucy-Richardson deconvolution algorithm deconvolution if
    % deconvolution is specified.  
    if deconvolve == true
        J = deconvlucy(Tempmatrix,PSF);
        maxj=max(max(J));
        maxtemp=max(max(Tempmatrix));
        % Scale this image to the same dynamic range as the original image.
        Tempmatrix=(J.*(maxtemp/maxj));
        % Also update the original movie with the deconvolved image so the
        % improvements will be evident in the tracked data.
        if headertop ~= false
            Tempframe1(height-dheight+toppixel:height-dheight+bottompixel,leftpixel:rightpixel)=uint8(Tempmatrix);
        else
            Tempframe1(toppixel:bottompixel,leftpixel:rightpixel)=uint8(Tempmatrix);
        end
    end

    % Denoising the image with a mean filter may help in some processing 
    % situations.
    if blur == 1
        Tempmatrix = imfilter(Tempmatrix,h,'replicate');
    end
    
    if trackgradient == false
          tempmaxval=max(max(Tempmatrix));
          tempminval=min(min(Tempmatrix));
          % If the intensity will be used to track, invert the image to
          % take advantage of the dark rings at the beads edges when
          % thresholding.
          FT=tempminval+tempmaxval-Tempmatrix;
    else
        % Otherwise take the gradient for tracking
        [FX,FY] = gradient(Tempmatrix);
        FT=sqrt((FX).^2+(FY).^2);
    end
    
    % Also try denoising by median filtering the gradient image if the 
    % user indicates denoising is necessary.
    if blur == 1
        FT = medfilt2(FT,[filtersizegradient filtersizegradient]);
    end
    
    % Scale the gradient image based on the maximum intensity.
    tempmaxval=max(max(FT));
    tempminval=min(min(FT));
    FT=(255/(tempmaxval-tempminval).*(FT-tempminval));
    FT=uint8(FT);

    % Find the threshold value for which there could potentially be enough
    % pixels to form a ring with the given diameter.
    Pixvals=sort(FT(:),1,'descend');
    % Start with the first threshold value that potentially gives and object
    % large enough to be the microbead.
    Pixvals=Pixvals(round(borderareamin):end);
    % Remove degenerate pixel values from the hitlist.  Begin by
    % counting the number of nonzero, nonrepeating elements
    [pixvalslength,pixvalswidth]=size(Pixvals);
    temppixval=-1;
    pixvalscounter2=0;
    for pixvalscounter1 = 1 : pixvalslength
        if (temppixval ~= Pixvals(pixvalscounter1)) && (Pixvals(pixvalscounter1)>0)
            temppixval=Pixvals(pixvalscounter1);
            pixvalscounter2=pixvalscounter2+1;
            Pixvalssel(pixvalscounter2,1)=temppixval;
        end
    end
    threshvalind=1;
    
    while threshvalind <= pixvalscounter2 && cellsizegood == false
        
        % Number of attempts at the current threshold value
        numberattempts=1;
        threshval=Pixvalssel(threshvalind);
        BW = ~(im2bw(FT,double(threshval-1)/255));
        % Iterate with current threshold value numberattemptsmax times to 
        % see if a good pixel can be selected that picks out the bead.
        while cellsizegood==false && numberattempts <= numberattemptsmax

            % Badloop is a logical check on the current loop with the 
            % current threshold value.  When it becomes true subsequent 
            % activities in the loop for the current threshold value stop  
            % and a new iteration is attempted.
            badloop=false;
            
            % If this isn't the first attempt at processing the frame,
            % search the local environment for the object
            if numberattempts < 2
                hortemprand=0;
                vertemprand=0;
            elseif numberattempts < 3 && counter > 2
                % First, take a guess based on the motion in the previous
                % frames.
                hortemprand=round((Centertracker(counter-1,3)-Centertracker(counter-2,3))*xres);
                vertemprand=round((Centertracker(counter-1,4)-Centertracker(counter-2,4))*yres);
                % Check to make sure these fall within the subsection read
                % in
                if hortemprand + nhorpointo >= rightpixel
                    hortemprand = rightpixel - 1 - nhorpointo;
                elseif hortemprand + nhorpointo <= leftpixel
                    hortemprand = leftpixel + 1 - nhorpointo;
                end
                if vertemprand + nverpointo >= bottompixel
                    vertemprand = bottompixel - 1 - nverpointo;
                elseif vertemprand + nverpointo <= toppixel
                    vertemprand = toppixel + 1 - nverpointo;
                end 
            else
                % Otherwise, search inside the local environment for the object
                hortemprand=round(randn*(rightpixel-leftpixel-2)/4);
                vertemprand=round(randn*(bottompixel-toppixel-2)/4);
                % Check to make sure these fall within the subsection read
                % in
                if hortemprand + nhorpointo >= rightpixel
                    hortemprand = rightpixel - 1 - nhorpointo;
                elseif hortemprand + nhorpointo <= leftpixel
                    hortemprand = leftpixel + 1 - nhorpointo;
                end
                if vertemprand + nverpointo >= bottompixel
                    vertemprand = bottompixel - 1 - nverpointo;
                elseif vertemprand + nverpointo <= toppixel
                    vertemprand = toppixel + 1 - nverpointo;
                end                
            end

            % horpointg and verpointg are relative to the boundaries of the
            % subselection that comprises the Tempmatrix and FT matrices
            horpointg=round(nhorpointo+hortemprand)-leftpixel+1;
            verpointg=round(nverpointo+vertemprand)-toppixel+1;

            % First check whether point is nonzero
            % If a middlering is in the image, set the initial guess pixel 
            % to one.
            if (middlering == true && BW(verpointg,horpointg)==0)
                BW(verpointg,horpointg)=1;
            end
            % If there is no middle ring and the guess is zero then the
            % guess must be bad.
            if (BW(verpointg,horpointg)==0) && (middlering == false)
                badloop=true;
            else
                % Otherwise try detecting an area.
                BW2 = bwselect(BW,horpointg,verpointg,4);
                BW3=bwfill(BW2,'holes');
                % if intensity values are used to track, add back in the
                % edges of the ring, which were selected out by
                % this process
                if trackgradient == false
                    BW3=BW3+~BW;
                    BW3 = bwselect(BW3,horpointg,verpointg,4);
                    BW3=bwfill(BW3,'holes');
                end
                STATS.Area = sum(sum(BW3));
            end
        
            % If it passes the nonzero point check, try the small ring 
            % check to see whether a strong gradient in the object that 
            % might cause the tracking program to use it rather than the 
            % outer edge.
            if (badloop == false) && (middlering == true) && ...
                    (STATS.Area > minringarea) && (STATS.Area < maxringarea)
                % Use the algorithm to get the outside of the outer 
                % rings and blank them out in the gradient image.
                BW=~BW;
                BW2=bwfill(BW,horpointg,verpointg);
                BW3 = bwselect(BW2,horpointg,verpointg,4);
                FT=immultiply(FT,~BW3);
                BW = ~(im2bw(FT,double(threshval-1)/255));
                badloop=true;
            end
        
            % Check to make sure there is not too much variation in area from
            % the previous frame.  If there is, it is likely due to an invalid
            % object.  Low shutter times may result in elongation of the 
            % object with increasing velocity. 
            if STATS.Area < .7*lastarea || STATS.Area >=dthresharea
                badloop = true;
            end
        
            if badloop == true
                numberattempts = numberattempts+1;
            else
                cellsizegood=true;
            end

        end
        % Go to the next value in the threshold vector.
        threshvalind=threshvalind+1;
    end
    
    % Store the radius of the detected object assuming it is truly round
    % Broader equation for an ellipse used to determine in case xres ~=
    % yres.
    Centertracker(counter,1)=sqrt(STATS.Area/(pi*xres*yres));
    
    % If the edge was found, record the data
    if cellsizegood == true
        lastarea = STATS.Area;
        STATS = regionprops(uint8(BW3),'Centroid');
        % If this is the first frame, remember the x,y coordinate so
        % subsequent positions will be remembered relative to it.
        % Record all coordinates relative to the original image,
        % not the subimage used for identifying the bead area.
        if counter == 1
            horpointo=STATS.Centroid(1)+leftpixel-1;
            verpointo=STATS.Centroid(2)+toppixel-1;
        end
        %Update the horizontal center
        nhorpointo=STATS.Centroid(1)+leftpixel-1;
        % Update the vertical centroid
        nverpointo=STATS.Centroid(2)+toppixel-1;
        % Store the current position in Centertracker.
        % Time in column 2
        Centertracker(counter,2)=(counter-1)/fps;
        % X-Coordinate in column 3
        Centertracker(counter,3)=(nhorpointo)/xres;
        % Y-Coordinate in column 4
        Centertracker(counter,4)=(nverpointo)/yres;
        ASPECTSTATS = regionprops(uint8(BW3),'BoundingBox');
        % Object width in column 5
        Centertracker(counter,5)=ASPECTSTATS.BoundingBox(3)/xres;
        % Object height in column 6
        Centertracker(counter,6)=ASPECTSTATS.BoundingBox(4)/yres;       
        % Search for the object orientation in degrees for column 7
        % First get the perimeter image
        BWperimeterimage = bwperim(BW3,8);
        [Perimeterpixelsrow,Perimeterpixelscolumn,Perimeterpixelsvalue]=find(BWperimeterimage);
        Perimeterpixelscolumn=Perimeterpixelscolumn+leftpixel-1;
        Perimeterpixelsrow=Perimeterpixelsrow+toppixel-1;
        % Now find the orientation angle & pixels
        [majoraxisangle,majoredgepixelx1,majoredgepixely1,majoredgepixelx2,majoredgepixely2,majoraxislength] = findorientation(Perimeterpixelsrow,Perimeterpixelscolumn,nhorpointo,nverpointo,xres,yres);
        Centertracker(counter,7)=majoraxisangle;
        Centertracker(counter,8)=majoraxislength;
        % Find the distances along the minor axis
        % Define the minor axis to be perpendicular to the major axis
        [minoredgepixelx1,minoredgepixely1,minoredgepixelx2,minoredgepixely2,minoraxislength] = findedgepointsonline(Perimeterpixelsrow,Perimeterpixelscolumn,nhorpointo,nverpointo,xres,yres,majoraxisangle+90);
        Centertracker(counter,9) = minoraxislength;
        morphologicalstrainindex=majoraxislength/minoraxislength;
    else
        % If the frame was a failure, the area will still be zero.  Don't
        % change the vertical or horizontal test coordinates for the next
        % frame since this is the best estimate.
        Centertracker(counter,2)=(counter-1)/fps;
        Centertracker(counter,3)=NaN;
        Centertracker(counter,4)=NaN;
        Centertracker(counter,5)=NaN;
        Centertracker(counter,6)=NaN;
        Centertracker(counter,7)=NaN;
        Centertracker(counter,8)=NaN;
        Centertracker(counter,9)=NaN;
        morphologicalstrainindex=1;
        ASPECTSTATS = regionprops(uint8(BW3),'BoundingBox');
        BW3=ones(bottompixel-toppixel+1,rightpixel-leftpixel+1);
        previousthresh=255;
    end

    % Save an overlay between the detected areas and the movie.
    % Tempframe1: Direct read of current frame including header data, ind
    % Movmatrix: Selection of entire current frame without header data, ind
    % BW3: Binary image of detected bead area within small subselection
    if binaryout == false 
        % RGBim1 will serve to highlight the tracked area in color. 
        RGBim1=uint8(BW3);
        % Reference colors to a jet colormap
        cmapmax=256;
        cmap=jet(cmapmax);
        % Index color according to morphological strain index
        if morphologicalstrainindex >= 1.5
            colorindex=cmapmax;
        else
            colorindex=round((morphologicalstrainindex-1)*(cmapmax-1)/(1.5-1));
        end
        if showshapedata == true
            RGBim1=colorindex*RGBim1;
        else
            RGBim1=(cmapmax*.74)*RGBim1;
        end
        RGBim1=ind2rgb((RGBim1),cmap);
        % Add in lines to show the strain axes.
        if showshapedata == true && cellsizegood == true
            for strainaxiscounter = 1 : 2
                if strainaxiscounter == 1
                    edgepixely1=majoredgepixely1;
                    edgepixelx1=majoredgepixelx1;
                    edgepixely2=majoredgepixely2;
                    edgepixelx2=majoredgepixelx2;
                    Linecolor=[.8,0,0];
                else
                    edgepixely1=minoredgepixely1;
                    edgepixelx1=minoredgepixelx1;
                    edgepixely2=minoredgepixely2;
                    edgepixelx2=minoredgepixelx2;
                    Linecolor=[0,0,.8];                    
                end
                [temprows,tempcolumns] = size(BW3);
                Linepicturered=zeros(temprows,tempcolumns,3);
                Linepicturenegative=ones(temprows,tempcolumns,3);
                slopeaxis=((edgepixely1-toppixel+1)-(nverpointo-toppixel+1))/((edgepixelx1-leftpixel+1)-(nhorpointo-leftpixel+1));
                interceptaxis=(nverpointo-toppixel+1)-slopeaxis*(nhorpointo-leftpixel+1);   
                % Create the image for the case where the line in the  
                % image is between +/- 45º
                columndecimal=(nhorpointo-leftpixel+1)-round(nhorpointo-leftpixel+1-.5);
                rowdecimal=(nverpointo-toppixel+1)-round(nverpointo-toppixel+1-.5);
                if slopeaxis < 1 && slopeaxis > -1
                    for columncounter = 1 : tempcolumns - 1
                        if round(slopeaxis.*(columncounter+columndecimal)+interceptaxis) > 0 && round(slopeaxis.*(columncounter+columndecimal)+interceptaxis) <= temprows
                            columnposition=columncounter+columndecimal;
                            Linepicturered(round(slopeaxis.*columnposition+interceptaxis),round(columnposition),:)=Linecolor;
                            Linepicturenegative(round(slopeaxis.*columnposition+interceptaxis),round(columnposition),:)=[0,0,0];
                        end
                    end
                    % Otherwise draw a more vertical line:
                else
                    for rowcounter = 1 : temprows - 1
                        if round(((rowcounter+rowdecimal)-interceptaxis)./slopeaxis) > 0 && round(((rowcounter+rowdecimal)-interceptaxis)./slopeaxis) <= tempcolumns
                            rowposition=rowcounter+rowdecimal;
                            Linepicturered(round(rowposition),round((rowposition-interceptaxis)./slopeaxis),:)=Linecolor;
                            Linepicturenegative(round(rowposition),round((rowposition-interceptaxis)./slopeaxis),:)=[0,0,0];              
                        end
                    end
                end
                RGBim1=RGBim1.*Linepicturenegative+Linepicturered;
                % Highlight edge points along axis
                RGBim1(round(edgepixely1-toppixel+1),round(edgepixelx1-leftpixel+1),:)=round(Linecolor./max(Linecolor));
                RGBim1(round(edgepixely2-toppixel+1),round(edgepixelx2-leftpixel+1),:)=round(Linecolor./max(Linecolor));
            end
        end
        % Highlight the approximate location of the centroid in white.
        RGBim1(round(nverpointo-toppixel+1),round(nhorpointo-leftpixel+1),:)=[1,1,1];       
        cmap=gray(256);
        % Finalframe will be the output image.  Begin with full frame 
        % including header data.
        Finalframe=ind2rgb(Tempframe1,cmap);    
        % Brightening parameter to make tracked area stand out.  If you 
        % need to brighten the tracked area, .25 is a good gammaval. 
        % 1 = no effect on tracked area brightness.
        gammaval=.5;
        RGBim2=imadjust(Tempframe1(toppixel+height-dheight:bottompixel+height-dheight,leftpixel:rightpixel),[0;1],[0;1],gammaval);
        RGBim2=ind2rgb(RGBim2,cmap);
        % Multiply the images to highlight the brightened tracked area with color.
        Tempframe2=immultiply(RGBim1,RGBim2);
        % Now reset the untracked area to background
        BW3cat=cat(3,BW3,BW3,BW3);
        if headertop ~= false
            Finalframe(toppixel+height-dheight:bottompixel+height-dheight,leftpixel:rightpixel,:)=Finalframe(toppixel+height-dheight:bottompixel+height-dheight,leftpixel:rightpixel,:).*~BW3cat+Tempframe2(1:bottompixel-toppixel+1,1:rightpixel-leftpixel+1,:).*BW3cat;
        else
            Finalframe(toppixel:bottompixel,leftpixel:rightpixel,:)=Finalframe(toppixel:bottompixel,leftpixel:rightpixel,:).*~BW3cat+Tempframe2(1:bottompixel-toppixel+1,1:rightpixel-leftpixel+1,:).*BW3cat;
        end

    % It is simpler for generating a binary movie.
    else
        Finalframe=zeros(height,width);
        % Paste the header info on
        if headertop ~= false
            Finalframe(1:height-dheight,1:width)=im2bw(Tempframe1(1:height-dheight,1:width),1/tempmaxval);  
        end
        Finalframe(toppixel + height-dheight : bottompixel + height-dheight,leftpixel : rightpixel)=BW3;
         % Highlight the approximate location of the centroid in black.
        if headertop ~= false
            Finalframe(round(nverpointo)+height-dheight,round(nhorpointo),:)=0;
        else
            Finalframe(round(nverpointo),round(nhorpointo),:)=0;
        end
    Finalframe=uint8(Finalframe);
    Finalframe=255*Finalframe;
    end
    
    % If only a partial band of the image is desired to save disk size, 
    % strip away undersired areas.
    if reducewindow == true
        if height > dheight
            if headertop ~= false
                Finalframe = cat(1, Finalframe(1:height-dheight,:,:),...
                Finalframe(max(height-dheight+1, round(verpointo-1.5*beadsize*yres)) : min(height,round(verpointo+1.5*beadsize*yres)), :,:));
            else
                Finalframe = cat(1, Finalframe(max(1, round(verpointo-1.5*beadsize*yres)):min(round(verpointo+1.5*beadsize*yres),dheight),:,:),...
                Finalframe(dheight+1:height, :,:));
            end
        else
            Finalframe = Finalframe(max(1, round(verpointo-1.5*beadsize*yres)) : min(height,round(verpointo+1.5*beadsize*yres)), :,:);
        end
    end

    if binaryout == false
        frame = im2frame(Finalframe);
    else
        % Need to add 1 since cmap goes on interval 1->256 rather than
        % 0-> 255 like uint8 data type.
        frame = im2frame(Finalframe+1,gray(256));        
    end
    aviobj = addframe(aviobj,frame);

    % Let the user know every 10 frames
    if mod(counter,10)==0
        fprintf('Currently processing frame %g of %g.  ', counter, nmovieframes)
        toc
        % Save the center tracker every 10 frames in case execution is terminated early. 
        save(centertrackerfile,'Centertracker')
    end

    %Check to make sure the object is still in frame, stop when it is not.
    if ASPECTSTATS.BoundingBox(1) <= .5 || ASPECTSTATS.BoundingBox(2) <= .5
        inframe = false;
    end
    
    counter=counter+1;
end
aviobj = close(aviobj);

% Now save the output of the centroid location
save(centertrackerfile,'Centertracker')
% This can be reopened with the load command
end

function [majoraxisangle,majoredgepixelx1,majoredgepixely1,majoredgepixelx2,majoredgepixely2,majoraxislength] = findorientation(Perimeterpixelsrow,Perimeterpixelscolumn,xcentroidpixel,ycentroidpixel,xres,yres)
    % This function tests the centroid against the perimeter pixels to find 
    % the one that is located the furthest.
    % Convert pixels to distances
    Xcoordinates=Perimeterpixelscolumn./xres;
    Ycoordinates=Perimeterpixelsrow./yres;
    xcentroid=xcentroidpixel/xres;
    ycentroid=ycentroidpixel/yres;
    % Calculate the angles each pixel makes with the centroid 
    Angles=180/pi*atan2((Ycoordinates-ycentroid),(Xcoordinates-xcentroid));
    % Initialize
    longestlength=0;
    % Test the pixel values

    [npixels,temp]=size(Perimeterpixelsrow);
    for counter = 1 : npixels
        angledifferencemin=180;        
        xcurrent=Xcoordinates(counter);
        ycurrent=Ycoordinates(counter);
        anglecurrent=Angles(counter);
        if anglecurrent >= 0
            oppositeangle = anglecurrent-180;
        elseif anglecurrent < 0
            oppositeangle = anglecurrent+180;
        end
        angledifferencemin=180;
        for anglecounter = 1 : npixels
            angledifference=abs(Angles(anglecounter)-oppositeangle);
            if angledifference < angledifferencemin
                oppositeindex=anglecounter;
                angledifferencemin=angledifference;
                xopposite=Xcoordinates(anglecounter);
                yopposite=Ycoordinates(anglecounter);
            end
        end         
        length=((xcurrent-xopposite)^2+(ycurrent-yopposite)^2)^.5;
        if length > longestlength
            longestlength=length;
            bestindex=counter;
            bestoppositeindex=oppositeindex;
        end        
    end
    % Calculate the angle 
     majoraxisangle=Angles(bestindex);    
    % Restrict this angle to between +/- 90º
    if majoraxisangle > 90
        majoraxisangle=majoraxisangle-180;
    elseif majoraxisangle<90
        majoraxisangle=majoraxisangle+180;
    end
    % Record the values for reporting back to the program
    majoredgepixelx1=Perimeterpixelscolumn(bestindex);
    majoredgepixely1=Perimeterpixelsrow(bestindex);
    majoredgepixelx2=Perimeterpixelscolumn(bestoppositeindex);
    majoredgepixely2=Perimeterpixelsrow(bestoppositeindex);
    majoraxislength=(((majoredgepixelx1-majoredgepixelx2)/xres)^2+((majoredgepixely1-majoredgepixely2)/yres)^2)^.5;
end
    
function [edgepixelx1,edgepixely1,edgepixelx2,edgepixely2,axislength] = findedgepointsonline(Perimeterpixelsrow,Perimeterpixelscolumn,xcentroidpixel,ycentroidpixel,xres,yres,angle)
    % Restrict the angle of search to between -180º and +180º
    if angle > 180
        angle=angle-180;
    elseif angle <= -180
        angle=angle+180;
    end
    % Convert pixels to distances
    Xcoordinates=Perimeterpixelscolumn./xres;
    Ycoordinates=Perimeterpixelsrow./yres;
    xcentroid=xcentroidpixel/xres;
    ycentroid=ycentroidpixel/yres;
    [npixels,temp]=size(Perimeterpixelsrow);
    % Calculate the angles each pixel makes with the centroid 
    Angles=180/pi*atan2((Ycoordinates-ycentroid),(Xcoordinates-xcentroid));
    % Find the corresponding edge pixels along the major axis on the less
    % stretched side
    angledifferencemin=180;
    for anglecounter = 1 : npixels
        angledifference=abs(Angles(anglecounter)-angle);
        if angledifference < angledifferencemin
            bestcounter=anglecounter;
            angledifferencemin=angledifference;
        end
    end 
    % Now look on the other side for the nearest edge pixel.
    if angle > 0
        oppositeangle=angle-180;
    elseif angle <= 0
        oppositeangle=angle+180;
    end
    angledifferencemin=180;
    for anglecounter = 1 : npixels
        angledifference=abs(Angles(anglecounter)-oppositeangle);
        if angledifference < angledifferencemin
            oppositecounter=anglecounter;
            angledifferencemin=angledifference;
        end
    end     
    % Record the values for reporting back to the program
    edgepixelx1=Perimeterpixelscolumn(bestcounter);
    edgepixely1=Perimeterpixelsrow(bestcounter);
    edgepixelx2=Perimeterpixelscolumn(oppositecounter);
    edgepixely2=Perimeterpixelsrow(oppositecounter);
    axislength=(((edgepixelx1-edgepixelx2)/xres)^2+((edgepixely1-edgepixely2)/yres)^2)^.5;
end
