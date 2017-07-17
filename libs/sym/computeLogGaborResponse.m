function EO = computeLogGaborResponse(img, wavParam)
%%
% Log-Gabor Filters are constructed in terms of two components : 
% 1) The radial component, which controls the frequency band that the filter
%    responds to
% 2) The angular component, which controls the orientation that the filter
%    responds to.
% The two components are multiplied together to construct the overall filter.
%
% Returned value :
% EO - A 2D cell array of complex valued convolution results
%%
% Copyright (c) 1996-2010 Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% peter.kovesi at uwa edu au
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
%
% References:
%
%     Peter Kovesi, "Image Features From Phase Congruency". Videre: A
%     Journal of Computer Vision Research. MIT Press. Volume 1, Number 3,
%     Summer 1999 http://mitpress.mit.edu/e-journals/Videre/001/v13.html
%
%     Peter Kovesi, "Phase Congruency Detects Corners and
%     Edges". Proceedings DICTA 2003, Sydney Dec 10-12
%%
rows = size(img,1); %- Number of rows inside an image
cols = size(img,2); %- Number of columns inside an image
minWaveLength = wavParam.minWaveLength; %- Wavelength of smallest scale filter.
mult = wavParam.mult; %- Scaling factor between successive filters.
%- Ratio of the standard deviation of the Gaussian describing 
%  the log Gabor filter's transfer function in the frequency domain to 
%  the filter center frequency.
radSigma = wavParam.radSigma;
%- Ratio of angular interval between filter orientations
%  and the standard deviation of the angular Gaussian
%  function used to construct filters in the frequency plane.
angSigma = wavParam.angSigma;
nAngs = wavParam.nAngs; %- Number of wavelet scales.
nScls = wavParam.nScls; %- Number of filter orientations.
%%
% Set up X and Y matrices with ranges normalised to +/- 0.5
% The following code adjusts things appropriately for odd and even values
% of rows and columns.

if mod(cols,2)
    xrange = [-(cols-1)/2:(cols-1)/2]/(cols-1);
else
    xrange = [-cols/2:(cols/2-1)]/cols; 
end

if mod(rows,2)
    yrange = [-(rows-1)/2:(rows-1)/2]/(rows-1);
else
    yrange = [-rows/2:(rows/2-1)]/rows; 
end
[x,y] = meshgrid(xrange, yrange);

radius = sqrt(x.^2 + y.^2); %- Matrix values contain *normalised* radius from centre.     
theta = atan2(-y,x); %- Matrix values contain polar angle.

%- Quadrant shift radius and theta so that filters are constructed 
%  with 0 frequency at the corners.
radius = ifftshift(radius);      
theta  = ifftshift(theta);

%- Get rid of the 0 radius value at the 0 frequency point (now at top-left corner)
% so that taking the log of the radius will not cause trouble.
radius(1,1) = 1;                 
                                 
sintheta = sin(theta);
costheta = cos(theta);
clear x; clear y; clear theta;
%%
%- Construct the radial filter component :
%  First construct a low-pass filter that is as large as possible, falls
%  away to zero at the boundaries.  All log Gabor filters are multiplied by
%  this to ensure no extra frequencies at the 'corners' of the FFT are
%  incorporated as this seems to upset the normalisation process when
%  calculating phase congrunecy.

lp = lowpassfilter([rows,cols],.45,15);
logGabor = cell(1,nScls);
for s = 1:nScls
    wavelength = minWaveLength*mult^(s-1);
    fo = 1.0/wavelength; %- Centre frequency of filter.
    logGabor{s} = exp((-(log(radius/fo)).^2) / (2 * log(radSigma)^2));
    logGabor{s} = logGabor{s}.*lp; %- Apply low-pass filter
    logGabor{s}(1,1) = 0; %- Set the value at the 0 frequency point of the filter
                          %  back to zero (undo the radius fudge).
end
%%
%- Construct the angular filter component :
%  For each point in the filter matrix calculate the angular distance from
%  the specified filter orientation.  To overcome the angular wrap-around
%  problem sine difference and cosine difference values are first computed
%  and then the atan2 function is used to determine angular distance.

thetaSigma = pi/nAngs/angSigma;
spread = cell(1,nAngs);
for o = 1:nAngs
  angl = (o-1)*pi/nAngs; %- Filter angle.
  ds = sintheta * cos(angl) - costheta * sin(angl); %- Difference in sine.
  dc = costheta * cos(angl) + sintheta * sin(angl); %- Difference in cosine.
  dtheta = abs(atan2(ds,dc)); %- Absolute angular distance.
  spread{o} = exp((-dtheta.^2) / (2*thetaSigma^2)); %- Calculate the angular filter component.
end
%%
%- Compute complex valued convolution results respect to input image

imagefft = fft2(img); %- Fourier transform of image
EO = cell(nScls, nAngs);
for o = 1:nAngs %- For each orientation.
    for s = 1:nScls %- For each scale.
        %- Multiply radial and angular components to get the filter.
        filter = logGabor{s} .* spread{o}; 
        %- Convolve image with even and odd filters returning the result in
        %  EO.
        EO{s,o} = ifft2(imagefft .* filter);
    end
end
end

