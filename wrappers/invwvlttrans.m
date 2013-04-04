function fr = invwvlttrans(c,g,shift,M,Ls,fb,tgtfl)
%INVWVLTTRANS  Wavelet frame synthesis
%   Usage:  fr = invwvlttrans(c,g,shift,M,Ls,fb,tgtfl)
%           fr = invwvlttrans(c,g,shift,M,Ls)
%           fr = invwvlttrans(c,g,shift,M,Ls)
%           fr = invwvlttrans(c,g,shift,M)
% 
%   Input parameters: 
%         c         : Cell array of Wavelet coefficients
%         g         : Cell array of Fourier transforms of the analysis 
%                     Wavelets
%         shift     : Vector of frequency shifts
%         M         : Number of time channels
%         fb	    : Frame bounds (vector)
%         tgtfl     : Tightflag (1 if frame is tight)
%   Output parameters:
%         fr        : Reconstructed signal
%
%   Given the cell array *c* and a painless Wavelet frame *g*, *shift*, 
%   *M*, this wrapper function performs the corresponding inverse Wavelet 
%   transform.
%
%   More information about the functions used can be found at:
%   http://univie.ac.at/nonstatgab/
%
%   See also:  wvlttrans, nsigtf, nsdual
%

% Author: Nicki Holighaus
% Date: 04.03.13

if nargin<6
    tgtfl = 0;
end

if nargin==6 && tgtfl==1
    fb=1;
    fprintf('No frame bounds provided, assuming them to be 1\n');
end

if nargin < 4
    error('Too few input arguments');
end

% Compute the dual frame
if tgtfl
    gd =g;
else
    gd = nsdual(g,shift,M);
end

% Inverse Wavelet transform
if nargin == 4
    fr = nsigtf(c,gd,shift);
else
    fr = nsigtf(c,gd,shift,Ls);
end

if tgtfl==1
    fr = fr/fb(1);
end