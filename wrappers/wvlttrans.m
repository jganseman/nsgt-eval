function [c,g,shift,M,Ls,fb,tgtfl] = wvlttrans(f,fmin,sr,bins,bw,winfun)
%WVLTTRANS  Wavelet frame transform
%   Usage: [c,g,shift,M,Ls,fb,tgtfl] = wvlttrans(f,fmin,sr,bins,bw,tgtfl)
%          [c,g,shift,M,Ls,fb,tgtfl] = wvlttrans(f,fmin,sr,bins,bw)
%          [c,g,shift,M,Ls,fb,tgtfl] = wvlttrans(f,fmin,sr,bins)
%          [c,g,shift,M,Ls,fb,tgtfl] = wvlttrans(f,fmin,sr)
%          [c,g,shift,M,Ls,fb,tgtfl] = wvlttrans(f,fmin)
%          [c,g,shift,M,Ls,fb] = wvlttrans(...)
%          [c,g,shift,M,Ls] = wvlttrans(...)
%          [c,g,shift,M] = wvlttrans(...)
%          c = wvlttrans(...)
%
%   Input parameters: 
%         f         : Input signal
%         fmin      : Desired minimum center frequency (in Hz)
%         bandwidth : Desired bandwidth in the first frequency band (in Hz)
%         bins      : Desired number of bins per octave
%         sr        : Sampling rate of f (in Hz)
%         winfun    : Window function to be used
%   Output parameters: 
%         c         : Cell array of Wavelet coefficients
%         g         : Cell array of Fourier transforms of the analysis
%                     Wavelets
%         shift     : Vector of frequency shifts
%         M         : Number of time channels
%         Ls        : Original signal length
%         fb	    : Frame bounds (vector)
%         tgtfl     : Tightflag (1 if frame is tight)
%
%   Given the function *f* and the necessary parameters, this wrapper 
%   function performs the corresponding Wavelet transform
%
%   More information about the functions used can be found at:
%   http://univie.ac.at/nonstatgab/
%

% Author: Nicki Holighaus, Christoph Wiesmeyr
% Date: 04.03.13


% set default window function
if nargin < 6
    winfun = @hannwin;
end

% Determine signal length and make column vectors
[Ls,w]=size(f);
if Ls<w
    f=f.';
    Ls=w;
end

if nargin == 4
    fac = 2^(2/bins)-2^(-2/bins);
    bw = fmin*fac;
elseif nargin == 3
    bins = 4;
    fac = 2^(2/bins)-2^(-2/bins);
    bw = fmin*fac;
elseif nargin == 2
    bins = 4;
    fac = 2^(2/bins)-2^(-2/bins);
    bw = fmin*fac;
    sr = Ls;
elseif nargin < 2
    error('Invalid number of input arguments');
end



% Compute Wavelet system
if nargout <= 5
    [g,shift,M] = nsgwvltwin(fmin,bw,bins,sr,Ls,winfun);
else
    [g,shift,M,fb] = nsgwvltwin(fmin,bw,bins,sr,Ls,winfun);
    tgtfl = (abs(fb(1)-fb(2))<10^(-12));
end

% Compute Wavelet coefficients
c = nsgtf(f,g,shift,M);