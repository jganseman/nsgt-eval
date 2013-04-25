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
%         sr        : Sampling rate of f (in Hz)
%         bins      : Desired number of bins per octave
%         bw        : Desired bandwidth in the first frequency band (in Hz)
%         winfun    : String containing the desired window function name
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
%   This is a wrapper function for the painless Wavelet transform via
%   nonstationary Gabor filterbank. Given a signal *f* and minimum
%   frequency *fmin*, a tight system with $4$ scales per octave is
%   constructed using logarithmically sampled Hann windows with 3/4
%   overlap. The additional parameters *sr*, *bins*, *bw* and *winfun* can
%   be specified to individually construct different Wavelet systems.
%
%   To construct systems with specific overlap factors $(n-1)/n$, choose
%   $bw = 2^{n/(2 bins)}-2^{-n/(2 bins)}$.
%
%   In addition to the Wavelet coefficients *c*, also the analysis system
%   *g*, *shift*, *M* can be returned, as can the length *Ls* of the input
%   signal *f*, the frame bounds of the system *g*, *shift*, *M* and a flag
%   indicating if a tight frame was used. These parameters are necessary to
%   perform reconstruction with the inverse Wavelet transform wrapper
%   |invwvlttrans|.
%
%   See also:  invwvlttrans, nsgtf, nsgwvltwin
%
%   References: badohojave11

% Author: Nicki Holighaus, Christoph Wiesmeyr
% Date: 25.04.13


% set default window function
if nargin < 6
    winfun = 'hann';
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