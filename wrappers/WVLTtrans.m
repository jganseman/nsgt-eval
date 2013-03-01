function [c,g,shift,M,Ls,fb,tgtfl] = WVLTtrans(f,fmin,sr,bins,bw,winfun)

% WVLTtrans.m - Nicki Holighaus, Christoph Wiesmeyr 15.05.12
%
% [c,g,shift,M,Ls] = WVLTtrans(f,fmin,sr,bins,bw)
% [c,g,shift,M,Ls] = WVLTtrans(f,fmin,sr,bins)
% [c,g,shift,M,Ls] = WVLTtrans(f,fmin,sr)
% [c,g,shift,M,Ls] = WVLTtrans(f,fmin)
%
% Given the function f and the necessary parameters described below, this
% wrapper function performs the corresponding Wavelet transform
%
% Input: 
%           f           : Function to analyze
%           fmin        : Desired minimum center frequency (in Hz)
%           bandwidth   : Desired bandwidth in the first frequency band 
%                         (in Hz)
%           bins        : Desired number of bins per octave
%           sr          : Sampling rate of f (in Hz)
%           winfun      : window function to be used
%
% Output: 
%           c           : Cell array of Wavelet coefficients
%           g           : Cell array of Fourier transforms of the analysis 
%                         Wavelets
%           shift       : Vector of frequency shifts
%           M           : Number of time channels
%           Ls          : Original signal length
%	    fb		: frame bounds (vector)
% 	    tgtfl 	: tightflag (1 if frame is tight)
%
%
% More information about the functions used can be found at:
% http://nuhag.eu/nonstatgab/

% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to 
% Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.


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