% invWVLTtrans.m - Nicki Holighaus 07.02.12
%
% f_rec = invWVLTtrans(c,g,shift,M,Ls,tgtfl)
%
% Given the cell array c and a painless Wavelet frame [g,shift,M], this
% wrapper function performs the corresponding inverse Wavelet transform
% 
% Input: 
%           c           : Cell array of Wavelet coefficients
%           g           : Cell array of Fourier transforms of the analysis 
%                         Wavelets
%           shift       : Vector of frequency shifts
%           M           : Number of time channels
%
% Output:
%           f_rec       : reconstructed signal
%
% More information about the functions used can be found at:
% http://nuhag.eu/nonstatgab/

% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to 
% Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

function fr = invWVLTtrans(c,g,shift,M,Ls,tgtfl,fb)

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