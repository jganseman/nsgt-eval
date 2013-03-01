% SPECFLUX.M - Nicki Holighaus 02.02.11
%
% This is a helper function for 'onsetdet' and not meant to
% be used individually.
%
% Uses routines from LTFAT 0.97 or higher, available at:
% http://ltfat.sourceforge.net/
%
% Computes the spectral flux onset-detection function
% of f with a Hann window of length win_length. 
% The STFT is taken with time shift parameter tgap
% and win_length frequency channels.
%
% External: DGT (LTFAT routine)

% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to 
% Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

function [SF,tgap,V0] = specflux(f,win_length,tgap)

% Check input arguments

if nargin < 3
    error('Not enough input arguments');
end

% Compute the Gabor transform (sampled STFT) of f

win=hannwin(win_length);

V0 = dgt(f,win,tgap,win_length);

% Compute the spectral flux 

VV = abs(V0);
VV = max(VV-circshift(VV,[0,1]),0);

SF = sum(VV);

% Normalize

SF = SF-mean(SF);
SF = SF./std(SF);