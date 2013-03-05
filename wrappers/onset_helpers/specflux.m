function [SF,V0] = specflux(f,win_length,tgap)
%SPECFLUX  Spectral flux onset detection function
%   Usage: [SF,V0] = specflux(f,win_length,tgap)
%          SF = specflux(f,win_length,tgap)
%
%   Input parameters: 
%         f         : Input signal
%         win_length: Desired window length for the STFT
%         tgap      : Time step for the STFT
%   Output parameters:
%         SF        : spectral flux of *f*
%         V0        : STFT coefficients of *f*
% 
%   This is a helper function for `onsetdet` and not meant to
%   be used individually.
%
%   Uses routines from LTFAT 0.97 or higher, available at:
%   http://ltfat.sourceforge.net/
%
%   Computes the spectral flux onset-detection function
%   of *f* with a Hann window of length *win_length*. 
%   The STFT is taken with time shift parameter *tgap*
%   and *win_length* frequency channels.
%
%   External: DGT (LTFAT routine)

% Author: Nicki Holighaus
% Date: 04.03.13

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