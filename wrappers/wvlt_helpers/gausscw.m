function g = gausscw(x)
%GAUSSCW  Gaussian window
%   Usage:  g = gausscw(x)
%
%   Input parameters: 
%         x         : Vector of sampling positions
%   Output parameters:
%         g         : Output window
%
%   Minimal implementation of a Gaussian window for usage in the Wavelet
%   wrapper.

% Author: Christoph Wiesmeyr
% Date: 04.03.13

g = exp(-18*x.^2);
g = g.';