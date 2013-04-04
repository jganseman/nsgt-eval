function g = wp2inp(x)
%WP2INP  Wavelet uncertainty minimizer (UnlocX workpage 2)
%   Usage:  g = wp2inp(x)
%
%   Input parameters: 
%         x         : Vector of sampling positions
%   Output parameters:
%         g         : Output window
%
%   Minimal implementation of the Wavelet uncertainty minimzier for usage
%   in the Wavelet wrapper.
%
%   Insert a reference here!
%
%   See also:  wvlttrans, invwvlttrans, nsgwvltwin
%

% Author: Christoph Wiesmeyr
% Date: 04.03.13

x = x.';
g = exp(exp(-2*x)*25.*(1+2*x)).*(abs(x)<=1/2);
g = g/max(g);