function g = hannwin(x,L)
%HANNWIN  Hann window 
%   Usage:  g = hannwin(x)
%           g = hannwin(N,L)
%           g = hannwin(N)
%
%   Input parameters: 
%         x         : Vector of sampling positions
%         N         : Window support (in samples)
%         L         : Output length (in samples)
%   Output parameters:
%         g         : Output window
%
%   Given a vector *x*, computes the values of the Hann window at the 
%   sampling points given. Given scalars *N* and *L*, computes a Hann 
%   window of length *N*, full-point centered on a vector of length *L*.
%

% Author: Nicki Holighaus
% Date: 04.03.13

if nargin < 1
    error('Not enough input arguments');
end

if numel(x) == 1
    N = x;
    if nargin < 2
        L = N;
    end
    if L<N
        error('Output length L must be larger than or equal to N');
    end
    x = [0:ceil(N/2)-1,-N*ones(1,L-N),-floor(N/2):-1]'/N;
end

if size(x,2) > 1
    x = x.';
end

g = .5 + .5*cos(2*pi*x);
g = g .* (x > -1/2) .* (x < 1/2);