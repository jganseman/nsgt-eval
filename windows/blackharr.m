function g = blackharr(x,L,mod)
%BLACKHARR  Blackman-Harris window 
%   Usage:  g = blackharr(x,mod)
%           g = blackharr(x)
%           g = blackharr(N,L,mod)
%           g = blackharr(N,L)
%           g = blackharr(N)
%
%   Input parameters: 
%         x         : Vector of sampling positions
%         N         : Window support (in samples)
%         L         : Output length (in samples)
%         mod       : Use a modified continuous version of the 
%                     Blackman-Harris window
%   Output parameters:
%         g         : Output window
%
%   Given a vector *x*, computes the values of the Blackman-Harris window
%   at the sampling points given. Given scalars *N* and *L*, computes a 
%   Blackman-Harris window of length *N*, full-point centered on a vector
%   of length *L*.
%
%   By default, a slightly modified version of the original Blackman-Harris
%   construction is used to yield a `continuous` window function. The
%   original construction can be used by setting the parameter *mod* to 0.
%
%   See also:  hannwin
%

% Author: Nicki Holighaus
% Date: 04.03.13

if nargin < 1
    error('Not enough input arguments');
end

if numel(x) == 1
    N = x;
    if nargin < 3
        mod = 1;
        if nargin < 2
            L = N;
        end
    end
    if L<N
        error('Output length L must be larger than or equal to N');
    end
    x = [0:ceil(N/2)-1,-N*ones(1,L-N),-floor(N/2):-1]'/N;
elseif nargin < 2
    mod = 1;
end

if size(x,2) > 1
    x = x.';
end

x = x+1/2;
if mod == 0
    g = 0.35875 - 0.48829*cos(2*pi*x) + 0.14128*cos(4*pi*x) - ...
        0.01168*cos(6*pi*x);
else
    g = 0.35872 - 0.48832*cos(2*pi*x) + 0.14128*cos(4*pi*x) - ...
        0.01168*cos(6*pi*x);
end
g = g .* (x > 0) .* (x < 1);