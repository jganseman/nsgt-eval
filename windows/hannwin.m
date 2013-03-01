function g = hannwin(x,L)

% HANNWIN.M - Nicki Holighaus 09.03.11
%
% g = hannwin(L)
%
% Computes a Hann window, centered around 
% sample 1, so that the FFT is real-valued.

% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to 
% Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

if nargin < 1
    error('Not enough input arguments');
end

if numel(x) == 1
    N = x;
    if nargin < 2
        L = N;
    end
    x = [0:ceil(N/2)-1,-N*ones(1,L-N),-floor(N/2):-1]'/N;
end

if size(x,2) > 1
    x = x.';
end

g = .5 + .5*cos(2*pi*x);
g = g .* (x > -1/2) .* (x < 1/2);