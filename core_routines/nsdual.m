function gd = nsdual(g,shift,M)
%NSDUAL  Canonical dual NSG frame (for painless systems)
%   Usage: gd = nsdual(g,shift,M)
%
%   Input parameters:
%         g         : Cell array of window functions
%         shift     : Vector of time/frequency shifts
%         M         : Number of frequency channels (vector/scalar)
%   Output parameters:
%         gd        : Dual window functions 
%
%   Computes (for the painless case) the canonical dual frame 
%   corresponding to a given non-stationary Gabor frame specified by the 
%   windows *g* and time shifts *shift*.
% 
%   Note, the time shifts corresponding to the dual window sequence is the
%   same as the original shift sequence and as such already given.
%
%   This routine's output can be used to achieve reconstruction of a signal 
%   from its non-stationary Gabor coefficients using the inverse 
%   non-stationary Gabor transform 'nsigt'.
% 
%   More information on Non-stationary Gabor transforms
%   can be found at:
%
%   http://univie.ac.at/nonstatgab/
%

% Author: Nicki Holighaus, Gino Velasco
% Date: 03.03.13

% Check input arguments

if nargin < 3
    for kk = 1:length(shift)
        M(kk) = length(g{kk}); M = M.';
    end
end

if nargin < 2
    error('Not enough input arguments');
end

if max(size(M)) == 1
    M = M(1)*ones(length(shift),1);
end

% Setup the necessary parameters
N = length(shift);

timepos = cumsum(shift);
Ls = timepos(N);
timepos = timepos-shift(1);

diagonal=zeros(Ls,1);
win_range = cell(N,1);

% Construct the diagonal of the frame operator matrix explicitly

for ii = 1:N
    Lg = length(g{ii});
    
    win_range{ii} = mod(timepos(ii)+(-floor(Lg/2):ceil(Lg/2)-1),Ls)+1;
    diagonal(win_range{ii}) = diagonal(win_range{ii}) + ...
        (fftshift(g{ii}).^2)*M(ii);
end

% Using the frame operator and the original window sequence, compute
% the dual window sequence

gd = g;

for ii=1:N
    gd{ii} = ifftshift(fftshift(gd{ii})./diagonal(win_range{ii}));
end