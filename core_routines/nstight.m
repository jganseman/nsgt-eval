function gt = nstight(g,shift,M)

% NSTIGHT.M - Nicki Holighaus 02.02.11
%
% gt = nstight(g,shift,M)
%
% Computes (for the painless case) the tight frame corresponding to a given 
% non-stationary Gabor frame specified by the windows 'g' and time shifts
% 'shift'.
% 
% Note, the time shifts corresponding to the tight window sequence is the
% same as the original shift sequence and as such already given.
%
% This routine's output can be used for decomposition and reconstruction of 
% a signal into, resp. from its non-stationary Gabor coefficients using the 
% 'nsgt' and 'nsigt'.
% 
% More information on Non-stationary Gabor transforms
% can be found at:
%
% http://univie.ac.at/nonstatgab/
%
% minor edit by Gino Velasco 23.02.11

% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to 
% Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

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
  diagonal(win_range{ii}) = diagonal(win_range{ii}) + (fftshift(g{ii}).^2)*M(ii);   
end

% Using the frame operator and the original window sequence, compute 
% the dual window sequence

gt = g;
  
for ii=1:length(timepos)
    gt{ii} = ifftshift(fftshift(gt{ii})./sqrt(diagonal(win_range{ii})));
end
