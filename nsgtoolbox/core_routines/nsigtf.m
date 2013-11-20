function fr = nsigtf(c,g,shift,Ls)
%NSIGTF  Nonstationary Gabor filterbank synthesis
%   Usage: fr = nsigtf(c,g,shift,Ls)
%          fr = nsigtf(c,g,shift)
%
%   Input parameters: 
%         c         : Cell array of nonstationary Gabor coefficients
%         g         : Cell array of synthesis filters
%         shift     : Vector of time shifts
%         Ls        : Length of the analyzed signal
%   Output parameters:
%         fr        : Synthesized real-valued signal (Channels are stored 
%                     in the columns)
%
%   Given the cell array *c* of nonstationary Gabor filterbank 
%   coefficients, a set of filters *g* and frequency shifts *shift*, this 
%   function computes the corresponding nonstationary Gabor filterbank
%   synthesis. Let `N=numel(g)` and $P(n)=\sum_{l=1}^{n} shift(l)$, then 
%   the synthesis formula reads:
%
%   ..               N-1 
%       fft(fr)(l) = sum sum c{n}(m)g{n}[l-P(n)]*exp(-2*pi*i*(l-P(n))*m/M(n)),
%                    n=0  m
%   
%   .. math::  \textbf{FFT}(fr)[l] = \sum_{n=0}^{N-1}\sum_{m} c\{n\}(m)g\{n\}[l-P(n)] e^{-2\pi i(l-P(n))m/M(n)},
%
%   for $l=0,\cdots,Ls-1$. The final reconstruction step then is 
%   `fr = ifft(fr)`. In practice, the synthesis formula is realized by `fft` 
%   and overlap-add.
% 
%   If a nonstationary Gabor frame was used to produce the coefficients 
%   and *g* is a corresponding dual frame, this function should perfectly 
%   reconstruct the originally analyzed signal to numerical precision.
%   
%   Multichannel output will save each channel in a column of *fr*.
%
%   See also:  nsgtf, nsdual, nstight
% 
%   References: badohojave11 dogrhove11 

% Author: Nicki Holighaus, Gino Velasco
% Date: 23.04.13

%Check input arguments
if nargin < 3
    error('Not enough input arguments');
end

if iscell(c) == 0 % If matrix format coefficients were used, convert to
    % cell
    [M,N,CH] = size(c);
    c = reshape(c,N*M,CH);
    c = mat2cell(c,M*ones(N,1),CH);
else
    N = length(c);
    CH = size(c{1},2);
end

posit = cumsum(shift);        % Calculate positions from shift vector
NN = posit(end);              % Reconstruction length before truncation
posit = posit-shift(1);     % Adjust positions

fr = zeros(NN,CH); % Initialize output

if nargin < 4
    Ls = NN; % If original signal length is not given do not truncate
end

% The overlap-add procedure including multiplication with the synthesis
% windows

for ii = 1:N
    Lg = length(g{ii});
    Lc = length(c{ii}(:,1));
    
    win_range = mod(posit(ii)+(-floor(Lg/2):ceil(Lg/2)-1),NN)+1;
    
    temp = fft(c{ii})*Lc;
    temp = temp(mod([end-floor(Lg/2)+1:end,1:ceil(Lg/2)]-1,Lc)+1,:);
    
    fr(win_range,:) = fr(win_range,:) + ...
        bsxfun(@times,temp,g{ii}([Lg-floor(Lg/2)+1:Lg,1:ceil(Lg/2)]));
end

fr = ifft(fr);
fr = fr(1:Ls,:); % Truncate the signal to original length (if given)

end