function fr = nsigt(c,gd,shift,Ls)
% NSIGT  Nonstationary Gabor synthesis
%   Usage: fr = nsigt(c,gd,shift,Ls)
%
%   Input parameters: 
%         c         : Cell array of non-stationary Gabor coefficients
%         gd        : Cell array of synthesis windows
%         shift     : Vector of time shifts
%         Ls        : Length of the analyzed signal
%
%   Output parameters:
%         fr        : Synthesized signal (Channels are stored in the
%                     columns)
%
%   Given the cell array *c* of non-stationary Gabor coefficients, and a 
%   set of windows and time shifts, this function computes the 
%   corresponding inverse non-stationary Gabor transform.
% 
%   If a non-stationary Gabor frame was used to produce the coefficients 
%   and 'gd' is a corresponding dual frame, this function should give 
%   perfect reconstruction of the analyzed signal (up to numerical errors).
% 
%   The inverse transform is computed by simple 
%   overlap-add. For each entry of the cell array *c*,
%   the coefficients of frequencies around a certain 
%   position in time, the inverse Fourier transform
%   is taken, giving 'time slices' of a signal.
%   These slices are added onto each other with an overlap
%   depending on the window lengths and positions, thus
%   (re-)constructing a signal. For multichannel signals, the overlap-add
%   procedure is done for each channel.
% 
%   More information can be found at:
%   http://univie.ac.at/nonstatgab/
%

% Author: Nicki Holighaus, Gino Velasco
% Date: 03.03.13

% some preparation

if nargin < 3
    error('Not enough input arguments');
end

if iscell(c) == 0 % If matrix format coefficients were used, convert to cell
    [M,N,CH] = size(c);
    c = reshape(c,N*M,CH);
    c = mat2cell(c,M*ones(N,1),CH);
else
    N = length(c);
    CH = size(c{1},2);
end

timepos = cumsum(shift);        % Calculate positions from shift vector
NN = timepos(end);              % Length of the reconstruction before truncation
timepos = timepos-shift(1);     % Adjust positions

fr = zeros(NN,CH); % Initialize output

if nargin < 4
    Ls = NN; % If original signal length is not given do not truncate
end

% The overlap-add procedure including multiplication with the synthesis
% windows

for ii = 1:N
    Lg = length(gd{ii});
    
    win_range = mod(timepos(ii)+(-floor(Lg/2):ceil(Lg/2)-1),NN)+1;
    
    temp = ifft(c{ii},[],1)*length(c{ii});
    temp = temp(mod([end-floor(Lg/2)+1:end,1:ceil(Lg/2)]-1,length(temp))+1,:);
    
    fr(win_range,:) = fr(win_range,:) + ...
        bsxfun(@times,temp,gd{ii}([Lg-floor(Lg/2)+1:Lg,1:ceil(Lg/2)]));
end

fr = fr(1:Ls,:); % Truncate the signal to original length (if given)

end