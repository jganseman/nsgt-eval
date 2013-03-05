function fr = nsigtf_real(c,gd,shift,Ls)
% NSIGTF_REAL  Nonstationary Gabor filterbank synthesis for real signals
%   Usage:  fr = nsigtf_real(c,gd,shift,M,Ls)
%
%   Input parameters: 
%         c         : Cell array of non-stationary Gabor coefficients
%         gd        : Cell array of synthesis filters
%         shift     : Vector of time shifts
%         M         : Number of time channels (vector/scalar)
%         Ls        : Length of the analyzed signal
%   Output parameters:
%         fr        : Synthesized real-valued signal (Channels are stored 
%                     in the columns)
%
%   Given the cell array *c* of non-stationary Gabor coefficients, and a 
%   set of filter, frequency shifts and time step parameters this function 
%   computes the corresponding inverse non-stationary Gabor transform.
%
%   This routine always assumes that the output is supposed to be
%   real-valued.
% 
%   If a non-stationary Gabor frame was used to produce the coefficients 
%   and *gd* is a corresponding dual frame, this function should give 
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

%Check input arguments
if nargin < 3
    error('Not enough input arguments');
end

if iscell(c) == 0 % If matrix format coefficients were used, convert to cell
   if ndims(c) == 2
       [N,chan_len] = size(c); CH = 1; 
       c = mat2cell(c.',chan_len,ones(1,N)).';
    else
       [N,chan_len,CH] = size(c);
       ctemp = mat2cell(permute(c,[2,1,3]),chan_len,ones(1,N),ones(1,CH));
       c = permute(ctemp,[2,3,1]); 
       clear ctemp;
   end
else
    [N,CH] = size(c);
end

timepos = cumsum(shift);        % Calculate positions from shift vector
NN = timepos(end);              % Length of the reconstruction before truncation
timepos = timepos-shift(1);   % Adjust positions

fr = zeros(NN,CH); % Initialize output

if nargin < 4
    Ls = NN; % If original signal length is not given do not truncate
end

% The overlap-add procedure including multiplication with the synthesis
% windows

for ii = 1:N
    Lg = length(gd{ii});
    
    win_range = mod(timepos(ii)+(-floor(Lg/2):ceil(Lg/2)-1),NN)+1;
    
    temp = fft(c{ii},[],1)*length(c{ii});
    temp = temp(mod([end-floor(Lg/2)+1:end,1:ceil(Lg/2)]-1,length(temp))+1,:);
    
    fr(win_range,:) = fr(win_range,:) + ...
        bsxfun(@times,temp,gd{ii}([Lg-floor(Lg/2)+1:Lg,1:ceil(Lg/2)]));
end

%fr(floor(Ls/2)+2:end,:) = conj(flipud(fr(2:ceil(Ls/2),:)));
%fr = real(ifft(fr));

fr = ifftreal(fr(1:floor(Ls/2)+1),Ls,1);

%fr = fr(1:Ls,:); % Truncate the signal to original length (if given)

end