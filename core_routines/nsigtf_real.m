function fr = nsigtf(c,gd,shift,Ls)

% NSIGTF.N - Gino Velasco, Nicki Holighaus 01.03.12
%
% fr = nsigtf(c,gd,shift,Ls)
%
% This is a modified version of nsigt.m for the case where the resolution 
% evolves over frequency.
%
% Given the cell array 'c' of non-stationary Gabor coefficients, and a set 
% of windows and frequency shifts, this function computes the corresponding 
% inverse non-stationary Gabor transform.
% 
% Input: 
%           c           : Cell array of non-stationary Gabor coefficients
%           gd          : Cell array of Fourier transforms of the synthesis 
%                         windows
%           shift       : Vector of frequency shifts
%           Ls          : Length of the analyzed signal
%
% Output:
%           fr          : Synthesized signal (Channels are stored in the
%                         columns)
% 
% If a non-stationary Gabor frame was used to produce the coefficients 
% and 'gd' is a corresponding dual frame, this function should give perfect 
% reconstruction of the analyzed signal (up to numerical errors).
% 
% The inverse transform is computed by simple 
% overlap-add. For each entry of the cell array c,
% the coefficients of frequencies around a certain 
% position in time, the Fourier transform
% is taken, giving 'frequency slices' of a signal.
% These slices are added onto each other with an overlap
% depending on the window lengths and positions, thus
% (re-)constructing the frequency side signal. For multichannel 
% signals, the overlap-add procedure is done for each channel. 
% In the last step, an inverse Fourier transform brings the 
% signal back to the time side.
% 
% More information can be found at:
% http://nuhag.eu/nonstatgab/
%
% This file was last updated for the Nonstationary Gabor Toolbox V0.02

% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to 
% Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

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