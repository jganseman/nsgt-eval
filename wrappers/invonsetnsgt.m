function fr = invonsetnsgt(c,g,shift,M,Ls)
%INVONSETNSGT  Onset-based nonstationary Gabor synthesis
%   Usage:  fr = invonsetnsgt(c,g,shift,M,Ls)
%           fr = invonsetnsgt(c,g,shift,M)  
%
%   Input parameters: 
%         c         : Cell array of transform coefficients
%         g         : Cell array of analysis windows
%         shift     : Vector of time shifts
%         M         : Number of time channels
%         Ls        : Original signal length
%   Output parameters: 
%         fr        : Reconstructed signal
%
%   This is a wrapper function for the inverse scaleframe nonstationary 
%   Gabor transform with onset detection based adaptation. It basically
%   just forwards the input to |nsdual| and |nsigt_real|. 
%
%   For more information see |onsetnsgt|.
%
%   See also:  onsetnsgt, nsigt_real, nsdual
%
%   References:  badohojave11

% Author: Nicki Holighaus
% Date: 25.04.13

if nargin < 5
    if nargin < 4
        error('Not enough input arguments');
    end
    Ls = sum(shift);
end

gd = nsdual(g,shift,M);

fr = nsigt_real(c,gd,shift,M,Ls);