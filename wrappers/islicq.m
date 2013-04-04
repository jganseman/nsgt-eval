function [fr,gd] = islicq(c,g,shift,M,Ls,sl_len,tr_area)
%ISLICQ  Sliced constant-Q/variable-Q synthesis
%   Usage:  [fr,gd] = islicq(c,g,shift,M,Ls,sl_len,tr_area)
%           fr = islicq(c,g,shift,M,Ls,sl_len,tr_area)
%
%   Input parameters:
%         c         : Cell array of coefficients
%         g         : Cell array of analysis windows
%         shift     : Vector of frequency shifts of windows
%         M         : Number of channels (vector/scalar)
%         Ls        : Original signal length
%         sl_len    : Slice length
%         tr_area   : Transition area length
%   Output parameters:
%         fr        : Reconstructed signal
%         gd        : Cell array of synthesis windows
%
%   This wrapper function implements the inverse constant-Q nonstationary
%   Gabor transform. Given the output of |slicq|, this routine 
%   synthesizes the signal *fr* from the coefficient array *c*. If *c* is 
%   unaltered and the system specified by *g*, *shift* and *M* is a 
%   painless frame, then *fr* will be equal to the signal analyzed by 
%   |slicq|.
%
%   For more information on the constant-Q nonstationary Gabor transform, 
%   see
%   http://www.univie.ac.at/nonstatgab/
%
%   For more information on sliced transforms, see
%   N. Holighaus, M. Doerfler, G. Velasco, and T. Grill, "A framework for 
%   invertible, real-time constant-q transforms," Audio, Speech, and 
%   Language Processing, IEEE Transactions on, vol. 21, no. 4, 
%   pp. 775-785, April 2013.
%
%   See also:  slicq, nsigtf, nsdual, unslicing
%

% Author: Nicki Holighaus
% Date: 04.03.13

if nargin < 7
    error('Not enough input arguments');
end

% Before reconstruction, the coefficients have to be rearranged

if iscell(c) == 0 % Matrix coefficients
    slices = size(c,3);
    A = [3*M(1)/4+1:M(1),1:3*M(1)/4];
    B = [M(1)/4+1:M(1),1:M(1)/4];
    
    c(:,:,1:2:slices) = c(B,:,1:2:slices);
    c(:,:,2:2:slices) = c(A,:,2:2:slices);
else % Cell array coefficients
    slices = size(c{1},2);
    for jj = 1:length(g) % Can this be faster somehow?
        A = [3*M(jj)/4+1:M(jj),1:3*M(jj)/4];
        B = [M(jj)/4+1:M(jj),1:M(jj)/4];
        
        c{jj}(:,1:2:end) = c{jj}(B,1:2:end);
        c{jj}(:,2:2:end) = c{jj}(A,2:2:end);
    end
end

gd = nsdual(g,shift,M); % Compute the dual frame

% The actual transform
fr = isliCQ_int(c,gd,shift,Ls,sl_len,tr_area,slices);

end

%% Internal computation of the transform

function fr = isliCQ_int(c,gd,shift,Ls,sl_len,tr_area,slices)

%% Compute the inverse CQ of each slice
frec_sliced = nsigtf(c,gd,shift,sl_len);
%% Glue the parts back together
fr = unslicing(real(frec_sliced),sl_len,tr_area,slices);
fr = fr(1:Ls);

end