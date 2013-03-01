function [f_rec,gd] = isliCQ(c,g,shift,M,Ls,sl_len,tr_area)

% ISLICQ.M - Nicki Holighaus 02.03.12
%
% This wrapper function implements the inverse constant-Q nonstationary
% Gabor transform. Given the output of 'sliCQ.m', this routine synthesizes
% the signal 'f_rec' from the coefficient array 'c'. If 'c' is unaltered
% and the system specified by 'g','shift' and 'M' is a painless frame, then
% 'f_rec' will be equal to the signal analyzed by 'sliCQ.m'.
%
% For more information on the constant-Q nonstationary Gabor transform, see
% http://www.univie.ac.at/nonstatgab/
% For more information on sliced transforms, see
% '(dohogrveXX - The sliCQ paper reference goes here)'
%
% Usage:
%       [f_rec,gd] = isliCQ(c,g,shift,M,Ls,sl_len,tr_area)
%
% Description text
%
%   Input parameters:
%         c        : Cell array of coefficients.
%         g        : Cell array of analysis windows.
%         shift    : Vector of frequency shifts of windows.
%         M        : Number of channels (vector or constant).
%         Ls       : Original signal length.
%         sl_len   : slice length.
%         tr_area  : transition area length.
%
%   Output parameters:
%         f_rec    : Reconstructed signal.
%         gd       : Cell array of synthesis windows.

% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to 
% Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
    
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
    f_rec = isliCQ_int(c,gd,shift,Ls,sl_len,tr_area,slices); 
        
end

%% Internal computation of the transform

function f_rec = isliCQ_int(c,gd,shift,Ls,sl_len,tr_area,slices)

    %% Compute the inverse CQ of each slice
    frec_sliced = nsigtf(c,gd,shift,sl_len);
    %% Glue the parts back together
    f_rec = unslicing(real(frec_sliced),sl_len,tr_area,slices);
    f_rec = f_rec(1:Ls);
      
end