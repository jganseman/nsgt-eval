function [c,g,shift,M,Ls,sl_len,tr_area] = sliCQ(f,fmin,fmax,bins,sl_len,tr_area,sr,M,min_win,Qvar)

% SLICQ.M - Nicki Holighaus 02.03.12
%
% This wrapper function computes the sliced constant-Q nonstationary Gabor
% transform of the signal 'f' given by the specified parameter set. The 
% parameters are identical to the parameters for 'nsgcqwin.m' plus some
% additional ones. Those additional parameters are the slice and
% transition area lengths as well as (optionally) the desired number of frequency
% channels 'M', if a fixed number of channels is desired, and the
% parameter 'Qvar' that allows variation of the overlap between frequency
% channels.
%
% The output parameters can be used to reconstruct the signal using the
% inverse sliced constant-Q 'isliCQ.m', provided the system specified by 
% 'g','shift' and 'M' is a painless frame.
%
% For more information on the constant-Q nonstationary Gabor transform, see
% http://www.univie.ac.at/nonstatgab/
% For more information on sliced transforms, see
% '(dohogrveXX - The sliCQ paper reference goes here)'
%
% Usage:
%       [c,g,shift,M,Ls,sl_len,tr_area] = sliCQ(f,fmin,fmax,bins,sl_len,tr_area,sr,M,min_win,Qvar)
%
% Description text
%
%   Input parameters:
%         f        : input signal.
%         fmin     : desired minimum frequency (in Hz).
%         fmax     : desired maximum frequency (in Hz).
%         bins     : bins per octave (constant or vector (for VQ)).
%         sl_len   : desired slice length (in samples).
%         tr_area  : transition area length (in samples, <= sl_len/2).
%         sr       : sampling rate (in Hz).
%         M        : desired number of time channels per slice, if set to 0,
%                    a channel vactor will be computed ('M' must be a 
%                    multiple of 4 or will be set to '4*ceil(M/4)').
%         min_win  : minimum window bandwidth (default 16 samples).
%         Qvar     : Factor varying the bandwidth. 'Qvar=X' leads to a Q-factor 
%                    of 'Q/X'.
%
%   Output parameters:
%         c        : Cell array of coefficients.
%         g        : Cell array of analysis windows.
%         shift    : Vector of frequency shifts of windows.
%         M        : Number of channels (vector or constant).
%         Ls       : Original signal length.
%         sl_len   : slice length.
%         tr_area  : transition area length.

% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to 
% Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

    if size(f,1) == 1
        f = f.';
    end

    Ls = length(f);
    
    if nargin < 10 
        Qvar = 1;
        if nargin < 9
            min_win = 16;
            if nargin < 8
                M = 0;
                if nargin < 7
                    sr = 1;
                    if nargin < 6
                        if nargin < 5
                            sl_len = 16384;
                            if nargin < 4
                                bins = 12;
                                if nargin < 3
                                    fmax = .5;
                                    if nargin < 2
                                        error('Too few input arguments');
                                    end
                                end
                            end
                        end
                        tr_area = round(sl_len/16);
                    end
                end
            end
        end
    end
    
    if numel(M) > 1
        warning('Number of channels must be a constant or 0, channel vector will be recomputed');
        M(1) = 0;
    end
    
    if M(1) == 0 
        % This is just a slightly modified version of nsgfwin
        [g,shift,M] = nsgcqwin(fmin,fmax,bins,sr,sl_len,'min_win',min_win,'Qvar',Qvar,...
            'wL_fac',4,'fractional',1,'winfun',@blackharr);
    else
        % This is just a slightly modified version of nsgfwin
        [g,shift] = nsgcqwin(fmin,fmax,bins,sr,sl_len,'min_win',min_win,'Qvar',Qvar,...
            'wL_fac',4,'fractional',1,'winfun',@blackharr);
        M = 4*ceil(M/4);
    end
    N = length(shift); % The number of filters
    
    %% Compute the slices
    f_sliced = slicing(f,sl_len,tr_area,Ls);
    
    %% Compute the CQ of each slice
    c = nsgtf(f_sliced,g,shift,M);

    %% Rearrange the coefficients such that they the slices are centered
     
    if iscell(c) == 0 % Matrix coefficients
        A = [3*M(1)/4+1:M(1),1:3*M(1)/4];
        B = [M(1)/4+1:M(1),1:M(1)/4];
        
        c(:,:,1:2:end) = c(A,:,1:2:end);
        c(:,:,2:2:end) = c(B,:,2:2:end);
    else % Cell array coefficients
         for jj = 1:length(g) % Can this be faster somehow?            
            A = [3*M(jj)/4+1:M(jj),1:3*M(jj)/4];
            B = [M(jj)/4+1:M(jj),1:M(jj)/4];
  
            c{jj}(:,1:2:end) = c{jj}(A,1:2:end);
            c{jj}(:,2:2:end) = c{jj}(B,2:2:end);
        end
     end