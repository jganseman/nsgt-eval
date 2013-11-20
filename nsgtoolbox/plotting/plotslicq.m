function plotslicq(c,shift,varargin)
%PLOTSLICQ  |plotnsgtf| wrapper for sliced transforms (sliCQ)
%   Usage:  plotslicq(c,shift,sr,fmin,fmax,bins,cutout,dynrange)
%           plotslicq(c,shift,sr,fmin,fmax,bins,cutout)
%           plotslicq(c,shift,sr,fmin,fmax,bins)
%           plotslicq(c,shift,sr,cutout,dynrange)
%           plotslicq(c,shift,sr,cutout)
%           plotslicq(c,shift,sr)
%           plotslicq(c,shift)
%
%   Input parameters:
%         c        : Array of coefficients.
%         shift    : Vector of frequency shifts of windows
%         sr       : signal sample rate in Hz (default 1 Hz)
%         fmin     : Minimum frequency used in the transform
%         fmax     : Maximum frequency used in the transform
%         bins     : Bins per octave (in constant or vector form)
%         cutout   : Desired part of the spectrogram, e.g.
%                    choice of `2` shows frequencies up to Nyquist
%                    (`X` shows the $number\_of\_bins/X$ lowest frequency 
%                    bins)
%         dynrange : Colorscale dynamic range in dB (default 60 dB)
%
%   This is a wrapper function for |plotnsgtf| that rearranges the 
%   coefficients of a sliced nonstationary Gabor filterbank, in particular
%   those of a sliced constant-Q nonstationary Gabor filterbank, to
%   resemble a full length transform.
%
%   For an explanation of the parameters, please refer to the help of
%   |plotnsgtf|.
%
%   See also:  plotnsgtf, slicq
%
%   References:  dogrhove12

% Author: Nicki Holighaus
% Date: 26.04.13

if nargin < 2
    error('Not enough input arguments');
end

% Rearrange the coefficients, summing up the slices with half overlap
% Recall that the first quarter of the first slice, as well as the last
% quarter of the final slice contain information corresponding to 'the
% other end' of the analyzed signal.

if iscell(c)==0 % Matrix coefficients
    
    [M,Lc,slices] = size(c);
    
    full_cell = zeros(M*slices/2,Lc);
    
    full_cell([end-M/4+1:end,1:3*M/4],:) = c(:,:,1);
    for jj = 1:slices-2
        full_cell(jj*M/2+[-M/4+1:3*M/4],:) = ...
            full_cell(jj*M/2+[-M/4+1:3*M/4],:) + c(:,:,jj+1);
    end
    full_cell([end-3*M/4+1:end,1:M/4],:) = ...
        full_cell([end-3*M/4+1:end,1:M/4],:) + c(:,:,slices);
    
else % Cell array coefficients
    
    Lc = size(c,1);
    slices = size(c{1},2);
    M = cellfun(@(x) length(x(:,1)),c);
    full_cell = cell(Lc,1);
    
    for kk = 1:Lc
        temp = zeros(slices*M(kk)/2,1);
        temp([end-M(kk)/4+1:end,1:3*M(kk)/4]) = c{kk}(:,1);
        for jj = 1:slices-2
            temp(jj*M(kk)/2+(-M(kk)/4+1:3*M(kk)/4)) = ...
                temp(jj*M(kk)/2+(-M(kk)/4+1:3*M(kk)/4)) + c{kk}(:,jj+1);
        end
        temp([end-3*M(kk)/4+1:end,1:M(kk)/4]) = ...
            temp([end-3*M(kk)/4+1:end,1:M(kk)/4]) + c{kk}(:,slices);
        full_cell{kk} = temp;
    end
end

% Plot the rearranged coefficients as spectrogram.

plotnsgtf(full_cell,shift*slices/2,varargin{:});