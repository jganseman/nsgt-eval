function plotnsgtf(c,shift,sr,fmin,fmax,bins,cutout,dynrange)
%PLOTNSGTF  Plot nonstationary Gabor filterbank coefficients
%   Usage:  plotnsgtf(c,shift,sr,fmin,fmax,bins,cutout,dynrange)
%           plotnsgtf(c,shift,sr,fmin,fmax,bins,cutout)
%           plotnsgtf(c,shift,sr,fmin,fmax,bins)
%           plotnsgtf(c,shift,sr,cutout,dynrange)
%           plotnsgtf(c,shift,sr,cutout)
%           plotnsgtf(c,shift,sr)
%
%   Input parameters:
%         c        : Array of coefficients.
%         shift    : Vector of frequency shifts
%         sr       : signal sample rate in Hz (default 1 Hz)
%         fmin     : Minimum frequency used in the transform
%         fmax     : Maximum frequency used in the transform
%         bins     : Bins per octave (in constant or vector form)
%         cutout   : Desired part of the spectrogram, e.g.
%                    choice of '2' shows frequencies up to Nyquist
%                    ('X' shows the 'number_of_bins/X' lowest frequency 
%                    bins)
%         dynrange : Colorscale dynamic range in dB (default 60 dB)
%
%   This variation of LTFATs plotndgt allows to plot only a part of the 
%   spectrogram obtained from the frequency side version of the
%   non-stationary Gabor transform and also accepting matrix array
%   input. If runs with the necessary input parameters (fmin,fmax,bins) of
%   the transform, it adjusts the plot labels according to those 
%   parameters.
%

% Author:  Gino Velasco, Nicki Holighaus and Radu C. Frunza
% Original code by: Florent Jaillet
% Date: 04.03.13

%   Remainder of the original help file:
%
%   plondgt plots the spectrogram from coefficients computed with the 
%   function ndgt. For more details on the format of the variables c and a 
%   format, please read the ndgt function help.
%
%   plotndgt uses a dB colorscale, and the dynrange value can be used to
%   specify the dynamic of this colorscale, as the produced image uses a 
%   colormap in the interval [chigh-dynrange,chigh], where chigh is the 
%   highest value in the plot.
%
%   Limitation: plotndgt only works for coefficients c obtained from a
%   monochannel signal.
%
%   SEE ALSO:  NDGT
%
% The original file is part of LTFAT version 0.97

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

ticklabels = 1;

if nargin < 8
    % Default value for colorscale dynamic.
    dynrange = 60;
    if nargin < 7 % Default value for frequency cutout.
        cutout = 2;
        if nargin < 6
            ticklabels = 0;
            if nargin == 5
                cutout = fmin;
                dynrange = fmax;
            elseif nargin < 4
                cutout = 2;
                if nargin < 3
                    % Default sampling rate
                    sr = 1;
                    if nargin < 2
                        error('Not enough input arguments');
                    end
                end
            else
                cutout = fmin;
            end
        end
    end
end

N = length(shift);

% For the coefficient output of nsgtf_real, adjust the frequeny
% channels shown, if necessary

if N > size(c,1) && cutout < 2
    cutout = cutout/2;
end

clf %Clear previous figures

timepos=cumsum(shift)-shift(1);

% Compute maximum of the representation for colorscale dynamic handling.
if iscell(c) == 1
    if size(c{1},2) > 1
        error(['Multichannel spectrograms are not supported. Please ',...
            'use ''cellfun(@(x) x(:,k),c,''UniformOutput'',0)'' to ',...
            'select the k-th channel.']);
    end
    temp=cell2mat(c);
else
    if size(c,3) > 1
        error(['Multichannel spectrograms are not supported. Please ',...
            'use ''c(:,:,k)'' to select the k-th channel.']);
    end
    temp=c;
end
ma=20*log10(max(abs(temp(:))));

% Plot the representation: as the sampling grid in the time frequency plane
% is irregular, the representation by done by plotting many images next to
% each other, with one image for each window
if iscell(c) == 1
    hold('on');
    for ii=1:floor((N-1)/cutout)+1
        temp = 20*log10(abs(c{ii})+eps);
        % +eps is here to avoid log of 0
        % Octave cannot plot images that are only one point wide, so we use
        % images that are to points wide
        imagesc([0,timepos(end)+shift(1)-1]/sr,[ii,ii+1],[temp,temp].',...
            [ma-dynrange,ma]);
    end
    hold('off');
    axis('tight');
else
    imagesc([0, (timepos(end)+shift(1)-1)/sr],[1 size(c,2)],...
        20*log10(abs(c)'+eps),[ma-dynrange, ma]);
    select=[0,(timepos(end)+shift(1)-1)/sr,1,round(size(c,2)/cutout)+1];
    axis(select);
    set(gca,'YDir','normal');
end

% Compute YTickLabels for the correct frequencies
if ticklabels == 1
    
    vfq = [0, fmin*2.^(0:log2(fmax/fmin))];
    
    if length(bins) < length(vfq) - 1
        xbins = [bins, bins(end)*ones(1,length(vfq)-length(bins)-1)];
    else
        xbins = bins(1:length(vfq)-1);
    end;
    
    sbins = [0,2,1+cumsum(xbins)];
    
    yTick = [2;2+2*cumsum(xbins(1:end-1))';(1+N/2)];
    yTick = unique([yTick(yTick <= floor((length(shift)-1)/cutout)+1);...
        floor((N-1)/cutout)+1]);
    
    yTickLabel = zeros(length(yTick),1);
    
    for kk = 1:length(yTick)
        ind = find(yTick(kk) < sbins,1);
        yTickLabel(kk) = vfq(ind-1)*...
            2^(((yTick(kk)-1)-sbins(ind-1))/(sbins(ind)-sbins(ind-1)));
    end;
    
    yTickLabel = num2str(round(yTickLabel),5);
    set(gca,'YTick',yTick,'YTickLabel',yTickLabel);
end