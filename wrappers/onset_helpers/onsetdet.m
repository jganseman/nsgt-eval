function [pos,V0] = onsetdet(f,win_length,thre,range,multi,shift,showplot)
%ONSETDET  Onset detection wrapper
%   Usage: [pos,V0] = onsetdet(f,win_length,thre,range,multi,shift,showplot)
%          [pos,V0] = onsetdet(f,win_length,thre,range,multi,shift)
%          [pos,V0] = onsetdet(f,win_length,thre,range,multi)
%          [pos,V0] = onsetdet(f,win_length,thre,range)
%          [pos,V0] = onsetdet(f,win_length,thre)
%          [pos,V0] = onsetdet(f,win_length)
%          pos = onsetdet(...)
%
%   Input parameters:
%         f         : The signal to be analyzed (single channel only)
%         win_length: Window length for the STFT analysis (in samples)
%         thre      : Peak-picking threshold
%         range     : Area of interest for the choice of local maxima 
%         multi     : Area of interest multiplication factor for the   
%                     peak-picking
%         shift     : Readjustment of the peaks 
%                     (in $shift*win_length/16$)
%         showplot  : Plot the results (0/1) 
%   Output parameters:
%         pos       : Onset sequence
%         V0        : Regular discrete Gabor transform of f
%
%   Uses routines from LTFAT 0.97 or higher, available at:
%   http://ltfat.sourceforge.net/
%
%   This routine produces a sequence of onsets using a straightforward
%   realization of the onset detection procedure described in 
%
%   S. Dixon. Onset Detection Revisited. Proc. DAFx 06, pages 133-137, 
%   September 2006.
%
%   Uses a Hann window STFT analysis and the spectral flux onset detection
%   function.
%
%   A local maximum of the onset detection function is chosen as onset
%   if it is larger than the local mean by at least *thre*. 
%   A time slice is considered a local maximum if its spectral flux value 
%   is larger than those of the surrounding slices on an area of $+-range$.
%   The local mean is computed as the mean value of the spectral flux
%   function on an area corresponding to $-multi*range$ to $+range$ of the 
%   current position.
%
%   External: DGT (LTFAT routine)

% Author: Nicki Holighaus
% Date: 04.03.13

% Check input arguments and set default arguments if necessary

if nargin < 7
    showplot = 0;
    if nargin < 6
        shift = 0;
        if nargin < 5
            multi = 3;
            if nargin < 4
                range = 3;
                if nargin < 3
                    error('Not enough input arguments');
                end
            end
        end
    end
end

[Ls,col] = size(f);

if min(Ls,col) > 1
    error('Right now, this routine supports only single channel signals');
end

if ( col ~= 1 && Ls == 1 )
    f = f.';
    Ls = col;
end

clear Ls col;

% Compute the spectral flux

[ODF,V0]=specflux(f,win_length,win_length/16);

% Select the significant paeks in the spectral
% flux

pos = peakpick(ODF,thre,range,multi);

% Shift the onset positions according by a fixed amount to be more precise
% (experimental, but improves the results on simple signals)

pos = onsets(pos,tgap,win_length,shift);

% Due to periodization and shifts, some onsets might appear after the
% end of the signal. Those are omitted

X = length(pos);
if ( X > 0 )
  while ( pos(X) >= length(f))
      X = X-1;
  end
end

% The first sample is always considered an onset

pos = [1,pos(1:X)].';

if showplot ~= 0    % Plot the results

    pos2 = floor(1+(pos-1)/tgap);
    g=max(ODF)*ones(length(ODF),1);
    g(pos2) = min(ODF)*ones(length(pos2),1);

    subplot(1,2,1);
    imagesc(20*log10(abs(V0(size(V0,1)*4/8+1:size(V0,1),:))+eps));
    vline(pos2,'k');
    subplot(1,2,2);
    plot(ODF,'k'); hold on; plot(g,'r+');
    axis tight; hold off; shg
end
    
end

% onsets should shift the detected onsets uniformly
% to another position. In many cases, this might improve
% results significantly.
% 
% The real onsets do not always exactly coincide with
% the chosen peaks, but for simple signals they should 
% be around (peak + some constant times [+-tgap]).
%
% Right now, this is just an experimental routine
% and still a work in progress. Work has to be done.

function pos = onsets(peaks,tgap,win_length,shift)

pos = 1+(peaks-1).*tgap + floor(shift*win_length/16);

end
