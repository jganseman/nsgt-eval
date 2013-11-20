function plot_wins(g,shift,normalize)
%PLOT_WINS  Plot nonstationary Gabor windows/filters
%   Usage:  plot_wins(g,shift)
%           plot_wins(g,shift,normalize)
%
%   Input parameters:
%         g         : Cell array of windows/filters
%         shift     : Vector of time/frequency shifts
%         normalize : Re-normalize the windows to have approximately 
%                     uniform height
%
%   This helper function plots the distribution of the windows/filters of a
%   nonstationary Gabor system/filterbank along the time/frequency axis. 
%   The shape of the windows/filters is determined from the cell array *g* 
%   and their position on the respective axis from the position vector 
%   *shift*.
%
%   See also:  nsgsclwin, nsgwvltwin, nsgerbwin, nsgcqwin

% Author:  Nicki Holighaus
% Date: 26.04.13

if nargin < 3
    normalize = 0;
end

N = length(shift);
posit = cumsum(shift)-shift(1);

% Every second window is plotted in red, the other ones in blue

color = ['b', 'r'];

% This loop just plots each window at its corresponding time position.

for ii = 1:N
    Lg = length(g{ii});
    
    win_range = posit(ii)+(-floor(Lg/2):ceil(Lg/2)-1);
    if normalize == 1 % If normalize is set to 1, normalize
        % the windows to have the same maximum
        plot(win_range, fftshift(g{ii}).*sqrt(length(g{ii})), ...
            color(rem(ii,2)+1));
    else
        plot(win_range, fftshift(g{ii}), ...
            color(rem(ii,2)+1));
    end
    hold on;
end
hold off; shg
