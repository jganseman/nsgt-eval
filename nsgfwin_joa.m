function [g,shift,M] = nsgfwin_joa(fmin,fmax,bins,sr,Ls,min_win)

% Create a set of windows and time shift positions for the Nonstationary
% Gabor Transform. This particular implementation creates a grid structure
% (which is overly redundant in most frequencies), aiming to obtain a
% perfectly rasterized rectangular matrix of coefficients.

% This routine is based on ngsfwin.m by Doerfler, Velasco and Holighaus,
% part of the NSGT toolbox, which has been released under CC
% Attribution-NonCommercial-ShareAlike 3.0 Unported license.

% INPUT : fmin,fmax . Minimum and Maximum frequency (in Hz)
%         bins ...... Vector consisting of the number of bins per octave
%         sr ........ Sampling rate (in Hz)
%         Ls ........ Length of signal (in samples)
%         min_win.... Minimum admissible window length (in samples) 
%
% OUTPUT : g ......... Cell array of window functions.
%          shift ..... Vector of shifts between the center frequencies.
%          M ......... Vector of lengths of the window functions.

% EDITED by Joachim Ganseman to create a set of windows
% effectuating a rasterized (oversampled) CQT transform
% except for the Nyquist component!

% The file is heavily commented while I tried to figure out what every line
% of code was actually doing.

if nargin < 6
    min_win = 4;
    if nargin < 5
	error('Not enough input arguments');
    end
end

nf = sr/2; 

if fmax > nf            % don't go higher than the Nyquist frequency
    fmax = nf; 
end

b = ceil(log2(fmax/fmin))+1;    % the number of octaves

if length(bins) == 1;           % create vector of octaves with bins/octave
    bins = bins*ones(b,1);
elseif length(bins) < b
    if size(bins,1) == 1
        bins=bins.'; 
    end
    bins(find(bins<=0)) = 1;
    bins = [bins ; min(bins)*ones(b-length(bins),1)];
end
                                % 'bins' is e.g. [48 48 48 48 48]'

fbas = [];          % fbas are all center frequencies from fmin to fmax

for kk = 1:length(bins); 
    fbas = [fbas;fmin*2.^(((kk-1)*bins(kk):(kk*bins(kk)-1)).'/bins(kk))]; 
end

                    % limit the frequency bins to the Nyquist frequency
if fbas(min(find(fbas>=fmax))) >= nf 
    fbas = fbas(1:max(find(fbas<fmax)));    
else
    fbas = fbas(1:min(find(fbas>=fmax)));
end

lbas = length(fbas);            % nr of bins
fbas = [0;fbas];                % append 0 frequency at beginning
fbas(lbas+2) = nf;              % append Nyquist frequency at end
fbas(lbas+3:2*(lbas+1)) = sr-fbas(lbas+1:-1:2);     % append mirror freqs
        % this line mirrors the frequencies already defined over the
        % Nyquist frequency, as to represent the bins' equivalent between
        % sr/2 and sr. It's the "redundant" part of the transform, which
        % will begin at sr-fmax and end at sr-fmin.
        
  % now 'fbas' is e.g. [0 20 40 80 ... 10000 22050 34100 ... 44060 44080]

fbas = fbas*(Ls/sr);      % multiply by length of signal. Question is: why?

% fill a vector M with 0, same length as the center frequency vector
M = zeros(length(fbas),1);
M(1) = 2*fmin*(Ls/sr);      % first element: 2 x fmin (x signallength)
M(2) = (fbas(2))*(2^(1/bins(1))-2^(-1/bins(1)));
        % second element: fmin x binwidth (x signallength).
        % binwidth goes from previous center frequency to next center frequency (there is overlap). 
        % The next loop fills up M with the widths of all frequency bins,
        % (still multiplied by signallength / samplerate ! )
for k = [3:lbas , lbas+2]
M(k ) = (fbas(k+1)-fbas(k-1));
end
        % do the same for fmax
M(lbas+1) = (fbas(lbas+1))*(2^(1/bins(end))-2^(-1/bins(end)));
        % repeat a similar mirroring procedure for Nyquist -> samplerate freqs
M(lbas+3:2*(lbas+1)) = M(lbas+1:-1:2);
M(end) = M(2);
        % round all of this to the nearest integer
M = round(M);

% I've got a feeling these numbers are the nr of time divisions for any
% frequency. To rasterize, we'd need to bring them all to the value for fmax!
% but: output says they're only windowlengths?



% create windows of these sizes. Now the strange thing is: the size of
% these windows still depends on signal length (since that's factored into M)?
for ii = 1:2*(lbas+1);
    
    if M(ii) < min_win; 
        M(ii) = min_win;
    end 
    g{ii} = hann(M(ii));
end

% create windows for 0 and Nyquist frequencies
for kk = [1,lbas+2]
    if M(kk) > M(kk+1);    
        g{kk} = ones(M(kk),1);     
        g{kk}((floor(M(kk)/2)-floor(M(kk+1)/2)+1):(floor(M(kk)/2)+...
        ceil(M(kk+1)/2))) = hann(M(kk+1));
    end
end

rfbas = round(fbas); % round frequency center values
shift = [mod(-rfbas(end),Ls); diff(rfbas)]; % differences between frequency centers.

% try this: use fmax for all values in M. This makes it seem that all
% windows have this length. Suppose: they'll be zero-padded if shorter???
    % EDIT: yes, they'll be zero-padded as required, thus acting according
    % to https://ccrma.stanford.edu/~jos/sasp/Zero_Phase_Zero_Padding.html

M(1:end) = (fbas(lbas+1))*(2^(1/bins(end))-2^(-1/bins(end)));
%re-add 0 frequency and Nyquist components
%M(1) = 2*fmin*(Ls/sr);
%M(lbas+1) = (fbas(lbas+1))*(2^(1/bins(end))-2^(-1/bins(end)));
M(lbas+2) = (fbas(lbas+3)-fbas(lbas+1));
M = round(M);

% this shift vector is defining starting positions of the windowed parts of
% the signal that compute the CQT. Quite important to know: those starting
% positions are used on the FFT of the entire signal, not on the signal
% itself... Take a look at the plot function to see how time is calculated

% thus: the coefficients of the NSGT are calculated through the FFT first!
% This complicates pinning things on time. Put a breakpoint in the imaging
% routine to see how the image is built up: line per line from bottom up.

% In fact we want the bottom lines to have same resolution as top lines, so
% they should be same window length. Going back, this means all g{ii}
% should have same length.
