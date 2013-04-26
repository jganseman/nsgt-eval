function [g,shift,M,fb] = nsgwvltwin(fmin,bw,bins,sr,Ls,winfun)
%NSGWVLTWIN  Wavelet dictionary generator
%   Usage:  [g,shift,M,fb] = nsgwvltwin(fmin,bw,bins,sr,Ls,winfun)
%           [g,shift,M,fb] = nsgwvltwin(fmin,bw,bins,sr,Ls)
%           [g,shift,M] = nsgwvltwin(...)
%
%   Input parameters: 
%         fmin      : Desired minimum center frequency (in Hz)
%         bw        : Desired bandwidth in the first frequency band
%                     (in Hz)
%         bins      : Desired number of bins per octave
%         sr        : Sampling rate of f (in Hz)
%         Ls        : signal length
%         winfun    : String containing the window fucntion name, e.g. the 
%                     following are available:
%                     'hann'      - Hann window (default),
%                     'blackharr' - Blackman-Harris window,
%                     'gauss'     - Truncated Gaussian window,
%                     'wp2inp'    - Uncertainty minimizer
%   Output parameters:
%         g         : Cell array of Fourier transforms of the analysis
%                     Wavelets
%         shift     : Vector of frequency shifts
%         M         : Number of time channels
%         fb        : Frame bounds of the resulting system
%
%   Given the parameter set *fmin*, *bw*, *bins*, *sr* and *Ls*, this
%   function constructs a painless system of bandlimited Wavelets spanning 
%   the range of frequencies from *fmin* to $2^{k/bins} fmin$ with $k$ such 
%   that the dilate centered at this frequency will be completely contained 
%   in the positive frequencies, but the Wavelet at $2^{k+1/bins} fmin$
%   would not be. The number of scales per octave is determined by the 
%   input parameter *bins*, while the Wavelet corresponding to the largest
%   (time-)scale will is constructed to have a bandwidth of *bw* Hz. The 
%   low and high frequencies will be spanned by a plateau-like filter each 
%   to ensure the frame property. 
%
%   If you are not familiar with Wavelet systems, please use |wvlttrans| 
%   instead.
%
%   See also:  wvlttrans, invwvlttrans, winfuns
%
%   References: badohojave11

% Author: Christoph Wiesmeyr, Nicki Holighaus
% Date: 25.04.13

% Check input parameters

if nargin < 6
    winfun = 'hann';
    if nargin < 5
	error('Not enough input arguments');
    end
end

%% Set preliminaries
a = 2^(1/bins);
fmin = fmin*Ls/sr;
bw = bw*Ls/sr;
delta = 1/fmin;

bwd = bw*delta;

gamma = 1/2*(bwd+sqrt(bwd^2+4));
k = .5 * log(a)/log(gamma);

bl=1/gamma;
br=gamma;

bl=bl/delta;
br=br/delta;
BL = bl;
BR = br;

% Determine number of scales
scales = floor(log2(Ls/2/br)*bins)+1;

% Prepare Wavelet system container
g = cell(2*scales+2,1);
M = zeros(2*scales+2,1);
posit = M;

%% Compute the Wavelet system parameters for all scales
pow2 = 2.^((0:scales-1)/bins)';
bl = bl*pow2;
br = br*pow2;
bw = bw*pow2;

% The translation factors are taken to be the bandwidth (in samples)
% rounded up

M(2:scales+1) = ceil(bw);
M(end:-1:scales+3) = M(2:scales+1);

% Cunstruct the Wavelets as log-warped functions of type 'winfun'

points = arrayfun(@(x,y) ceil(x)*delta:delta:y*delta,bl,br,...
    'UniformOutput',0);
points = arrayfun(@(x) k*(log(points{x})/log(a)-x+1),1:scales,...
    'UniformOutput',0);

g(2:scales+1) = cellfun(@(x) winfuns(winfun,x),points,'UniformOutput',0);
g(2:scales+1) = arrayfun(@(x) g{x}/sqrt(M(x)),2:scales+1,...
    'UniformOutput',0);
g(end:-1:scales+3) = cellfun(@(x) flipud(x),g(2:scales+1),...
    'UniformOutput',0);
g = cellfun(@(x) ifftshift(x),g,'UniformOutput',0);

Lg = cellfun(@length,points)';

% Construct positions and relative shifts of the Wavelets

posit(1) = 1;
posit(scales+2) = floor(Ls/2)+1;
posit(2:scales+1) = ceil(bl) + floor(Lg/2)+1;
posit(end:-1:scales+3) = Ls-posit(2:scales+1) ...
                                + 1 + mod(Lg-1,2);

shift = [Ls-mod(posit(end)-1,Ls);diff(posit)];

%% Compute 'scaling' or 'padding' functions for the 0- and
%% Nyquist-frequency

% Construct the diagonal of the 'incomplete' frame operator matrix
% explicitly

diagonal=zeros(Ls,1);
for ii = 2:scales+1
  win_range = mod(posit(ii)+...
      (-floor(length(g{ii})/2):ceil(length(g{ii})/2)-1)-1,Ls)+1;
  diagonal(win_range) = diagonal(win_range) + M(ii)*fftshift(g{ii}).^2 ;   
end
for ii = scales+3:2*scales+2
  win_range = mod(posit(ii)+...
      (-floor(length(g{ii})/2):ceil(length(g{ii})/2)-1)-1,Ls)+1;
  diagonal(win_range) = diagonal(win_range) + M(ii)*fftshift(g{ii}).^2 ;
end

% Compute the padding functions ('scaling functions'), this functions will
% be chosen such that they work well with a window set having 3/4 overlap.
% In this case, the result should be a tight frame.
% In other cases of non-zero overlap, the padding windows will produce a
% frame, but might be non-smooth.
 
 maxX = max(diagonal);
 
 LowLim = ceil(BL*2^(3/bins));
 UppLim = ceil(Ls/2-BR*2^((scales-4)/bins));
 
 M(1) = 2*LowLim;
 g{1} = sqrt((maxX-diagonal([1:LowLim,end-LowLim+1:end]))/M(1));
 
 M(scales+2) = 2*UppLim;
 g{scales+2} = sqrt((maxX-diagonal(floor(Ls/2)+[1:UppLim,-UppLim+1:0]))...
     /M(scales+2));
 
 %% Compute the frame bounds if output parameter 'fb' is desired
 if nargout == 4
 
    % Finalize the construction of the diagonal 
    Lg = cellfun(@length,g);
    
    for ii = [1,scales+2]
      range = mod(posit(ii)+(-floor(Lg(ii)/2):ceil(Lg(ii)/2)-1)-1,Ls)+1;
      diagonal(range) = diagonal(range) + (fftshift(g{ii}).^2)*M(ii);
    end
    
    % Since the system is painless, the minimum and maximum of the diagonal
    % equal the frame bounds
    
    fb = [min(diagonal) max(diagonal)];
 end