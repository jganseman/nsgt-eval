function [g,shift,M] = nsgcqwin(fmin,fmax,bins,sr,Ls,varargin)
%NSGCQWIN  Constant-Q/Variable-Q dictionary generator
%   Usage:  [g,shift,M] = nsgcqwin(fmin,fmax,bins,sr,Ls,varargin)
%           [g,shift,M] = nsgcqwin(fmin,fmax,bins,sr,Ls)
%
%   Input parameters:
%         fmin      : Minimum frequency (in Hz)
%         fmax      : Maximum frequency (in Hz)
%         bins      : Vector consisting of the number of bins per octave
%         sr        : Sampling rate (in Hz)
%         Ls        : Length of signal (in samples)
%         varargin  : Optional input pairs (see table below)
%   Output parameters: 
%         g         : Cell array of window functions
%         shift     : Vector of shifts between the center frequencies
%         M         : Vector of lengths of the window functions
%
%   Create a nonstationary Gabor filterbank with constant or varying 
%   Q-factor and relevant frequency range from *fmin* to *fmax*. To allow
%   for perfect reconstruction, the frequencies outside that range will be
%   captured by 2 additional filters placed on the zero and Niyquist
%   frequencies, respectively.
%
%   The Q-factor (quality factor) is the ratio of center frequency to
%   bandwidth `cent_freq/bandwidth`.
%
%   To create a constant-Q filterbank with a fixed number of bins per 
%   octave, use a scalar parameter *bins*. The default parameters serve to
%   set up a filter sequence with approximately 1/2 overlap and only
%   approximately constant Q-factor (up to 1 sample deviation). The 
%   optional switch *fractional* can be set to 1 to allow for fractional 
%   sampling and exactly constant Q-factor.
%
%   Alteratively, a vector *bins* can be supplied. In this case, successive
%   octaves can have different numbers of filters regularly spaced on a
%   logarithmic scale, e.g. *bins(1)* filters will be placed between `fmin`
%   and `2*fmin`, *bins(2)* filters between `2*fmin` and `4*fmin` and so
%   on.
%
%   For more details on the construction of the constant-Q nonstationary 
%   Gabor filterbank, please check the reference.
%   
%   Optional input arguments arguments can be supplied like this::
%
%       nsgcqwin(fmin,fmax,bins,sr,Ls,'min_win',min_win)
%
%   The arguments must be character strings followed by an
%   argument:
%
%     'min_win',min_win  Minimum admissible window length (in samples) 
%
%     'Qvar',Qvar              Bandwidth variation factor
%
%     'bwfac',bwfac            Channel numbers *M* are rounded to multiples 
%                              of this
%
%     'fractional',fractional  Allow fractional shifts and bandwidths
%
%     'winfun',winfun          String containing the desired window 
%                              function name
%
%   See also:  nsgtf, nsgtf_real, winfuns
%
%   References: dogrhove11 dogrhove12

% Authors: Nicki Holighaus, Gino Velasco, Monika Doerfler
% Date: 25.04.13

% Set defaults
Qvar = 1;
bwfac = 1;
min_win = 4;
fractional = 0;
winfun = 'hann';

% Check input arguments

if nargin < 5
    error('Not enough input arguments');
end

if nargin >= 6
    Lvar = length(varargin);
    if mod(Lvar,2)
        error('Invalid input argument');
    end
    for kk = 1:2:Lvar
        if ~ischar(varargin{kk})
            error('Invalid input argument');
        end
        switch varargin{kk}
            case {'min_win'}
                min_win = varargin{kk+1};
            case {'Qvar'}
                Qvar = varargin{kk+1};
            case {'bwfac'}
                bwfac = varargin{kk+1};
            case {'fractional'}
                fractional = varargin{kk+1};
            case {'winfun'}
                winfun = varargin{kk+1};
            otherwise
                error(['Invalid input argument: ', varargin{kk}]);
        end
    end
end

nf = sr/2;

if fmax > nf
    fmax = nf;
end

b = ceil(log2(fmax/fmin))+1;

if length(bins) == 1;
    bins = bins*ones(b,1);
elseif length(bins) < b
    if size(bins,1) == 1
        bins=bins.';
    end
    bins( bins<=0 ) = 1;
    bins = [bins ; bins(end)*ones(b-length(bins),1)];
end

fbas = zeros(sum(bins)+1,1);

ll = 0;
for kk = 1:length(bins);
    fbas(ll+(1:bins(kk)+1)) = ...
        fmin*2.^(((kk-1)*bins(kk):(kk*bins(kk))).'/bins(kk));
    ll = ll+bins(kk);
end

temp = find(fbas >= fmax,1);
if fbas(temp-1) + (fbas(temp) - fbas(temp-2))/2 >= nf
    fbas = fbas(1:temp-2);
elseif fbas(temp) + (fbas(temp+1) - fbas(temp-1))/2 >= nf
    fbas = fbas(1:temp-1);
else 
    fbas = fbas(1:temp); 
end

Lfbas = length(fbas);

fbas = [0;fbas];
fbas(Lfbas+2) = nf;
fbas(Lfbas+3:2*(Lfbas+1)) = sr-fbas(Lfbas+1:-1:2);

fbas = fbas*(Ls/sr);

% Set bandwidths

bw = zeros(2*Lfbas+2,1);

bw(1) = 2*fmin*(Ls/sr);
bw(2) = (fbas(2))*(2^(1/bins(1))-2^(-1/bins(1)));

for k = [3:Lfbas , Lfbas+2]
    bw(k) = (fbas(k+1)-fbas(k-1));
end

bw(Lfbas+1) = (fbas(Lfbas+1))*(2^(1/bins(end))-2^(-1/bins(end)));
bw(Lfbas+3:2*Lfbas+2) = bw(Lfbas+1:-1:2);

posit = zeros(size(fbas));
posit(1:Lfbas+2) = floor(fbas(1:Lfbas+2));
posit(Lfbas+3:end) = ceil(fbas(Lfbas+3:end));

shift = [mod(-posit(end),Ls); diff(posit)];

bw = Qvar*bw;

if fractional
    corr_shift = fbas-posit;
    M = ceil(bw+1);
else
    bw = round(bw);
    M = bw;
end

for ii = 1:2*(Lfbas+1)
    if bw(ii) < min_win;
        bw(ii) = min_win;
        M(ii) = bw(ii);
    end
end

if fractional
    g = arrayfun(@(x,y,z) ...
        winfuns(winfun,([0:ceil(z/2),-floor(z/2):-1]'-x)/y)/sqrt(y),corr_shift,...
        bw,M,'UniformOutput',0);
else
    g = arrayfun(@(x) winfuns(winfun,x)/sqrt(x),...
        bw,'UniformOutput',0);
end

M = bwfac*ceil(M/bwfac);

% Setup Tukey window for 0- and Nyquist-frequency
for kk = [1,Lfbas+2]
    if M(kk) > M(kk+1);
        g{kk} = ones(M(kk),1);
        g{kk}((floor(M(kk)/2)-floor(M(kk+1)/2)+1):(floor(M(kk)/2)+...
            ceil(M(kk+1)/2))) = winfuns('hann',M(kk+1));
        g{kk} = g{kk}/sqrt(M(kk));
    end
end