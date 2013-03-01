function [g,shift,M] = nsgcqwin(fmin,fmax,bins,sr,Ls,varargin)

% NSGCQWIN.M
%---------------------------------------------------------------
% [g,shift,M]=nsgcqwin(fmin,fmax,bins,sr,Ls) creates a set of windows whose
% centers correspond to center frequencies to be
% used for the nonstationary Gabor transform with varying Q-factor. 
%---------------------------------------------------------------
%
% INPUT : fmin ...... Minimum frequency (in Hz)
%	      fmax ...... Maximum frequency (in Hz)
%         bins ...... Vector consisting of the number of bins per octave
%         sr ........ Sampling rate (in Hz)
%         Ls ........ Length of signal (in samples)
%         min_win.... Minimum admissible window length (in samples) 
%
% OUTPUT : g ......... Cell array of window functions.
%          shift ..... Vector of shifts between the center frequencies.
%          M ......... Vector of lengths of the window functions.
%
%---------------------------------------------------------------
% If min_win is not specified, it is set to 4. Note that while large 
% values for min_win might lead to changing Q-factors among the lower
% bins, small values significantly increase time aliasing.
%---------------------------------------------------------------
%
% AUTHOR(s) : Monika Dörfler, Gino Angelo Velasco, Nicki Holighaus, 2011
%
% EXTERNALS : hannwin 

% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to 
% Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

% Set defaults 
Qvar = 1;
wL_fac = 1;
min_win = 4;
fractional = 0;
winfun = @hannwin;

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
            case {'wL_fac'}
                wL_fac = varargin{kk+1};
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

fbas = zeros(sum(bins),1);

ll = 0;
for kk = 1:length(bins); 
    fbas(ll+(1:bins(kk))) = fmin*2.^(((kk-1)*bins(kk):(kk*bins(kk)-1)).'/bins(kk)); 
    ll = ll+bins(kk);
end

temp = find(fbas>=fmax,1);
if fbas(temp) >= nf 
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

rfbas = zeros(size(fbas));
rfbas(1:Lfbas+2) = floor(fbas(1:Lfbas+2));
rfbas(Lfbas+3:end) = ceil(fbas(Lfbas+3:end));

shift = [mod(-rfbas(end),Ls); diff(rfbas)];

bw = Qvar*bw;

if fractional
    corr_shift = fbas-rfbas;
    M = ceil(bw+1);
else
    bw = round(bw);
    M = bw;
end

N = length(shift);
g = cell(N,1);

for ii = 1:2*(Lfbas+1)    
    if bw(ii) < min_win; 
        bw(ii) = min_win;
        M(ii) = bw(ii);
    end 
end

if fractional 
    g = arrayfun(@(x,y,z) winfun(([0:ceil(z/2),-floor(z/2):-1]'-x)/y)/sqrt(y),...
        corr_shift,bw,M,'UniformOutput',0);
else
    g = arrayfun(@(x) winfun(x)/sqrt(x),...
        bw,'UniformOutput',0);
end

M = wL_fac*ceil(M/wL_fac);

% Setup Tukey window for 0- and Nyquist-frequency
for kk = [1,Lfbas+2]
    if M(kk) > M(kk+1);    
        g{kk} = ones(M(kk),1);     
        g{kk}((floor(M(kk)/2)-floor(M(kk+1)/2)+1):(floor(M(kk)/2)+...
            ceil(M(kk+1)/2))) = hannwin(M(kk+1));
        g{kk} = g{kk}/sqrt(M(kk));
    end
end