function [g,shift,M] = nsgerbwin(bins,sr,Ls,varargin)
%NSGERBWIN  ERBlet dictionary generator
%   Usage:  [g,shift,M]=nsgerbwin(bins,sr,Ls,varargin)
%           [g,shift,M]=nsgerbwin(bins,sr,Ls)
%
%   Input parameters: 
%         bins      : Desired bins per ERB
%         sr        : Sampling rate of f (in Hz)
%         Ls        : Signal length
%         varargin  : Optional input pairs (see table below)
%   Output parameters: 
%         g         : Cell array of ERBlet filters
%         shift     : Vector of frequency shifts
%         M         : Number of time channels
%
%   Create a nonstationary Gabor filterbank composed of filters regularly 
%   spaced on the ERB frequency scale and having constant Equivalent 
%   Rectangular Bandwidth. 
%
%   The conversion formula of Hz to ERB number is given by
%
%   .. ERBnum(x) = 9.2645*sign(x)*log(1+abs(x)*0.00437).
%
%   .. math:: ERB_{num}(x) = 9.2645\operatorname{sgn}(x)\log(1+0.00437|x|).
%
%   The Equivalent Rectangular Bandwidth at frequency `x` is 
%
%   .. ERB(x) = 24.7*(0.00437*x + 1).
%
%   .. math:: ERB(x) = 24.7(1 + 0.00437x).
%
%   The filters are chosen symmetrically around the zero frequency and
%   finally a symmetric filter is placed on the Nyquist frequency.
%
%   The result can serve as input parameters for |nsgtf| to obtain the
%   ERBlet analysis coefficients or |nsigtf| to synthesize from
%   coefficients, as well as their counterparts for real-valued signals.
%
%   Optional input arguments arguments can be supplied like this::
%
%       nsgerbwin(bins,sr,Ls,'Qvar',Qvar)
%
%   The arguments must be character strings followed by an
%   argument:
%
%     'Qvar',Qvar              Bandwidth variation factor
%
%     'bwfac',bwfac            Channel numbers *M* are rounded to multiples 
%                              of this
%
%     'winfun',winfun          String containing the desired window 
%                              function name
%
%   See also:  nsgtf, nsgtf_real, winfuns
%
%   References: badohojave11 nebahoso13 

% Author: Thibaud Necciari, Nicki Holighaus
% Date: 25.04.13

% Set defaults
Qvar = 1;
bwfac = 1;
winfun = 'modblackharr';

% Check input arguments

if nargin < 3
    error('Not enough input arguments');
end

if nargin >= 4
    Lvar = length(varargin);
    if mod(Lvar,2)
        error('Invalid input argument');
    end
    for kk = 1:2:Lvar
        if ~ischar(varargin{kk})
            error('Invalid input argument');
        end
        switch varargin{kk}
            case {'Qvar'}
                Qvar = varargin{kk+1};
            case {'bwfac'}
                bwfac = varargin{kk+1};
            case {'winfun'}
                winfun = varargin{kk+1};
            otherwise
                error(['Invalid input argument: ', varargin{kk}]);
        end
    end
end

df = sr/Ls; % frequency resolution in the FFT

fmin = 0;
fmax = sr/2;

% Convert fmin and fmax into ERB
erblims = 9.2645*sign([fmin,fmax]).*log(1+abs([fmin,fmax])*0.00437);

% Determine number of freq. channels
Nf = bins*ceil(erblims(2)-erblims(1));

% Determine center frequencies
erbs = linspace(erblims(1),erblims(2),Nf)';
fc = (1/0.00437)*sign(erbs).*(exp(abs(erbs)/9.2645)-1);
fc(1)=fmin;
fc(end)=fmax;

% Concatenate "virtual" frequency positions of negative-frequency windows
fc = [fc ; flipud(fc(1:end-1))];

gamma = 24.7*(4.37*fc*1E-3 +1); % ERB scale

% Convert center frequencies in Hz into samples

posit = round(fc/df);% Positions of center frequencies in samples
posit(Nf+1:end) = Ls-posit(Nf+1:end);% Extension to negative freq.
shift = [Ls-posit(end); diff(posit)];% Hop sizes in samples

% Compute desired essential (Gaussian) support for each filter
Lwin = 4*round(gamma/df)*Qvar;

% Blackman-Harris windows are slightly broader than Gaussians, this is
% offset by the factor 1.1

M = round(Lwin/1.1);

% Compute cell array of analysis filters
g = arrayfun(@(x) winfuns(winfun,x)/sqrt(x),M,'UniformOutput',0);

M = bwfac*ceil(M/bwfac);

g{1}=1/sqrt(2)*g{1};
g{end}=1/sqrt(2)*g{end};