function [fr,res,Nit]=nsgsiterf(c,g,shift,M,varargin)
%NSGSITERF  Iterative nonstationary Gabor filterbank synthesis
%   Usage:  [fr,res,Nit]=nsgsiterf(c,g,shift,M,Ls,varargin)
%           [fr,res,Nit]=nsgsiterf(c,g,shift,M,varargin)
%           [fr,res,Nit]=nsgsiterf(c,g,shift,M,Ls)
%           [fr,res,Nit]=nsgsiterf(c,g,shift,M)
%           [fr,res]=nsgsiterf(...)
%           fr=nsgsiterf(...)
%
%   Input parameters:
%         c         : Nonstationary Gabor coefficients
%         g         : Cell array of filters
%         shift     : Vector of shifts between the center frequencies
%         M         : Vector of lengths of the filters
%         Ls        : Original signal length
%         varargin  : Optional input pairs (see table below)
%   Output parameters: 
%         fr        : Synthesized output signal
%         res       : Vector of relative residuals
%         Nit       : Number of iterations
%
%   Help text goes here.
%
%   Optional input arguments arguments can be supplied like this::
%
%       nsgsiter(c,g,shift,M,'tol',tol)
%
%   The arguments must be character strings followed by an
%   argument:
%
%     'tol',tol      Error tolerance
%
%     'Mit',Mit      Maximum number of iterations
%
%     'prec',prec    Preconditioning switch
%

% Author: Nicki Holighaus
% Date: 04.03.13

if nargin < 2
    error('Not enough input arguments');
end

% Set default parameters
tol = 10^-10;   % Error tolerance
Mit = 200;      % Maximum number of iterations
prec = 0;

if nargin >= 3
    if isnumeric(varargin{1}) && numel(varargin{1}) == 1
        Ls = varargin{1};
        varargin = varargin(2:end);
    end
    Lvar = length(varargin);
    if mod(Lvar,2)
        error('Invalid input argument');
    end
    for kk = 1:2:Lvar
        if ~ischar(varargin{kk})
            error('Invalid input argument');
        end
        switch varargin{kk}
            case {'tol'}
                tol = varargin{kk+1};
            case {'Mit'}
                Mit = varargin{kk+1};
            case {'prec'}
                prec = varargin{kk+1};
            otherwise
                error(['Invalid input argument: ', varargin{kk}]);
        end
    end
end

N = length(shift);
posit = cumsum(shift);
L = posit(end);
posit = posit-shift(1);

frmop = @(x) nsigtf(nsgtf(x,g,shift,M),g,shift,L);

fr = nsigtf(c,g,shift,L);

if prec == 0
    [fr,tmp1,tmp2,Nit,res] = pcg(frmop,fr,tol,Mit);
else 
    % Construct the diagonal of the frame operator matrix explicitly
    diagonal=zeros(L,1);
    for ii = 1:N
        Lg = length(g{ii});

        win_range = mod(posit(ii)+(-floor(Lg/2):ceil(Lg/2)-1),L)+1;
        diagonal(win_range) = diagonal(win_range) + ...
            (fftshift(g{ii}).^2)*M(ii);   
    end
    D = spdiags(diagonal,0,L,L);
    [fr,tmp1,tmp2,Nit,res] = pcg(frmop,fr,tol,Mit,D);    
end

if  exist('Ls','var')
    fr = fr(1:Ls,:);
end

if nargout>1
   res=res/norm(fr(:));
end