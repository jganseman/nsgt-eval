function [c,Ls,res,Nit]=nsgaiterf(f,g,shift,M,varargin)
%NSGAITERF  Iterative nonstationary Gabor filterbank analysis
%   Usage:  [c,Ls,res,Nit]=nsgaiterf(f,g,shift,M,varargin)
%           [c,Ls,res,Nit]=nsgaiterf(f,g,shift,M)
%           [c,Ls,res]=nsgaiterf(...)
%           [c,Ls]=nsgaiterf(...)
%           c=nsgaiterf(...)
%
%   Input parameters:
%         f         : Input signal
%         g         : Cell array of filters
%         shift     : Vector of shifts between the center frequencies
%         M         : Vector of lengths of the filters
%         varargin  : Optional input pairs (see table below)
%   Output parameters: 
%         c         : Filterbank coefficients
%         Ls        : Input signal length
%         res       : Vector of relative residuals
%         Nit       : Number of iterations
%
%   Help text goes here.
%
%   Optional input arguments arguments can be supplied like this::
%
%       nsgaiterf(f,g,shift,M,'tol',tol)
%
%   The arguments must be character strings followed by an
%   argument:
%
%     'tol',tol      Error tolerance
%
%     'maxit',maxit      Maximum number of iterations
%
%     'prec',prec    Preconditioning switch
%
%   See also:  nsigtf, nsgsiterf
%

% Author: Nicki Holighaus
% Date: 04.03.13

if nargin < 2
    error('Not enough input arguments');
end

% Set default parameters
tol = 10^-10;   % Error tolerance
maxit = 200;      % Maximum number of iterations
prec = 0;

if nargin >= 3
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
            case {'maxit'}
                maxit = varargin{kk+1};
            case {'prec'}
                prec = varargin{kk+1};
            otherwise
                error(['Invalid input argument: ', varargin{kk}]);
        end
    end
end

Ls = length(f);
N = length(shift);

frmop = @(x) nsigtf(nsgtf(x,g,shift,M),g,shift,Ls);

if prec == 0
    [f,tmp1,tmp2,Nit,res] = pcg(frmop,f,tol,maxit);
else 
    % Construct the diagonal of the frame operator matrix explicitly
    diagonal=zeros(Ls,1);
    for ii = 1:N
        Lg = length(g{ii});

        win_range = mod(timepos(ii)+(-floor(Lg/2):ceil(Lg/2)-1),Ls)+1;
        diagonal(win_range) = diagonal(win_range) + ...
            (fftshift(g{ii}).^2)*M(ii);   
    end
    D = spdiags(diagonal,0,Ls,Ls);
    [f,tmp1,tmp2,Nit,res] = pcg(frmop,f,tol,maxit,D);    
end

c = nsgtf(f,g,shift,M);

if nargout>1
   res=res/norm(f(:));
end