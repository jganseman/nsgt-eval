function [c,res,Nit]=nsgaiterf(f,g,shift,M,varargin)
%NSGAITERF  Iterative nonstationary Gabor filterbank analysis
%   Usage:  [c,res,Nit]=nsgaiterf(f,g,shift,M,varargin)
%           [c,res,Nit]=nsgaiterf(f,g,shift,M)
%           [c,res]=nsgaiterf(...)
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
%         res       : Vector of relative residuals
%         Nit       : Number of iterations
%
%   Optional input parameters:
%         ['tol',tol]               : Error tolerance
%         ['Mit',Mit]               : Maximum number of iterations
%         ['prec',prec]             : Preconditioning switch
%
%   Help text goes here.
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

Ls = length(f);
N = length(shift);

frmop = @(x) nsigtf(nsgtf(x,g,shift,M),g,shift,Ls);

if prec == 0
    [f,tmp1,tmp2,Nit,res] = pcg(frmop,f,tol,Mit);
else 
    % Construct the diagonal of the frame operator matrix explicitly
    diagonal=zeros(Ls,1);
    for ii = 1:N
        Lg = length(g{ii});

        win_range = mod(timepos(ii)+(-floor(Lg/2):ceil(Lg/2)-1),Ls)+1;
        diagonal(win_range) = diagonal(win_range) + (fftshift(g{ii}).^2)*M(ii);   
    end
    D = spdiags(diagonal,0,Ls,Ls);
    [f,tmp1,tmp2,Nit,res] = pcg(frmop,f,tol,Mit,D);    
end

c = nsgtf(f,g,shift,M);

if nargout>1
   res=res/norm(f(:));
end