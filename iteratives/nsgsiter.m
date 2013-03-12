function [fr,res,Nit]=nsgsiter(c,g,shift,M,varargin)
%NSGSITER  Iterative nonstationary Gabor synthesis
%   Usage:  [fr,res,Nit]=nsgsiter(c,g,shift,M,varargin)
%           [fr,res,Nit]=nsgsiter(c,g,shift,M)
%           [fr,res]=nsgsiter(...)
%           fr=nsgsiter(...)
%
%   Input parameters:
%         c         : Nonstationary Gabor coefficients
%         g         : Cell array of window functions
%         shift     : Vector of shifts between the window positions
%         M         : Vector of lengths of the window functions
%         varargin  : Optional input pairs (see table below)
%   Output parameters: 
%         fr        : Synthesized output signal
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

if nargin < 4
    error('Not enough input arguments');
end

% Set default parameters
tol = 10^-10;   % Error tolerance
Mit = 200;      % Maximum number of iterations
prec = 0;

if nargin >= 4
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

frmop = @(x) nsigt(nsgt(x,g,shift,M),g,shift,L);

fr = nsigt(c,g,shift,L);

if prec == 0
    [fr,tmp1,tmp2,Nit,res] = pcg(frmop,fr,tol,Mit);
else 
    % Construct the diagonal of the frame operator matrix explicitly
    diagonal=zeros(Ls,1);
    for ii = 1:N
        Lg = length(g{ii});

        win_range = mod(posit(ii)+(-floor(Lg/2):ceil(Lg/2)-1),Ls)+1;
        diagonal(win_range) = diagonal(win_range) + (fftshift(g{ii}).^2)*M(ii);   
    end
    D = spdiags(diagonal,0,Ls,Ls);
    [fr,tmp1,tmp2,Nit,res] = pcg(frmop,fr,tol,Mit,D);    
end

if  exist('Ls','var')
    fr = fr(1:Ls,:);
end

if nargout>1
   res=res/norm(fr(:));
end