function g = blackharr(x,L,mod)

if nargin < 1
    error('Not enough input arguments');
end

if numel(x) == 1
    N = x; 
    if nargin < 3
        mod = 1;
        if nargin < 2
           L = N;
        end
    end    
    x = [0:ceil(N/2)-1,-N*ones(1,L-N),-floor(N/2):-1]'/N;
elseif nargin < 2
    mod = 1;
end

if size(x,2) > 1
    x = x.';
end

x = x+1/2;
if mod == 0
    g = 0.35875 - 0.48829*cos(2*pi*x) + 0.14128*cos(4*pi*x) - 0.01168*cos(6*pi*x);
else
    g = 0.35872 - 0.48832*cos(2*pi*x) + 0.14128*cos(4*pi*x) - 0.01168*cos(6*pi*x);
end
g = g .* (x > 0) .* (x < 1);