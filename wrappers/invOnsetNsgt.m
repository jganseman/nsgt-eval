function fr = invOnsetNsgt(c,g,shift,M,Ls)

if nargin < 5
    if nargin < 4
        error('Not enough input arguments');
    end
    Ls = sum(shift);
end
    
gd = nsdual(g,shift,M);

fr = nsigt_real(c,gd,shift,M,Ls);