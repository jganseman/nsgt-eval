% allows to pad a window obtained from the nsgt

function f = window_pad(g,pos,i,Ls)

f = zeros(Ls,1);

win_range= pos(i)+(-length(g{i})/2:length(g{i})/2-1);
if mod(length(g{i}),2)~=0
    win_range = win_range+.5;
end
f(win_range) = fftshift(g{i});

f = fftshift(ifft(f));