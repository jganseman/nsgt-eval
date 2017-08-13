function cq=constQ(x, noyauatrous)
% x est un vecteur ligne
cq=fft(x, size(noyauatrous,1))*noyauatrous/size(noyauatrous,1);