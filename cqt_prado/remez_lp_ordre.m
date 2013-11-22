function [lgfiltre, fo, ao, w] = remez_lp_ordre(f_caract, amplitude, ondulations, f_ech)
% Rabiner & Gold, Theory and Applications of Digital Signal Processing
%  fréquences normalisées
f_caract = f_caract/f_ech;       
a1_a3 = [5.309e-03 7.114e-02 -4.761e-01];
a4_a6 = [-2.660e-03 -5.941e-01 -4.278e-01];
d1 = log10(ondulations(1));
d2 = log10(ondulations(2));
D1 = [d1*d1 ; d1 ; 1];
D = (d2*a1_a3 + a4_a6)*D1;
fd1d2 =  11.01217 + 0.51244*(d1-d2);
df = f_caract(2) - f_caract(1);
% longueur du filtre
lgfiltre = ceil( D/df - fd1d2*df + 1 );
%=== Reformatage des données pour remez_lp.m
fo = [0 f_caract .5];
ao = amplitude;
w = ones(size(ondulations))*max(ondulations)./ondulations;
return;