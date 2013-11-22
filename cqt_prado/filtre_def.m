function [b, bandwidth] = filtre_def()
% Calcul du filtre RIF de décimation par 2.
rp = 0.1;          % Passband ripple
rs = 120;          % Stopband ripple
bandwidth=0.95;     % Pourcent de la largeur de bande conservée
f = [bandwidth*0.25 0.25];  % Cutoff frequencies
a = [1 0];        % Desired amplitudes
% Compute deviations
dev = [(10^(rp/20)-1)/(10^(rp/20)+1)  10^(-rs/20)]; 
[n,fo,ao,w] = remez_lp_ordre(f,a,dev,1);
% % filtre de longueur impaire nécessaire pour un retard entier
% if rem(n,2)~=1;n=n+1;end;
b = remez_lp(n,fo,ao,w);
% len_b=length(b);
% n_len=2^(nextpow2(len_b))+1;
% z_pad=(n_len-len_b)/2;
% b=[zeros(1,z_pad) b];