function [cqt,t,f_cal, fs_new, noyau, Nfft, R, avance_reelle]=cqt_resamp(b,bandwidth,x,fs,minfreq,bins,nbo,avance,fen_type)
% En entrée 
%   b : coefficients du filtre de décimation
%   bandwidth : pourcent de bande conservée
% (b, bandwidth et n_len sont renvoyés par la fonction filtre_def.m)
%   x : signal à traiter
%   fs : fréquence d'échantillonnage
%   minfreq : fréquence du début de la bande à traiter
%   bins : nombre de bins par octave
%   nbo : nombre d'octaves
%   avance : pas d'analyse en secondes
% En sortie
%   cqt : tableau des coefficients calculés
%   t : vecteur des instants de calcul
%   f_cal : vecteurs des fréquences correspondantes aux cqt
freqmin=minfreq*(2^(nbo-1));
freqmax=2*freqmin;
% Nombre de coef cqt à calculer
K=bins*nbo;
%Q=1/(2^(1/bins)-1);
y=x(:)'; %vecteur ligne
% Définition du sous échantillonnage possible sur le signal initial
%nk=floor(log2(bandwith*fs)-log2(minfreq))-nbo-1;
nk=floor(log2(bandwidth*fs)-log2(freqmax))-1;
if(nk<0)
    cqt.sig=0;
    t=0;
    f_cal=0;
    fs_new=0;
    noyau=0;
    Nfft=0;
    R=0;
    avance_reelle=0;
    return;
end
fs_new=fs/(2^nk);
% Calcul du noyau qui agit sur le dernier octave
[noyau, Nfft]=noyauatrous_res(freqmin,bins,fs_new);
% Avance en nombre d'échantillons
R=2^nextpow2(fs_new*avance);
% Vérifie que l'on peut faire R/2 -> R/(2^(nbo-1))
while R/(2^(nbo-1))<1
    R=2*R;
end
avance_reelle=R/fs_new;
% On fait le pré-traitement décimation avant calcul des CQT
% On applique nk fois un sous échantillonnage par 2
for nk1=1:nk
    y=filtfilt(b,1,y);
    y=downsample(y,2);
end;
% fin du pré-traitement
% On calcule les cqt
x_oct=struct([]);
x_oct(nbo).sig=y; % pour la derniere octave
for nb_sub=nbo-1:-1:1
    y=filtfilt(b,1,x_oct(nb_sub+1).sig);
    x_oct(nb_sub).sig=downsample(y,2);
end
%Recalage des signaux pour la synchronisation
cqt_non_sync=[];val=2^nbo;
for nb_sub=1:nbo
    val=val/2;
    xx=x_oct(nb_sub).sig;
    x_oct(nb_sub).sig=buffer(xx,Nfft,Nfft-R/val,'nodelay');
    for n_buf=1:size(x_oct(nb_sub).sig,2)
        cqt_non_sync((nb_sub-1)*bins+1:nb_sub*bins,n_buf)=constQ(x_oct(nb_sub).sig(:,n_buf)',noyau);
    end
end
%Resynchronise les coefficients cqt aux instants de calcul
%par rapport au signal
[cqt, t]=resync(cqt_non_sync,Nfft,nbo,bins,R,fs_new);
% Fin du calcul des cqt
%-------------------------------------------------------------------
f_cal=minfreq*2.^((0:K-1)/bins);