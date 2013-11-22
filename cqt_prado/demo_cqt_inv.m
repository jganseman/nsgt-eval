[xx,fs]=wavread('test_small_mono.wav');
x=xx(:,1)';
%_____________________
minfreq=55; % fréquence de début d'analyse
bins=72; % nombre de bins par octave
nbo=8;  % nombre d'octaves.
maxfreq=minfreq*(2^nbo);
% On va filtrer pour comparer les mêmes bandes de fréquences.
rp = 3;           % Passband ripple
rs = 70;          % Stopband ripple
f = [maxfreq min(1.1*maxfreq,fs/2)];    % Cutoff frequencies
a = [1 0];        % Desired amplitudes
dev = [(10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)]; 
[n,fo,ao,w] = remez_lp_ordre(f,a,dev,fs);
bdp = remez_lp(n,fo,ao,w);
x=filtfilt(bdp,1,x);
% pas d'analyse en seconde demandé
avance=0.005;
% Filtre de décimation pour le calcul de la cqt
[b,bandwidth] = filtre_def();
[cqt,t, f_cal, fs_new, noyau, Nfft, R,avance_reelle]=cqt_resamp(b,bandwidth,x,fs,minfreq,bins,nbo,avance);
if (R >= Nfft/2)
    disp('Le pas d''avance est trop grand par rapport à la taille de la FFT.');
    disp('La reconstruction risque de ne pas être correcte.');
    disp('Vous devriez augmenter le nombre de bins par octave,');
    disp('ou diminuer le pas d''avance.');
end
if avance_reelle==0
    disp('Trop d''octaves demandées!');
    return;
end
txt=['L''avance utilisée est : ' num2str(avance_reelle) ' secondes'];
disp(txt);
figure(1);
mesh(t,f_cal,abs(cqt));title('CQT resample');
% Reconstruction du signal
% type_fen = 1 : fenêtre de Hamming
% type_fen = 2 : fenêtre rectangulaire
% type_fen = 3 : fenêtre de Tukey
type_fen=1;
y_rec_ham=cqt_inv(cqt,Nfft,nbo,bins,R,minfreq,b,fs_new,fs,type_fen);
type_fen=2;
y_rec_rect=cqt_inv(cqt,Nfft,nbo,bins,R,minfreq,b,fs_new,fs,type_fen);
type_fen=3;
y_rec_tuk=cqt_inv(cqt,Nfft,nbo,bins,R,minfreq,b,fs_new,fs,type_fen);
% Comparaison de reconstruction suivant la fenêtre utiliséée.
figure(2);
subplot(311);plot(x(124001:132192)); hold on; plot(y_rec_ham(124001:132192),'r');hold off;
axis([1 8192 -1.3 1.3]);hleg1 = legend('Original','Reconstruit');
set(hleg1,'Location','NorthEast'); set(hleg1,'Interpreter','none');
xlabel('Avec fenêtre de Hamming');
subplot(312);plot(x(124001:132192)); hold on; plot(y_rec_rect(124001:132192),'r');hold off;
axis([1 8192 -1.3 1.3]);hleg1 = legend('Original','Reconstruit');
set(hleg1,'Location','NorthEast'); set(hleg1,'Interpreter','none');
xlabel('Avec fenêtre rectangulaire');
subplot(313);plot(x(124001:132192)); hold on; plot(y_rec_tuk(124001:132192),'r');hold off;
axis([1 8192 -1.3 1.3]);hleg1 = legend('Original','Reconstruit');
set(hleg1,'Location','NorthEast'); set(hleg1,'Interpreter','none');
xlabel('Avec fenêtre de Tukey');



