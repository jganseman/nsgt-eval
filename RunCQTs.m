% This script is related to the paper:
% Improving PLCA-based Score-Informed Source Separation 
% with Invertible Constant-Q Transforms
% J. Ganseman, P. Scheunders and S. Dixon
% EUSIPCO 2012

% Author: Joachim Ganseman, University of Antwerp, all rights reserved
% Adaptations and updates: december 2013

% This is a demo that demonstrates how the different implementations
% of (almost) invertible CQTs need to be run on data.
  

%%
clear;
addpath(genpath('./'));     % add subdirectories to path

createfigs=0;       % boolean: create figures for paper

%% Parameters

fftsize = 1024;
hopsize = fftsize / 4;

bins_per_octave = 48;       % make this a multiple of 12 for music


% Read file
datadir = getDataDirectory();       % directory with example files
[origMix,samplerate] = wavread([datadir 'pianoclip4notes.wav']);

  % if stereo, make mono
  if size(origMix, 2) > 1
     origMix = sum(origMix, 2) ./2 ; 
  end
  
  
%% PADDING: the exact amount needed depends on the largest kernel element and overlap parameters
% For now, just add a second before and after the signal. That should be
% enough to even take care of any windows of 11025 samples with 75% overlap.
origMix = [zeros(samplerate, 1); origMix; zeros(samplerate, 1)];
 

  
%% Test the Klapuri CQT transform

% set min and max frequencies available
fs=samplerate;
fmax = fs/3;     %center frequency of the highest frequency bin .
fmin = fmax/512; %lower boundary for CQT (lowest frequency bin will be immediately above this): fmax/<power of two> 
% at 44100, max=14700, min=28.7
% FYI, a piano goes from 27.5 to 4186 Hz. (we want some higher harmonics too though)

% prepare signal: cut out frequencies outside range
%x = [zeros(500,1); origMix; zeros(500,1)];
x=origMix;
w1 = 2*(fmin/(fs/2)); w2 = 0.8*(fmax/(fs/2));
[B,A] = butter(6,[w1 w2]); x = filtfilt(B,A,x); 

% Calculate CQT
KlapuriCQT = cqt_11(x,fmin,fmax,bins_per_octave,fs);     % more options available for hop size, tresholding, window function used
KlapuriReverted = icqt_11(KlapuriCQT);
%KlapuriReverted = KlapuriReverted(501:size(KlapuriReverted)-500);

% display!
plotCQT(KlapuriCQT,fs,0.6,'surf');

% interpolate, and display that one too
intCQT = getCQT(KlapuriCQT,'all','all');
figure;
spectrum_plot(intCQT); colormap(jet);


% calculate the perfect rasterized version
KlapuriRast = cqtPerfectRast(x,fmin,fmax,bins_per_octave,fs);
KlapuriRevertedRast = icqt_11(KlapuriRast);
%KlapuriRevertedRast = KlapuriRevertedRast(501:size(KlapuriRevertedRast)-500);

s = origMix;
s_r = KlapuriReverted;
rec_err = norm(s-s_r)/norm(s);
fprintf(['Klapuri general error:'...
    '   %e \n'],rec_err);

s = origMix;
s_r = KlapuriRevertedRast;
rec_err = norm(s-s_r)/norm(s);
fprintf(['Klapuri rasterized error:'...
    '   %e \n'],rec_err);



%% Test the Prado CQT transform

minfreq=fmin; % fréquence de début d'analyse
bins=bins_per_octave; % nombre de bins par octave
nbo=9;  % nombre d'octaves.
maxfreq=minfreq*(2^nbo);

assert(maxfreq == fmax);

%x = [zeros(500,1); origMix; zeros(500,1)];
x=origMix;

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
avance=0; %0.0005;           % MAKE LOW ENOUGH for reconstruction
% Filtre de décimation pour le calcul de la cqt
[b,bandwidth] = filtre_def();
[Pradocqt,t, f_cal, fs_new, noyau, Nfft, R,avance_reelle]=cqt_resamp(b,bandwidth,x,fs,minfreq,bins,nbo,avance);

if (R >= Nfft/2)
    disp('Step size too big w.r.t. FFT size.');
    disp('The reconstruction risks to be inaccurate.');
    disp('Either increase the nr of bins per octave,');
    disp('or diminish the stepsize.');
end
if avance_reelle==0
    disp('You are asking for too many octaves!');
    return;
end
txt=['Stepsize used was : ' num2str(avance_reelle) ' seconds'];
disp(txt);

%DISPLAY:
%figure;
%mesh(t,f_cal,abs(cqt));title('Prado CQT resample');
figure;
imagesc(abs(Pradocqt)); axis xy; colormap(jet);title('Prado CQT resample');

%RECONSTRUCT using several windows
% type_fen=1;
% y_rec_ham=cqt_inv(Pradocqt,Nfft,nbo,bins,R,minfreq,b,fs_new,fs,type_fen);
% y_rec_ham=y_rec_ham(501:length(origMix)+500);
type_fen=2;
y_rec_rect=cqt_inv(Pradocqt,Nfft,nbo,bins,R,minfreq,b,fs_new,fs,type_fen);
    %y_rec_rect=y_rec_rect(501:length(origMix)+500);
type_fen=3;
% y_rec_tuk=cqt_inv(Pradocqt,Nfft,nbo,bins,R,minfreq,b,fs_new,fs,type_fen);
% y_rec_tuk=y_rec_tuk(501:length(origMix)+500);

% I've added a couple of windows
% type_fen=4;
% y_rec_han=cqt_inv(Pradocqt,Nfft,nbo,bins,R,minfreq,b,fs_new,fs,type_fen);
% y_rec_han=y_rec_han(501:length(origMix)+500);
% type_fen=5;
% y_rec_blha=cqt_inv(Pradocqt,Nfft,nbo,bins,R,minfreq,b,fs_new,fs,type_fen);
% y_rec_blha=y_rec_blha(501:length(origMix)+500);


% Comparaison de reconstruction suivant la fenêtre utiliséée.
% figure;
% subplot(311);plot(x); hold on; plot(y_rec_ham,'r');hold off;
% axis([1 8192 -1.3 1.3]);hleg1 = legend('Original','Reconstruit');
% set(hleg1,'Location','NorthEast'); set(hleg1,'Interpreter','none');
% xlabel('Avec fenêtre de Hamming');
% subplot(312);plot(x); hold on; plot(y_rec_rect,'r');hold off;
% axis([1 8192 -1.3 1.3]);hleg1 = legend('Original','Reconstruit');
% set(hleg1,'Location','NorthEast'); set(hleg1,'Interpreter','none');
% xlabel('Avec fenêtre rectangulaire');
% subplot(313);plot(x); hold on; plot(y_rec_tuk,'r');hold off;
% axis([1 8192 -1.3 1.3]);hleg1 = legend('Original','Reconstruit');
% set(hleg1,'Location','NorthEast'); set(hleg1,'Interpreter','none');
% xlabel('Avec fenêtre de Tukey');

% s = origMix;
% s_r = y_rec_ham;
% rec_err = norm(s-s_r)/norm(s);
% fprintf(['Prado rast Hamming error:'...
%     '   %e \n'],rec_err);

s = origMix;
s_r = y_rec_rect(1:length(s));
rec_err = norm(s-s_r)/norm(s);
fprintf(['Prado rast Rectangular error:'...
    '   %e \n'],rec_err);
% 
% s = origMix;
% s_r = y_rec_tuk;
% rec_err = norm(s-s_r)/norm(s);
% fprintf(['Prado rast Tukey error:'...
%     '   %e \n'],rec_err);
% 
% s = origMix;
% s_r = y_rec_han;
% rec_err = norm(s-s_r)/norm(s);
% fprintf(['Prado rast Hann error:'...
%     '   %e \n'],rec_err);
% 
% s = origMix;
% s_r = y_rec_blha;
% rec_err = norm(s-s_r)/norm(s);
% fprintf(['Prado rast Blackman-Harris error:'...
%     '   %e \n'],rec_err);

% --> Seems like best results are gotten with the Rectangular win
% Really wondering: is the signal also windowed before the forward
% transform? Otherwise using a window for the backward transform seems
% pretty pointless.


%% Test the Constant-Q Nonstationary Gabor Transform

fmin = fmin; % Minimum desired frequency (in Hz)
fmax = fmax;              
bins = bins_per_octave; %bins = [12; 24; 36; 48; 12]; 

%[s,fs] = wavread('glockenspiel.wav'); name = 'Glockenspiel';
s=origMix; 
fs=samplerate;
Ls = length(s); % Length of signal (in samples)

% Window design
%  Define a set of windows for the nonstationary Gabor transform with
%  resolution evolving over frequency. In particular, the centers of the
%  windows correspond to geometrically spaced center frequencies.
[g,shift,M] = nsgcqwin(fmin,fmax,bins,fs,Ls);

% Compute corresponding dual windows.
gd = nsdual(g,shift,M);

% Plot the windows and the corresponding dual windows
figure;
subplot(211); plot_wins(g,shift, 1);        %1: normalize display
subplot(212); plot_wins(gd,shift, 1);

% Calculate the coefficients
c = nsgtf(s,g,shift,M);

% Plot the resulting spectrogram
figure;
plotnsgtf(c,shift,fs,2,60);

% Test reconstruction
s_r = nsigtf(c,gd,shift,Ls);

% Print relative error of reconstruction.
rec_err = norm(s-s_r)/norm(s);
fprintf(['NSGTF Relative reconstruction error:'...
    '   %e \n'],rec_err);



%% For a rasterized representation, re-set M (vector of hops per freq)
% to a vector of constant values, and use the value of the highest freq
maxwidth = M(size(M,1)/2);
rastc = nsgtf(s,g,shift,maxwidth);
figure;
plotnsgtf(rastc,shift,fs,2,60);


% Test reconstruction
gd = nsdual(g,shift,maxwidth);
s_r = nsigtf(rastc,gd,shift,Ls);
rec_err = norm(s-s_r)/norm(s);
fprintf(['Rasterized NSGTF Relative reconstruction error:'...
    '   %e \n'],rec_err);



%% Try the sliCQ representation instead:

% it calculates a multiple of 4 frames at a time, it seems. As the
% slice length, use the maximum hopsize (from previously)
maxhop = ceil(length(s)/M(2)/4)*4 ;
framesperslice = ceil(M(size(M,1)/2) / M(2)) /4 ;


sr= samplerate;
%slice length (samples):
sl_len = maxhop;       % put to same length as fftsize, 1024 or something
%transition area length (<= sl_len/2) (?)
tr_area = 16;       % can be kept very low
% desired number of time steps per slice
%M=framesperslice;        % if set to 0, this will compute a 1xN vector, 1 for every freq
            % its max value is 5464 for this particular signal
            % If set to a constant, a rectangular matrix is created.
[Sc,Sg,Sshift,SM,SLs,Ssl_len,Str_area] = slicq(s,fmin,fmax,bins,sl_len,tr_area,sr,framesperslice);
% note: sliCQ has different defaults than nsgtf alone! e.g. window is
% modified blackman harris! nsgcqwin function here called as:
% nsgcqwin(fmin,fmax,bins,sr,sl_len,'min_win',min_win,...
%        'Qvar',Qvar,'bwfac',4,'fractional',1,'winfun','modblackharr');

s_r = islicq(Sc,Sg,Sshift,SM,SLs,Ssl_len,Str_area);

rec_err = norm(s-s_r)/norm(s);
fprintf(['SliCQ Relative reconstruction error:'...
    '   %e \n'],rec_err);

%plot
plotslicq(Sc,Sshift,sr,fmin,fmax,bins,2,60)


%% Use own windowfunction, to rasterize

fmin = fmin; % Minimum desired frequency (in Hz)
fmax = fmax;              
bins = bins_per_octave; %bins = [12; 24; 36; 48; 12]; 

%[s,fs] = wavread('glockenspiel.wav'); name = 'Glockenspiel';
s=origMix; 
fs=samplerate;
Ls = length(s); % Length of signal (in samples)

% Window design
[gjoa,shiftjoa,Mjoa] = nsgfwin_joa(fmin,fmax,bins,fs,Ls);

% Compute corresponding dual windows.
gdjoa = nsdual(gjoa,shiftjoa,Mjoa);
% Note: can result in NaN values in the dual windows. Replace them by 0:
% EDIT: NaN resulted from the window generator having been changed: it now
% produces windows that are normalized, so no need for this line:
% for ii = 1:length(gd),gd{ii}(isnan(gd{ii})) = 0;end 


% Plot the windows and the corresponding dual windows
figure;
subplot(211); plot_wins(gjoa,shiftjoa);
subplot(212); plot_wins(gdjoa,shiftjoa);

% Calculate the coefficients
cjoa = nsgtf(s,gjoa,shiftjoa,Mjoa);

% Plot the resulting spectrogram
figure;
plotnsgtf(cjoa,shiftjoa,fs,2,60);

% remove Nyquist frequency band, replace by zeroes
maxbin = size(cjoa,1)/2;
cjoa{maxbin+1} = zeros(length(origMix), 1);

% Test reconstruction
s_r = nsigtf(cjoa,gdjoa,shiftjoa,Ls);

% Print relative error of reconstruction.
rec_err = norm(s-s_r)/norm(s);
fprintf(['NSGTF Rast-JG Relative reconstruction error:'...
    '   %e \n'],rec_err);




% 
% %% Create figure for paper
% 
% %%
% 
% myimage = flipud( abs(KlapuriRast.spCQT) );
% myimage = myimage(:,1500:size(myimage,2)-3000); 
% 
% 
% % perform plca
% 
% %[w h] = plca(myimage, 4, 50);
% 
% figure;
% imagesc(myimage); colormap(1-gray)
% 
% figure;
% 
% subplot( 2, 2, 2), imagesc( w*h), title('resynthesized CQT'), colormap(1-gray);
% %subplot( 2, 2, 3), stem( z), axis tight, title('components')
% subplot( 2, 2, 1), multiplot( fliplr(w')), view( [-90 90]), title('spectra')
% subplot( 2, 2, 4), multiplot( h), title('timelines')
% drawnow
% 
