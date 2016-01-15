% Chapter 2: creating the images for Chapter 2 of my PhD Thesis.
% (c) Joachim Ganseman, january 2014

% NOTE: AFTER RUNNING, EXTERNALIZE THE DATA TABLE READING in 2 steps:
% 1. add the following code before subimporting the created tex file:
%     \pgfplotstableread{./images/chapter2/audiosignalexcerpt-1.tsv}\loadedtable
% 2. in the tex file itself, change the following:
%     table[row sep=crcr,format=inline]{\loadedtable}; 

clear;

%%%% PART ONE: the DFT %%%%


PRINTTOFILE=0;

%%
addpath(genpath('./'));     % add subdirectories to path

% Read file
datadir = getDataDirectory();       % directory with example files
[origMix,samplerate] = wavread([datadir 'pianoclip4notes.wav']);

  % if stereo, make mono
  if size(origMix, 2) > 1
     origMix = sum(origMix, 2) ./2 ; 
  end

%%  First: a figure of an entire audio example
  
% plot file
figure;
plot(origMix);
xlabel('sample');
% however, the dataset is too big for TikZ. Save as image plot instead:

%enlarge the figure first:
%currentposition = get(gcf, 'position');
%newposition = currentposition .* [1 1 2 2];     % double width and height)
%set(gcf, 'position', newposition);

% set white background color
set(gcf, 'Color', 'white')

% save axes for later use
oldgca=get(gca);
% set axes off and get only the image data
set(gca,'position',[0 0 1 1],'units','normalized') %this maximizes the fig
axis off 
F=getframe(gcf);
colormap(F.colormap);
% create a new figure with the image data
figure;
imagesc(oldgca.XLim, oldgca.YLim, F.cdata);
set(gcf, 'Color', 'white');
xlabel('sample');

%filename='audiosignalfull';
%matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', [filename '.tex'])
%movefile(filename, '../../thesis/images/chapter2/')
%
%for i=1:37
%    curfile = [filename '-' num2str(i) '.tsv'];
%    movefile(curfile, '../../thesis/images/chapter2/');
%end
if PRINTTOFILE
filename='../../thesis/images/chapter2/audiosignalfull.tex';
matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end

%% Now: a figure of only a (1024-point?) part of an audio example

excerpt = origMix(2000:2000+1023);
figure;
plot(excerpt);
xlabel('sample');
set(gcf, 'Color', 'white')

if PRINTTOFILE
filename='../../thesis/images/chapter2/audiosignalexcerpt.tex';
matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.', 'floatFormat', '%.4g', 'encoding','UTF-8')
end

%% The DFT of such a signal

% just plot the 70 first components, the rest is near 0
excerptfft = fft(excerpt);
figure;
%first create a full plot, for axis boundaries
plot(real(excerptfft));
hold on;
plot(imag(excerptfft), 'r');
xlabel('frequency');
set(gcf, 'Color', 'white')
oldgca=get(gca);

%now plot only the 70 first components. This is where most of it happens
figure;
plot(real(excerptfft(1:70)));
hold on;
plot(imag(excerptfft(1:70)), 'r');
ylim(oldgca.YLim);
xlabel('frequency bin');
set(gcf, 'Color', 'white')

if PRINTTOFILE
filename='../../thesis/images/chapter2/fftbegin.tex';
matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end

%% to show symmetry, show the last 70 components too
% get hold of the axes first so the scale is the same
endpos=length(excerptfft)
figure;
plot([endpos-70:endpos], real(excerptfft(end-70:end)));
hold on;
plot([endpos-70:endpos], imag(excerptfft(end-70:end)), 'r');
ylim(oldgca.YLim);
xlabel('frequency bin');
set(gcf, 'Color', 'white')

if PRINTTOFILE
filename='../../thesis/images/chapter2/fftend.tex';
matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end

%% now plot the absolut value, to show the harmonics
figure;
plot(abs(excerptfft(1:100)));
xlabel('frequency bin');
set(gcf, 'Color', 'white')

if PRINTTOFILE
filename='../../thesis/images/chapter2/fftabs.tex';
matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end


%%
PRINTTOFILE=0

%%%% PART TWO: the STFT %%%%

signal = ones(1, 4096/4);
% define an FFT length and hopsize
fftlen = 1024/4;
hopsize = fftlen / 4;
% do an appropriate amount of zero padding
paddedsignal = [zeros(1, fftlen-hopsize) signal  zeros(1, fftlen-hopsize+mod(size(signal, 2),hopsize))];
figure;
plot(paddedsignal)
set(gcf, 'Color', 'white')
ylim([-0.5,1.5])
xlim([0, length(paddedsignal)])

if PRINTTOFILE
filename='../../thesis/images/chapter2/blocksignal.tex';
matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end

%%

% now chop this up in pieces and pack it into a matrix
% inspiring this code on the stft implementation by ellis
nrframes = length(paddedsignal) / hopsize - (fftlen/hopsize)+1;
matrix = zeros(fftlen, nrframes);
currentframe=1;
for i = 0:hopsize:(length(paddedsignal)-fftlen)
   matrix(:,currentframe) = paddedsignal(i+1:i+fftlen)';
   currentframe = currentframe+1;
end

%let's define a window for our matrix
window = hann(fftlen, 'periodic');
% apply to every column in this matrix
windowmatrix=matrix;
for i=1:nrframes
    windowmatrix(:,i) = windowmatrix(:,i) .*window;
end

% let's take the fft of our windowmatrix
fftmatrix = fft(windowmatrix);
%plot this thing
figure
imagesc(20*log10(abs(windowmatrix)));
set(gcf, 'Color', 'white');

if PRINTTOFILE
filename='../../thesis/images/chapter2/blockstft.tex';
matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end

%% now let's modify this signal for a bit. e.g., let's take out the 4th and
% 1021th component's real values:

fftmatrix(4,:) = zeros(1, nrframes);% + imag(fftmatrix(4,:))*sqrt(-1);
fftmatrix(fftlen-3,:) = zeros(1, nrframes);% + imag(fftmatrix(fftlen-3,:))*sqrt(-1);

% and run the IFFT
newmatrix = ifft(fftmatrix);

% now, when we overlap-add (perfect ISTFT:
reconstructsignal = zeros(1, length(paddedsignal));
for i=1:nrframes
   index = (i-1)*hopsize;
   reconstructsignal(index+1:index+fftlen) = reconstructsignal(index+1:index+fftlen)+newmatrix(:,i)';
end
figure
plot(abs(reconstructsignal)/2)
set(gcf, 'Color', 'white');
ylim([0.85,1.15])
xlim([500/4,1800/4])

if PRINTTOFILE
filename='../../thesis/images/chapter2/ISTFTperfect.tex';
matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end

%% smooth ISTFT
window2matrix=newmatrix;
for i=1:nrframes
    window2matrix(:,i) = window2matrix(:,i) .*window;
end

reconstructsignal = zeros(1, length(paddedsignal));
for i=1:nrframes
   index = (i-1)*hopsize;
   reconstructsignal(index+1:index+fftlen) = reconstructsignal(index+1:index+fftlen)+window2matrix(:,i)';
end
figure
plot(abs(reconstructsignal)/1.5)
set(gcf, 'Color', 'white');
ylim([0.85,1.15])
xlim([500/4,1800/4])

if PRINTTOFILE
filename='../../thesis/images/chapter2/ISTFTsmooth.tex';
matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end



%% Compute STFT of an entire signal

datadir = getDataDirectory();       % directory with example files
[origMix,samplerate] = wavread([datadir 'pianoclip4notes.wav']);

  % if stereo, make mono
  if size(origMix, 2) > 1
     origMix = sum(origMix, 2) ./2 ; 
  end
  
fftsize=1024;
hopsize=fftsize/4;
  
  % Do STFT
SmarSTFT = stft(origMix, fftsize, hopsize, 0, 'hann');
figure;
imagesc(20*log(abs(SmarSTFT)));
colormap(jet);
axis xy;
xlabel('signal slice')
ylabel('frequency bin')

if PRINTTOFILE
filename='../../thesis/images/chapter2/stftplot.tex';
matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end


%% Compute CQT, Schorkhuber/Klapuri version
% set sample rate and minimum/maximum frequencies
fs=samplerate;
fmax = fs/3;     %center frequency of the highest frequency bin .
fmin = fmax/512;
bins_per_octave = 48; 

x=origMix;
w1 = 2*(fmin/(fs/2)); w2 = 0.8*(fmax/(fs/2));
[B,A] = butter(6,[w1 w2]); x = filtfilt(B,A,x); 

KlapuriCQT = cqt(x,fmin,fmax,bins_per_octave,fs); 
plotCQT(KlapuriCQT,fs,0.25,'surf');        % use 'image' to display only coordinates
xlabel('time [sec]')
ylabel('frequency [Hz]')

% this gets saved as data. Instead, display an image
oldgca=get(gca);
% set axes off and get only the image data
set(gca,'position',[0 0 1 1],'units','normalized') %this maximizes the fig
axis off 
F=getframe(gcf);
colormap(F.colormap);
% create a new figure with the image data
figure;
imagesc([-0.5, 3.85], oldgca.YLim, F.cdata);     % hard set limits here
set(gcf, 'Color', 'white');
xlabel('time (sec)')
ylabel('frequency')
%reset Y axis labeling, inspired on PlotCQT from Schorkhuber/Klapuri CQT
    % did it times 2
set(gca,'YTick',1:KlapuriCQT.bins:KlapuriCQT.octaveNr*KlapuriCQT.bins);
h = get(gca); yTick = flipud(h.YTick');     % note: flipped axis here!
yTickLabel = num2str(round(KlapuriCQT.fmin*2.^((yTick-1)/KlapuriCQT.bins)),5);
set(gca,'YTickLabel',yTickLabel);

    
if PRINTTOFILE
filename='../../thesis/images/chapter2/cqtplot.tex';
matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end

%% compute NSGT
s=origMix; 
fs=samplerate;
Ls = length(s);
bins=bins_per_octave;

[g,shift,M] = nsgcqwin(fmin,fmax,bins,fs,Ls);
% Compute corresponding dual windows.
gd = nsdual(g,shift,M);

% plot the windows
figure;
plot_wins(g,shift, 1);        %1: normalize
ylim([-0.2, 1.2]);
xlim([0, Ls/2]);

c = nsgtf(s,g,shift,M);
figure;
plotnsgtf(c,shift,fs,fmin, fmax, bins, 2,60);
set(gcf, 'Color', 'white');
xlabel('time [sec]')
ylabel('frequency [Hz]')

% This is displayed as an entire series of images in plotnsgtf.
% Convert this to a single image
oldgca=get(gca);
% set axes off and get only the image data
set(gca,'position',[0 0 1 1],'units','normalized') %this maximizes the fig
axis off 
F=getframe(gcf);
colormap(F.colormap);
% create a new figure with the image data
figure;
imagesc(oldgca.XLim, oldgca.YLim, F.cdata);
set(gcf, 'Color', 'white');
xlabel('time [sec]')
ylabel('frequency [Hz]')

%reset Y axis labeling, inspired on plotnsgtf.m
    N = length(shift);
    L=sum(shift);
    vfq = [0, fmin*2.^(0:log2(fmax/fmin))];
    cutout=2;   % second-to-last parameter in plotnsgtf
    if length(bins) < length(vfq) - 1
        xbins = [bins, bins(end)*ones(1,length(vfq)-length(bins)-1)];
    else
        xbins = bins(1:length(vfq)-1);
    end;
    sbins = [0,2,1+cumsum(xbins)];
    yTick = [2;2+2*cumsum(xbins(1:end-1))';(1+N/2)];
    yTick = unique([yTick(yTick <= floor((length(shift)-1)/cutout)+1);...
        floor((N-1)/cutout)+1]);
    yTickLabel = zeros(length(yTick),1);
    for kk = 1:length(yTick)
        ind = find(yTick(kk) < sbins,1);
        yTickLabel(kk) = vfq(ind-1)*...
            2^(((yTick(kk)-1)-sbins(ind-1))/(sbins(ind)-sbins(ind-1)));
    end;
    yTickLabel = num2str(flipud(round(yTickLabel)),5);      % added flipud
    set(gca,'YTick',yTick,'YTickLabel',yTickLabel);

    
if PRINTTOFILE
filename='../../thesis/images/chapter2/nsgtplot.tex';
matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end

%% My own Window design for oversampled rectangular matrix
% pad signal so it better corresponds to previous images
s = [zeros(fs/2,1); origMix; zeros(fs/2,1)];
Ls= length(s);
[g,shift,M] = nsgfwin_joa(fmin,fmax,bins,fs,Ls);
% Calculate the coefficients
c = nsgtf(s,g,shift,M);
figure;
plotnsgtf(c,shift,fs,fmin, fmax, bins, 2,60);
xlabel('time [sec]')
ylabel('frequency [Hz]')

% This is displayed as an entire series of images in plotnsgtf.
% Convert this to a single image
oldgca=get(gca);
% set axes off and get only the image data
set(gca,'position',[0 0 1 1],'units','normalized') %this maximizes the fig
axis off 
F=getframe(gcf);
colormap(F.colormap);
% create a new figure with the image data
figure;
imagesc(oldgca.XLim, oldgca.YLim, F.cdata);    % +[-0.5,0.5]
set(gcf, 'Color', 'white');
xlabel('time [sec]')
ylabel('frequency [Hz]')

%reset Y axis labeling, inspired on plotnsgtf.m
    N = length(shift);
    L=sum(shift);
    vfq = [0, fmin*2.^(0:log2(fmax/fmin))];
    cutout=2;   % second-to-last parameter in plotnsgtf
    if length(bins) < length(vfq) - 1
        xbins = [bins, bins(end)*ones(1,length(vfq)-length(bins)-1)];
    else
        xbins = bins(1:length(vfq)-1);
    end;
    sbins = [0,2,1+cumsum(xbins)];
    yTick = [2;2+2*cumsum(xbins(1:end-1))';(1+N/2)];
    yTick = unique([yTick(yTick <= floor((length(shift)-1)/cutout)+1);...
        floor((N-1)/cutout)+1]);
    yTickLabel = zeros(length(yTick),1);
    for kk = 1:length(yTick)
        ind = find(yTick(kk) < sbins,1);
        yTickLabel(kk) = vfq(ind-1)*...
            2^(((yTick(kk)-1)-sbins(ind-1))/(sbins(ind)-sbins(ind-1)));
    end;
    yTickLabel = num2str(flipud(round(yTickLabel)),5);      % added flipud
    set(gca,'YTick',yTick,'YTickLabel',yTickLabel);

if PRINTTOFILE
filename='../../thesis/images/chapter2/nsgtresampplot.tex';
matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end