% Dabbling with the STFT - a few experiments by Joachim Ganseman
% to better understand the nature of the (windowed) STFT and overlap-add
% reconstruction

% To get a good grasp of the STFT, I highly recommend reading
% http://www.katjaas.nl/FFTwindow/FFTwindow.html
% http://www.katjaas.nl/FFTwindow/FFTwindow&filtering.html
% Where some of the issues are illustrated. 
% This file was made to tests out a few things in Matlab

% first, let's create a signal. Any signal will do

sigsize = 4096;
signal = sin([1:1:sigsize]) + sin([1:1:sigsize]./4) + sin([1:1:sigsize]./20);
plot(signal)

allones = ones(1, 4096);        % for ease of understanding, take e.g. this
signal=allones;


% define an FFT length and hopsize
fftlen = 1024;
hopsize = fftlen / 2;

% DEFINE A WINDOW. To play with COLA condition checks, change hopsize and window here.
window = hann(fftlen, 'periodic');

% do an appropriate amount of zero padding
paddedsignal = [zeros(1, fftlen-hopsize) signal  zeros(1, fftlen-hopsize+mod(size(signal, 2),hopsize))];
figure
plot(paddedsignal)

% now chop this up in pieces and pack it into a matrix
% inspiring this code on the stft implementation by ellis

nrframes = length(paddedsignal) / hopsize - (fftlen/hopsize)+1;
matrix = zeros(fftlen, nrframes);
currentframe=1;

for i = 0:hopsize:(length(paddedsignal)-fftlen)
   matrix(:,currentframe) = paddedsignal(i+1:i+fftlen)';
   currentframe = currentframe+1;
end

% now do an overlap-add of this matrix

reconstructsignal = zeros(1, length(paddedsignal));
for i=1:nrframes
   index = (i-1)*hopsize;
   reconstructsignal(index+1:index+fftlen) = reconstructsignal(index+1:index+fftlen)+matrix(:,i)';
end



%% windowing and overlap-add


% apply to every column in this matrix
windowmatrix=matrix;
for i=1:nrframes
    windowmatrix(:,i) = windowmatrix(:,i) .*window;
end

% now suppose we do an FFT, and IFFT, resulting in the identical matrix.
% there's !!! no need to undo windowing (divide by the window) !!!
% at overlap-add, because the window is a partition of unity

reconstructsignal = zeros(1, length(paddedsignal));
for i=1:nrframes
   index = (i-1)*hopsize;
   reconstructsignal(index+1:index+fftlen) = reconstructsignal(index+1:index+fftlen)+windowmatrix(:,i)';
end

figure
plot(reconstructsignal)

% so in fact, at reconstruction, we can use the matrix as is 
% - that is, apply a rectangular synthesis window
% however, the picture changes when the STFT is changed and coefficients
% are modified:


%% windowing and modified STFTs

% let's take the fft of our windowmatrix
fftmatrix = fft(windowmatrix);
%plot this thing
figure
imagesc(20*log10(abs(windowmatrix)));

% now let's modify this signal for a bit. e.g., let's take out the 4th and
% 1021th component's real values:

fftmatrix(4,:) = zeros(1, nrframes) + imag(fftmatrix(4,:))*sqrt(-1);
fftmatrix(fftlen-3,:) = zeros(1, nrframes) + imag(fftmatrix(fftlen-3,:))*sqrt(-1);

% and run the IFFT
newmatrix = ifft(fftmatrix);

% now, when we overlap-add:
reconstructsignal = zeros(1, length(paddedsignal));
for i=1:nrframes
   index = (i-1)*hopsize;
   reconstructsignal(index+1:index+fftlen) = reconstructsignal(index+1:index+fftlen)+newmatrix(:,i)';
end
figure
plot(abs(reconstructsignal))

% a temporal 'ripple' was introduced because of the change of a component.
% this ripple is largest at the edges of the slices:
% (see also http://www.katjaas.nl/FFTwindow/FFTwindow&filtering.html)
% therefore, cutoffs appear at frame boundaries when overlap-adding.
% to undo this, smoothen all slices by windowing AGAIN:

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
plot(abs(reconstructsignal))

% at least, now our signal does not have these abrupt changes anymore.
% however, we've windowed twice, therefore with the square of the window,
% and for the Constant Overlap-Add property to hold, we need to have an
% overlap of at least 75% between frames (at least, for a Hann window)

%% Lesson Learned:

% if you're going to change coefficients of the STFT and then resynthesize,
% at least make sure you have 75% overlap, and a (Hann-)window applied both
% at synthesis AND resynthesis.

