function y=cqt_inv(cqt,Nfft,nbo,bins,R,minfreq,b,fs_new,fs,type_fen)
cqt_non_sync=resync_inv(cqt,Nfft,nbo,bins,R);
nblocks = size(cqt_non_sync,2);
lensig=R*nblocks+Nfft-R;
freqmin=minfreq*(2^(nbo-1));
noyau_inv=noyauatrous_inv(freqmin,bins,fs_new,type_fen);
val=2^nbo;
x_oct=struct([]);
y_rec=struct([]);
for k=1:nbo
    x_oct(k).sig_rec=zeros(1,lensig);
end
y=[];
for nb_sub=1:nbo
    val=val/2;
    dec=R/val;
    yf=dec*noyau_inv*cqt_non_sync((nb_sub-1)*bins+1:nb_sub*bins,:);
    y_oct_temp = ifft(yf);
    y_oct = 2*real(y_oct_temp);
    lensig = Nfft + (nblocks-1)*dec;
    x_oct(nb_sub).sig = zeros(lensig,1);
    for n = 1:nblocks
        x_oct(nb_sub).sig((n-1)*dec+1:((n-1)*dec)+Nfft) = y_oct(:,n) + x_oct(nb_sub).sig((n-1)*dec+1:((n-1)*dec)+Nfft); %addition recouvrement
    end
    if(nb_sub~=1)
        x_oct(nb_sub).sig=x_oct(nb_sub).sig+y_rec(nb_sub-1).sig(1:length(y_rec(nb_sub-1).sig)-Nfft);
    end
    if(nb_sub~=nbo)
         y_rec(nb_sub).sig=upsample(x_oct(nb_sub).sig,2);
         y_rec(nb_sub).sig = 2*filtfilt(b,1,y_rec(nb_sub).sig);
    else
        y=x_oct(nb_sub).sig;
    end
end
kend=nextpow2(fs/fs_new);
for k=1:kend
    y=upsample(y,2);
    y=2*filtfilt(b,1,y);
end
