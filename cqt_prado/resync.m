function [cqt_sync, t]=resync(cqt,Nfft,nbo,bins,R,fs_new)
nl=size(cqt,1);nc=size(cqt,2);
%nc_ns=2^(nbo-1)+(nc-1)+1;
decmin=Nfft/(2*R);
nc_ns=nc+Nfft*(2^(nbo-2)/R);
cqt_sync=zeros(nl,nc_ns);
for k=nbo:-1:1
    ind=(nbo-k)*bins;
    dec=Nfft*(2^(k-2)/R);
    if(k==1) && (decmin==0.5)
        dec=0;
    end
    cqt_sync(ind+1:ind+bins,1+dec:dec+nc)=cqt(ind+1:ind+bins,:);
end
t=((Nfft/2-decmin)+R*(0:nc_ns-1))/fs_new;
