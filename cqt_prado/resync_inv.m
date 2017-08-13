function cqt_non_sync=resync_inv(cqt,Nfft,nbo,bins,R)
nl=size(cqt,1);
%nc=1+(size(cqt,2)-2^(nbo-1)-1);
decmin=Nfft/(2*R);
nc=size(cqt,2)-Nfft*(2^(nbo-2)/R);
cqt_non_sync=zeros(nl,nc);
for k=nbo:-1:1
    ind=(nbo-k)*bins;
    dec=Nfft*(2^(k-2)/R);
    if(k==1) && (decmin==0.5)
        dec=0;
    end
    cqt_non_sync(ind+1:ind+bins,:)=cqt(ind+1:ind+bins,1+dec:dec+nc);
end
