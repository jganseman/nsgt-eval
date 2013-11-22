function  [noyauatrous, Nfft]=noyauatrous_res(minfreq,bins,fs)
seuil = 0.0075; %seuil pour Hamming
Q=1/(2^(1/bins)-1);
Nfft=2^nextpow2(ceil(Q*fs/minfreq));
noyauatrous=[];
noyautmp=zeros(Nfft,1);
im=sqrt(-1);
for k=bins:-1:1;
    N=ceil(Q*fs/(minfreq*2^((k-1)/bins)));
    if rem(N,2)==1
        N=N-1;
    end
    n0=Nfft/2-N/2;
            noyautmp=[zeros(n0,1) ;...
                hamming(N,'periodic')/sum(hamming(N,'periodic')) ;...
                zeros(n0,1)];
    noyautmp=...
        noyautmp.*exp(2*pi*im*(minfreq*2^((k-1)/bins))*((0:Nfft-1)-n0)'/fs);
    specnoyau=fft(noyautmp);
    specnoyau(abs(specnoyau)<=seuil)=0;
    noyauatrous=sparse([specnoyau noyauatrous]);
end
noyauatrous=conj(noyauatrous);
