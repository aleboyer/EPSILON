function [k,P1,P2,P12]=get_cospectrum(a1,a2,k)
    Lscan=length(a1);
    dk=k(1);
    window = hanning(Lscan);
    wc2=1/mean(window.^2);            % window correction factor
    a1  = window.'.*(a1-mean(a1));
    a2  = window.'.*(a2-mean(a2));
    P1  = fft(a1);P2  = fft(a2);
    P12 = conj(P1).*P2./Lscan^2/dk*wc2;
    
    if rem(Lscan,2)==0
        k=-Lscan/2*dk:dk:Lscan/2*dk-dk;
    else
        kp=dk:dk:dk*floor(Lscan/2);
        k=[fliplr(-kp) 0 kp];
    end
