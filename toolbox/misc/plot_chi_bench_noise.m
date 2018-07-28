
SPROUL=load('/Volumes/aleboyer/ARNAUD/SCRIPPS/NUWC/READ/NOISE_ANALYSIS/SPROUL/d2/raw/SD_sproul_FPO7_noise_spectrum.mat')
Profile=load('/Volumes/aleboyer/ARNAUD/SCRIPPS/NUWC/READ/NOISE_ANALYSIS/SPROUL/d2/raw/newMAD_noISO_newMAP2_raw.mat');
% get and split the Profile
tscan=3;
Profile.nbsample=Profile.time;
[f1raw,P11raw]=calc_raw_epsi(Profile,tscan,f);

w=0.5;
G=34; % conversion V to Celsius
h_freq=get_filters_MADRE('MADRE2.1',f1raw);
TFtemp=h_freq.FPO7(w);

FPO7noisef_sinc4=34^2*squeeze(P11raw(1,100,:))./h_freq.electFPO7.';
FPO7noisef=34^2*squeeze(P11raw(1,100,:))./TFtemp.';
FPO7noisek_sinc4=FPO7noisef_sinc4./w;
FPO7noisek=FPO7noisef./w;
k=f1raw./w;

FPO7noisek_k_sinc4=(2*pi*k.').^2.*FPO7noisek_sinc4;
FPO7noisek_k=(2*pi*k.').^2.*FPO7noisek;


figure(3)
eps_bin=-10:.5:-6;
chi_bin=-10:.5:-6;
kvis=1.1508e-06;
ktemp= 1.4759e-07;
fig=figure;
hold on
for i=1:length(eps_bin)
    eps=10^(eps_bin(i));
    for j=1:length(chi_bin)
        chi=10^(chi_bin(j));
            [kbatch,Pbatch]=batchelor(eps,chi,kvis,ktemp);
            loglog(kbatch,Pbatch,'Color',[1 1 1]*.2)
    end
end

set(gca,'Xscale','log','Yscale','log')


hold on
loglog(k,FPO7noisek_k,'r')
loglog(k,FPO7noisek_k_sinc4,'g')

i=1;
for j=100:150
    figure(1)
    loglog(MS{i}.k,squeeze(MS{i}.PphiT_k(j,:,2)))
    hold on
    [kbatch,Pbatch]=batchelor(MS{i}.epsilon(j,2),MS{i}.chi(j,2),MS{i}.kvis(j),MS{i}.ktemp(j));
    loglog(kbatch,Pbatch)
    plot([1 1].*MS{i}.fc_index(j,2),[min(squeeze(MS{i}.PphiT_k(j,:,2))) max(squeeze(MS{i}.PphiT_k(j,:,2)))],'g')
    title(sprintf('Pr=%3.1f. log_{10} \\epsilon=%2.2f  log_{10} \\chi=%2.2f', MS{i}.pr(j),log10(MS{i}.epsilon(j,2)),log10(MS{i}.chi(j,2))))
    pause
    clf
end



