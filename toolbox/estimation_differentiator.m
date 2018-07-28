NOISE=load('/Users/aleboyer/ARNAUD/SCRIPPS/EPSILOMETER/NOISE_ANALYSIS/SPROUL/d2/raw/comparison_temp_granite_sproul.mat');
FPO7diff=load('toolbox/FILTER/FPO7_coeffilt.mat','coef_filt','freq');

WWpath='/Volumes/Ahua/data_archive/WaveChasers-DataArchive/EPSI_SPROUL//EPSI_SPROUL/EPSI/d2/L1/';
load([WWpath 'Turbulence_Profiles.mat'],'MS','EPSI_Profiles')


Chi_class=calc_binned_chi(MS);

[F1,F2,legend_string]=plot_binned_chi(Chi_class,'SPROUL',10:20);


H=get_filters_MADRE('MADRE2.1',NOISE.k_granite);
diffTF=interp1(FPO7diff.freq,FPO7diff.coef_filt,NOISE.k_granite);
noise=NOISE.spec_granite./(H.FPO7(.6) .* diffTF.^2);
noise=interp1(NOISE.k_granite/.6,noise,Chi_class.k);


hold on
lnoise1=loglog(Chi_class.k,(2*pi*Chi_class.k).^2 .*.8^2.*noise/.6,'g');
legend_string{length(legend_string)+1}='differentiator noise';

print(F1,'Differentiator_noise_estimation.png','-dpng2')


