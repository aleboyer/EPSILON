if exist('/Volumes/aleboyer/ARNAUD/SCRIPPS/PLUMEX18/','dir')
    root_data='/Volumes/aleboyer/ARNAUD/SCRIPPS/PLUMEX18/';
    root_script='/Volumes/aleboyer/ARNAUD/SCRIPPS/EPSILON/';
else
    root_data='/Users/aleboyer/ARNAUD/SCRIPPS/PLUMEX18/';
    root_script='/Users/aleboyer/ARNAUD/SCRIPPS/EPSILON/';
end
Cruise_name='Plumex_Feb2018'; % 
WW_name='EPSIFISH'; % 
deployement='mar3';



WWpath=sprintf('%s/%s/%s/%s/L1/',root_data,Cruise_name,WW_name,deployement);
epsipath=sprintf('%s/%s/%s/%s/epsi/',root_data,Cruise_name,WW_name,deployement);


load([WWpath 'Turbulence_Profiles.mat'],'MS','EPSI_Profiles')
MS_5s=load([WWpath 'Turbulence_Profiles_5s.mat'],'MS','EPSI_Profiles');


i=1;j=100;
h_freq=get_filters_MADRE('MADRE2.1',MS{i}.f);


F=figure(1);
clf
loglog(MS{i}.f,squeeze(MS{i}.Pf(2,j,:)).*h_freq.FPO7(MS{i}.w(j)).','b','linewidth',2)
hold on
loglog(MS{i}.f,squeeze(MS{i}.Pf(2,j,:))./h_freq.electFPO7.','c')
loglog(MS{i}.f,squeeze(MS{i}.Pf(2,j,:))./h_freq.FPO7(MS{i}.w(j)).','m')

loglog(MS{i}.f,h_freq.electFPO7.','o-c')
loglog(MS{i}.f,h_freq.FPO7(MS{i}.w(j)).','o-m')
legend('FPO7(?C)','FPO7/sinc4','FPO7/(sinc4*probe response)','sinc4','probe response','location','southwest')
xlabel('Hz','fontsize',15)
ylabel('?C^2 /Hz','fontsize',15)
title('FPO7 2 - 3s spectra','fontsize',15)
grid on
xlim([MS{i}.f(1) MS{i}.f(end)])
set(gca,'fontsize',15)
F.PaperPosition=[0 0 10 5];
print('-dpng',[WWpath 'rawdeg.png'])



dataP=squeeze(MS_5s.MS{1}.PphiT_k(:,:,2));
mink=MS_5s.MS{1}.k(2);
maxk=MS_5s.MS{1}.k(find(~isnan(nansum(dataP,1)),1,'last'));
minP=min(dataP(:));
maxP=max(dataP(:));

F=figure(2);
clf
ax(1)=subplot(211);
loglog(MS{i}.k,squeeze(MS{i}.PphiT_k(100:120,:,2)),'linewidth',.5)
xlabel('cpm','fontsize',15)
ylabel('?C^2.m^{-2} /cpm','fontsize',15)
title('FPO7 2 - 3s spectra','fontsize',15)
xlim([mink maxk])
grid on
set(gca,'fontsize',15)
ax(2)=subplot(212);
loglog(MS_5s.MS{i}.k,squeeze(MS_5s.MS{i}.PphiT_k(50:60,:,2)),'linewidth',.5)
xlabel('cpm','fontsize',15)
ylabel('?C^2.m^{-2} /cpm','fontsize',15)
title('FPO7 2 - 5s spectra','fontsize',15)
grid on
xlim([mink maxk])
set(gca,'fontsize',15)
linkaxes(ax)

F.PaperPosition=[0 0 10 5];
print('-dpng',[WWpath 'Pdegk_3s_5s.png'])



Chi_class=calc_binned_chi(MS);
Chi_class5s=calc_binned_chi(MS_5s.MS);

indplot=[1 Nchi*Neps];
[Nchi,Neps,L]=size(Chi_class.Pbatch22);
dataP=reshape(Chi_class.mPphiT22,[Nchi*Neps,L]);
inddata=find(nansum(dataP,2)>0);
indplot1=max([indplot(1) 1]):min([indplot(end) length(inddata)]);
inddata=inddata(indplot1);
[chi,epsilon]=meshgrid(Chi_class.chibin,Chi_class.epsibin);
chi=chi(:);epsilon=epsilon(:);
databatch=reshape(Chi_class.Pbatch22,[Nchi*Neps,L]);
databatch(databatch==0)=nan;
minP=min(dataP(:));
maxP=max(dataP(:));

F3=figure(3);
cmap = parula(length(inddata)+1);
title_string='mar3'
ax(1)=subplot(211);
bb = loglog(Chi_class.kbatch,databatch(inddata,:),'color',[1 1 1]*.7,'linewidth',1);
hold on
b1 = loglog(Chi_class.kbatch,squeeze(Chi_class.Pbatch22(1,1,:)),'o-','color',[1 1 1]*.25,'linewidth',2);
b2 = loglog(Chi_class.kbatch,squeeze(Chi_class.Pbatch22(1,end,:)),'+-','color',[1 1 1]*.4,'linewidth',2);
b3 = loglog(Chi_class.kbatch,squeeze(Chi_class.Pbatch22(end,1,:)),'d-','color',[1 1 1]*.5,'linewidth',2);
b4 = loglog(Chi_class.kbatch,squeeze(Chi_class.Pbatch22(end,end,:)),'s-','color',[1 1 1]*.1,'linewidth',2);
for i=1:length(inddata)
    j=inddata(i);
    ll(i)=loglog(Chi_class.k,squeeze(dataP(j,:)),'Color',cmap(i,:),'linewidth',2);
end
xlim([mink maxk])
ylim([minP maxP])
grid on
xlabel('k [cpm]')
ylabel('[($\phi^T_k)^2$  / $^{\circ}C^2$ . $m^{-2}$ / cpm]','interpreter','latex')
set(gca,'fontsize',20)
title([title_string '- $\phi^{T22}_k$ 3second'],'interpreter','latex')



dataP=smoothdata(squeeze(MS_5s.MS{1}.PphiT_k(:,:,2)),1,'movmean',5);
mink=MS_5s.MS{1}.k(2);
maxk=MS_5s.MS{1}.k(find(~isnan(nansum(dataP,1)),1,'last'));
minP=min(dataP(dataP>0));
maxP=max(dataP(:));


ax(2)=subplot(212);
bb = loglog(Chi_class.kbatch,databatch(inddata,:),'color',[1 1 1]*.7,'linewidth',1);
hold on
b1 = loglog(Chi_class.kbatch,squeeze(Chi_class.Pbatch22(1,1,:)),'o-','color',[1 1 1]*.25,'linewidth',2);
b2 = loglog(Chi_class.kbatch,squeeze(Chi_class.Pbatch22(1,end,:)),'+-','color',[1 1 1]*.4,'linewidth',2);
b3 = loglog(Chi_class.kbatch,squeeze(Chi_class.Pbatch22(end,1,:)),'d-','color',[1 1 1]*.5,'linewidth',2);
b4 = loglog(Chi_class.kbatch,squeeze(Chi_class.Pbatch22(end,end,:)),'s-','color',[1 1 1]*.1,'linewidth',2);
i=1
loglog(MS_5s.MS{i}.k,dataP)
xlim([mink maxk])
ylim([minP maxP])
grid on
xlabel('k [cpm]')
ylabel('[($\phi^T_k)^2$  / $^{\circ}C^2$ . $m^{-2}$ / cpm]','interpreter','latex')
set(gca,'fontsize',20)
title([title_string '- $\phi^{T22}_k$ - 5second'],'interpreter','latex')

F3.PaperPosition=[0 0 15 10];
print('-dpng',[WWpath 'Pdegk_3s_5s_batch.png'])










