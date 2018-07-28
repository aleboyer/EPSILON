filename{1}='/Volumes/DataDrive/NISKINE/NISKINE/epsifish2/d3/L1/Turbulence_Profiles.mat';
filename{2}='/Volumes/DataDrive/NISKINE/NISKINE/epsifish1/d1/L1/Turbulence_Profiles.mat';
filename{3}='/Volumes/DataDrive/NISKINE/NISKINE/epsifish2/d4/L1/Turbulence_Profiles.mat';
filename{4}='/Volumes/DataDrive/NISKINE/NISKINE/epsifish1/d2/L1/Turbulence_Profiles.mat';
filename{5}='/Volumes/DataDrive/NISKINE/NISKINE/epsifish2/d5/L1/Turbulence_Profiles.mat';
filename{6}='/Volumes/DataDrive/NISKINE/NISKINE/epsifish1/d3/L1/Turbulence_Profiles.mat';
filename{7}='/Volumes/DataDrive/NISKINE/NISKINE/epsifish2/d6/L1/Turbulence_Profiles.mat';
filename{8}='/Volumes/DataDrive/NISKINE/NISKINE/epsifish1/d4/L1/Turbulence_Profiles.mat';
filename{9}='/Volumes/DataDrive/NISKINE/NISKINE/epsifish2/d8/L1/Turbulence_Profiles.mat';

zaxis=0:.5:1000;
tot_time=[];
tot_epsilon1=[];
tot_epsilon2=[];
tot_t=[];
tot_s=[];
tot_maxpr=[];

for f=1:length(filename)
    disp(filename{f})
    load(filename{f});
    MSempty=cellfun(@isempty,MS);
    Map_pr=cellfun(@(x) (x.pr),MS(~MSempty),'un',0);
    Map_epsilon1=cellfun(@(x) interp1(x.pr,x.epsilon(:,1),zaxis),MS(~MSempty),'un',0);
    Map_epsilon2=cellfun(@(x) interp1(x.pr,x.epsilon(:,2),zaxis),MS(~MSempty),'un',0);
    Map_t=cellfun(@(x) interp1(x.pr,x.t,zaxis),MS(~MSempty),'un',0);
    Map_s=cellfun(@(x) interp1(x.pr,x.s,zaxis),MS(~MSempty),'un',0);
    Map_time=cell2mat(cellfun(@(x) mean(x.time),MS(~MSempty),'un',0));
    Max_pr=cell2mat(cellfun(@(x) max(x.pr),MS(~MSempty),'un',0));
    Map_epsilon1=cell2mat(Map_epsilon1.');
    Map_epsilon2=cell2mat(Map_epsilon2.');
    Map_t=cell2mat(Map_t.');
    Map_s=cell2mat(Map_s.');
    tot_time=[tot_time Map_time];
    tot_maxpr=[tot_maxpr Max_pr];
    tot_epsilon1=[tot_epsilon1; Map_epsilon1];
    tot_epsilon2=[tot_epsilon2; Map_epsilon2];
    tot_t=[tot_t; Map_t];
    tot_s=[tot_s; Map_s];
end

tot_epsilon1(log10(tot_epsilon1)>-4)=nan;
tot_epsilon2(log10(tot_epsilon2)>-4)=nan;
sorttime=unique(sort(tot_time));
dt=nanmin(diff(sorttime));
timeaxis=min(tot_time)-dt:dt:max(tot_time)+dt;
epsilon1=zeros(length(zaxis),length(timeaxis)).*nan;
epsilon2=zeros(length(zaxis),length(timeaxis)).*nan;
temp=zeros(length(zaxis),length(timeaxis)).*nan;
sig =zeros(length(zaxis),length(timeaxis)).*nan;
for t=1:length(tot_time)
    ind=find(timeaxis>tot_time(t),1,'first');
    epsilon1(:,ind)=tot_epsilon1(t,:);
    epsilon2(:,ind)=tot_epsilon2(t,:);
    temp(:,ind)=tot_t(t,:);
    sig(:,ind)=tot_s(t,:);
end
smoothtemp=smoothdata(fillmissing(temp,'linear',2),'movmean',5);
%flagtemp=epsilon.*0;

figure;
colormap jet
ax(1)=subplot(211);
pcolor(timeaxis,zaxis,log10(epsilon1));shading flat;
hold on
contour(timeaxis,zaxis,smoothtemp,50,'k')
hold off
set(ax(1),'Xtick',timeaxis(1:20:end))
set(ax(1),'XtickLabel',datestr(timeaxis(1:20:end).'))
set(ax(1),'XTickLabelRotation',45)
axis ij
cax=colorbar;
caxis([-9.5 -7])
ylabel(cax,'W.kg^{-1}','fontsize',15)


colormap jet
ax(2)=subplot(212);
pcolor(timeaxis,zaxis,log10(epsilon2));shading flat;
hold on
contour(timeaxis,zaxis,smoothtemp,50,'k')
hold off
set(ax(2),'Xtick',timeaxis(1:20:end))
set(ax(2),'XtickLabel',datestr(timeaxis(1:20:end).'))
set(ax(2),'XTickLabelRotation',45)
axis ij
cax=colorbar;
caxis([-9.5 -7])
ylabel(cax,'W.kg^{-1}','fontsize',15)

linkaxes(ax)
