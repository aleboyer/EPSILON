% plot shear spectra and panchev as a function of epsilon level
root_data='/Users/aleboyer/ARNAUD/SCRIPPS/SHEAR_PROBE/';
root_script='/Users/aleboyer/ARNAUD/SCRIPPS/WireWalker/WW_process_profile/';
Cruise_name='EPSI_SPROUL'; % 
WW_name='EPSI'; % 
deployement='d2';
%% add the needed toobox 
%addpath /Users/aleboyer/ARNAUD/SCRIPPS/WireWalker/scripts/mixing_library/mixing_library/private1/seawater
addpath toolbox/
addpath toolbox/CTD
addpath toolbox/FILTER
addpath toolbox/process
%% define path
WWpath=sprintf('%s/%s/%s/%s/L1/',root_data,Cruise_name,WW_name,deployement);
epsipath=sprintf('%s/%s/%s/%s/epsi/',root_data,Cruise_name,WW_name,deployement);
name_rbr=[WW_name '_ctd_' deployement];

load([WWpath 'EpsiProfile.mat'],'Epsilon','chi','Profiles')

maxk=0;mink=1e10;
%for j=1:length(Epsilon)
clear Epsi;clear J;clear JJ;
for j=1:length(Epsilon)
    Epsi{j}   = Epsilon{j}.epsilon1;
    J{j}      = j+Epsilon{j}.epsilon1*0;
    JJ{j}     = 1:length(Epsilon{j}.epsilon1);
    A         = [Epsilon{j}.k{:}];
    maxk      = max([maxk,nanmax(A(:))]);
    mink      = min([mink,nanmin(A(:))]);
end

Epsi1 = [Epsi{:}];
J1    = [J{:}];
JJ1   = [JJ{:}];

[sEpsi1,I]=sort(Epsi1);
sJ1=J1(I);
sJJ1=JJ1(I);

maxEpsi=sEpsi1(end);
meanEpsi=nanmean(sEpsi1);

ind_high = find(log10(sEpsi1)>=-7.5 & log10(sEpsi1)<=-7);
ind_mid  = find(log10(sEpsi1)>=-9 & log10(sEpsi1)<=-8.5);
ind_low  = find(log10(sEpsi1)<=-9.5);

kaxis=mink:nanmean(diff(Epsilon{1}.k{1})):maxk;

H=length(ind_high);
M=length(ind_mid);
L=length(ind_low);

Pshear_H = zeros(H,length(kaxis));
Ppan_H   = zeros(H,length(kaxis));
test   = zeros(1,length(kaxis));
for i=1:length(ind_high)
    A  = Epsilon{sJ1(ind_high(i))}.Pshear1{sJJ1(ind_high(i))};
    B  = Epsilon{sJ1(ind_high(i))}.k{sJJ1(ind_high(i))};
    C  = Epsilon{sJ1(ind_high(i))}.Ppan1{sJJ1(ind_high(i))};
    D  = Epsilon{sJ1(ind_high(i))}.kpan1{sJJ1(ind_high(i))};
    test(i) = Epsilon{sJ1(ind_high(i))}.epsilon1(sJJ1(ind_high(i)));
    Pshear_H(i,:) = interp1(B,A,kaxis);
    Ppan_H(i,:)   = interp1(D,C,kaxis);
    
end
[kpan,Ppan_H]=panchev(10^(-7.25),10^(-6),kaxis);

Pshear_M = zeros(M,length(kaxis));
Ppan_M   = zeros(M,length(kaxis));
test   = zeros(1,length(kaxis));
for i=1:length(ind_mid)
    A  = Epsilon{sJ1(ind_mid(i))}.Pshear1{sJJ1(ind_mid(i))};
    B  = Epsilon{sJ1(ind_mid(i))}.k{sJJ1(ind_mid(i))};
    C  = Epsilon{sJ1(ind_mid(i))}.Ppan1{sJJ1(ind_mid(i))};
    D  = Epsilon{sJ1(ind_mid(i))}.kpan1{sJJ1(ind_mid(i))};
    test(i) = Epsilon{sJ1(ind_mid(i))}.epsilon1(sJJ1(ind_mid(i)));
    Pshear_M(i,:) = interp1(B,A,kaxis);
    Ppan_M(i,:)   = interp1(D,C,kaxis);
    
end
[kpan,Ppan_M]=panchev(10^(-8.75),10^(-6),kaxis);



Pshear_L = zeros(L,length(kaxis));
Ppan_L   = zeros(L,length(kaxis));
test   = zeros(1,length(kaxis));
for i=1:length(ind_low)
    A  = Epsilon{sJ1(ind_low(i))}.Pshear1{sJJ1(ind_low(i))};
    B  = Epsilon{sJ1(ind_low(i))}.k{sJJ1(ind_low(i))};
    C  = Epsilon{sJ1(ind_low(i))}.Ppan1{sJJ1(ind_low(i))};
    D  = Epsilon{sJ1(ind_low(i))}.kpan1{sJJ1(ind_low(i))};
    test(i) = Epsilon{sJ1(ind_low(i))}.epsilon1(sJJ1(ind_low(i)));
    Pshear_L(i,:) = interp1(B,A,kaxis);
    Ppan_L(i,:)   = interp1(D,C,kaxis);
    
end

[kpan,Ppan_L]=panchev(10^(-9.75),10^(-6),kaxis);

maxPshear_H=max(Pshear_H);
minPshear_H=min(Pshear_H);

SEM = std(Pshear_H)/sqrt(size(Pshear_H,2));     % Standard Error
ts = tinv([0.025;0.975],size(Pshear_H,2)-1);    % T-Score
CI_H = nanmean(Pshear_H) + ts.*SEM;             % Confidence Intervals
SEM = std(Pshear_M)/sqrt(size(Pshear_M,2));     % Standard Error
ts = tinv([0.025;0.975],size(Pshear_M,2)-1);    % T-Score
CI_M = nanmean(Pshear_M) + ts.*SEM;             % Confidence Intervals
SEM = std(Pshear_L)/sqrt(size(Pshear_L,2));     % Standard Error
ts = tinv([0.025;0.975],size(Pshear_L,2)-1);    % T-Score
CI_L = nanmean(Pshear_L) + ts.*SEM;             % Confidence Intervals



close all
loglog(kaxis,mean(Pshear_H,1),'r','linewidth',2)
hold on
loglog(kaxis,mean(Pshear_M,1),'g','linewidth',2)
loglog(kaxis,mean(Pshear_L,1),'b','linewidth',2)

loglog(kaxis,CI_H,'r--')
loglog(kaxis,CI_M,'g--')
loglog(kaxis,CI_L,'b--')

% loglog(kaxis,mean(Ppan_H,1),'ro','linewidth',2)
% loglog(kaxis,mean(Ppan_M,1),'go','linewidth',2)
% loglog(kaxis,mean(Ppan_L,1),'bo','linewidth',2)
loglog(kpan,Ppan_H,'ro','linewidth',2)
loglog(kpan,Ppan_M,'go','linewidth',2)
loglog(kpan,Ppan_L,'bo','linewidth',2)


hold off
xlabel('k (cpm)','fontsize',20)
ylabel('\phi_{du/dz} (s^{-2}/cpm)','fontsize',20)
set(gca,'fontsize',15)
legend('-7.5<log_{10}(\epsilon)<-7',...
       '-9<log_{10}(\epsilon)<-8.5',...
       'log_{10}(\epsilon)<-9.5')

% ind_high = find(log10(sEpsi1)>=-7.5 & log10(sEpsi1)<=-7);
% ind_mid  = find(log10(sEpsi1)>=-9 & log10(sEpsi1)<=-8.5);
% ind_low  = find(log10(sEpsi1)<=-9.5);

   
   
   
xlim([1,1e3]);
print('-dpng2','../FIGURE/EPSILON/HML_range_EPSISPROUL_EPSISBE.png')

return
Map_pr=cellfun(@(x) (x.pr),Epsilon,'un',0);
zaxis=min([Map_pr{:}]):.5:max([Map_pr{:}]);
Map_epsilon=cellfun(@(x) interp1(x.pr,x.epsilon1,zaxis),Epsilon,'un',0);
Map_epsilon=cell2mat(Map_epsilon.');
RBRgrid.Epsilon=interp1(zaxis,Map_epsilon.',RBRgrid.z);


level=min(RBRgrid.rho(:)):.01:max(RBRgrid.rho(:));
L=length(level);

close all
figure;
colormap('jet')
pcolor(RBRgrid.time,RBRgrid.z,log10(RBRgrid.Epsilon));shading interp;axis ij
hold on
contour(RBRgrid.time,RBRgrid.z,RBRgrid.rho,[level(1:5:L-200) level(L-200:L)],'k')
hold off
caxis([-9,-6])
set(gca,'XTickLabelRotation',45)
set(gca,'XTick',RBRgrid.time(1:5:end))
datetick
xlim([RBRgrid.time(1) RBRgrid.time(110)])
cax=colorbar;
xlabel(['Start date :' datestr(RBRgrid.time(1),'mm-dd-yyyy')],'fontsize',15)
set(gca,'fontsize',15)
ylabel(cax,'log_{10}(\epsilon)','fontsize',20)
ylabel('Depth (m)','fontsize',20)
%xlim([RBRgrid.time(1) RBRgrid.time(end)])
ylim([2 48])



















