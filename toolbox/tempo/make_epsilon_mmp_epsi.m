load([WWpath 'EpsiProfile.mat'],'Epsilon','Epsilon_mmp','Profiles')


fig=figure('Position',[100,100,800,1000]);
v = VideoWriter('../FIGURE/EPSILON/MakeEpsilon_EPSISPROUL_MMP_EPSI.avi');
v.FrameRate=10;
open(v)
i=1;
for j=1:length(Epsilon{i}.Pshear1)
    ref1(1,1)=max(Epsilon{i}.Pshear1{j}); ref1(2,1)=min(Epsilon{i}.Pshear1{j});
    ref2(1,1)=max(Epsilon_mmp{i}.Pshear1{j}); ref2(2,1)=min(Epsilon_mmp{i}.Pshear1{j});
    kcref1(1,1)=Epsilon{i}.kc1(j); kcref1(2,1)=Epsilon{i}.kc1(j);
    kcref2(1,1)=Epsilon_mmp{i}.kc1(j); kcref2(2,1)=Epsilon_mmp{i}.kc1(j);
    kmaxref(1,1)=Epsilon{i}.kmax(j); kmaxref(2,1)=Epsilon{i}.kmax(j);
    loglog(Epsilon{i}.k{j},Epsilon{i}.Pshear1{j},'b','linewidth',2);
    hold on
    loglog(Epsilon_mmp{i}.k{j},Epsilon_mmp{i}.Pshear2{j},'k','linewidth',2);
    loglog(kcref1,ref1,'b'); %loglog(kcref,ref,'x')
    loglog(kcref2,ref2,'k'); %loglog(kcref,ref,'x')
%    loglog(kmaxref,ref1,'g'); %loglog(kmaxref,ref,'+')
    loglog([50 50],ref1,'g'); %loglog(kmaxref,ref,'+')
    title(['Scan=' int2str(j) ', pr=' ...
        num2str(100*Epsilon{i}.pr(j)) ', epsilon EPSI=' num2str(Epsilon{i}.epsilon1(j)) ...
        ', epsilon MMP=' num2str(Epsilon_mmp{i}.epsilon1(j))]);
    xlabel('k (cpm)','fontsize',15)
    ylabel('\phi(k) / (s^{-2}/cpm)','fontsize',15)
    set(gca,'fontsize',15)
    loglog(Epsilon{i}.kpan1{j},Epsilon{i}.Ppan1{j},'bo','linewidth',2);
    loglog(Epsilon_mmp{i}.kpan1{j},Epsilon_mmp{i}.Ppan1{j},'ko','linewidth',2);
    legend('EPSI','MMP')
    set(gca,'fontsize',15)
    %ylim([0.1*min(Epsilon.Ppan1{j}) 1e3*max(Epsilon.Ppan1{j})])
    xlim([1e-1 1e3])
    ylim([1e-9 1e-1])
    hold off
    frame=getframe(fig);
    writeVideo(v,frame)
    clf;
end
close(v)
close all

                
