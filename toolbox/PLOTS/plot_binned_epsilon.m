function [F1,F2]=plot_binned_epsilon(Epsilon_class,title_string)

F1 = figure(1);clf
set(F1,'position',[100 50 1300 700])
axes('position',[.1 .1 .55 .82])
cmap = parula(length(Epsilon_class.bin)+1);

pp = loglog(Epsilon_class.kpan.',Epsilon_class.Ppan.','color',[1 1 1]*.7,'linewidth',1);
ll=pp;

mink=Epsilon_class.k(find(~isnan(nansum(Epsilon_class.mPshear1,1)),1,'first'));
maxk=Epsilon_class.k(find(~isnan(nansum(Epsilon_class.mPshear1,1)),1,'last'));
minP=min(Epsilon_class.mPshear1(:));
maxP=max(Epsilon_class.mPshear1(:));

hold on
for i=1:length(Epsilon_class.bin)
    indk=find(~isnan(Epsilon_class.Ppan(i,:)),1,'first');
    if Epsilon_class.kpan(i,indk)<=mink
        text(mink, ...
            Epsilon_class.Ppan(i,indk),...
            sprintf('10^{%1.1f}',log10(Epsilon_class.bin(i))),...
            'fontsize',20)
    else
        text(Epsilon_class.kpan(i,indk), ...
            Epsilon_class.Ppan(i,indk),...
            sprintf('10^{%1.1f}',log10(Epsilon_class.bin(i))),...
            'fontsize',20)
    end
    ll(i)=loglog(Epsilon_class.k,Epsilon_class.mPshear1(i,:),'Color',cmap(i,:),'linewidth',2);
end

xlim([mink maxk])
ylim([minP maxP])
grid on
xlabel('k [cpm]')
ylabel('[s^{-2} / cpm]')
set(gca,'fontsize',20)
title([title_string '- Shear1'])

nanind=find(nansum(Epsilon_class.mPshear1,2)>0);
legend_string=arrayfun(@(x) sprintf('log_{10}(\\epsilon)~%2.1f',x),log10(Epsilon_class.bin),'un',0);
hl=gridLegend(ll(nanind),2,legend_string{nanind},'location','northwest');
set(hl,'fontsize',10)

axes('position',[.7 .1 .25 .82])

plot(log10(Epsilon_class.bin),Epsilon_class.nbin1,'-+k')
set(gca,'fontsize',20)
xlabel('log_{10} \epsilon')
xlim([min(log10(Epsilon_class.bin)) max(log10(Epsilon_class.bin))])
title('PDF \epsilon_1')



F2 = figure(2);clf
set(F2,'position',[100 50 1300 700])
axes('position',[.1 .1 .55 .82])
cmap = parula(length(Epsilon_class.bin)+1);

pp = loglog(Epsilon_class.kpan.',Epsilon_class.Ppan.','color',[1 1 1]*.7,'linewidth',1);
ll=pp;

mink=Epsilon_class.k(find(~isnan(nansum(Epsilon_class.mPshear2,1)),1,'first'));
maxk=Epsilon_class.k(find(~isnan(nansum(Epsilon_class.mPshear2,1)),1,'last'));
minP=min(Epsilon_class.mPshear2(:));
maxP=max(Epsilon_class.mPshear2(:));

hold on
for i=1:length(Epsilon_class.bin)
    indk=find(~isnan(Epsilon_class.Ppan(i,:)),1,'first');
    if Epsilon_class.kpan(i,indk)<=mink
        text(mink, ...
            Epsilon_class.Ppan(i,indk),...
            sprintf('10^{%1.1f}',log10(Epsilon_class.bin(i))),...
            'fontsize',20)
    else
        text(Epsilon_class.kpan(i,indk), ...
            Epsilon_class.Ppan(i,indk),...
            sprintf('10^{%1.1f}',log10(Epsilon_class.bin(i))),...
            'fontsize',20)
    end
    ll(i)=loglog(Epsilon_class.k,Epsilon_class.mPshear2(i,:),'Color',cmap(i,:),'linewidth',2);
end

xlim([mink maxk])
ylim([minP maxP])
grid on
xlabel('k [cpm]')
ylabel('[s^{-2} / cpm]')
set(gca,'fontsize',20)
title([title_string '- Shear2'])

nanind=find(nansum(Epsilon_class.mPshear2,2)>0);
legend_string=arrayfun(@(x) sprintf('log_{10}(\\epsilon)~%2.1f',x),log10(Epsilon_class.bin),'un',0);
hl=gridLegend(ll(nanind),2,legend_string{nanind},'location','northwest');
set(hl,'fontsize',10)

axes('position',[.7 .1 .25 .82])

plot(log10(Epsilon_class.bin),Epsilon_class.nbin2,'-+k')
set(gca,'fontsize',20)
xlabel('log_{10} \epsilon')
xlim([min(log10(Epsilon_class.bin)) max(log10(Epsilon_class.bin))])
title('PDF \epsilon_2')


