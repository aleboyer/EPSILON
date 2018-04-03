F1=figure;
set(F1,'position',[100 50 1300 700])

ax(1)=subplot(211);
loglog(Chi_class.k,Chi_class.PphiT22{2,2} )
hold(ax(1),'on')
l1=loglog(Chi_class.k,squeeze(Chi_class.mPphiT22(2,2,:)),'c','linewidth',3 );
l2=loglog(Chi_class.kbatch,squeeze(Chi_class.Pbatch22(2,2,:)),'k','linewidth',3);

ax(2)=subplot(212);
loglog(Chi_class.k,Chi_class.PphiT22{3,8} )
hold(ax(2),'on')
l3=loglog(Chi_class.k,squeeze(Chi_class.mPphiT22(3,8,:)),'c','linewidth',3 );
l4=loglog(Chi_class.kbatch,squeeze(Chi_class.Pbatch22(3,8,:)),'k','linewidth',3);


title(ax(1),sprintf('log_{10} \\epsilon= %2.2f log_{10} \\chi= %2.2f',log10(Chi_class.epsibin(2)),log10(Chi_class.chibin(2))),'fontsize',20)
title(ax(2),sprintf('log_{10} \\epsilon= %2.2f log_{10} \\chi= %2.2f',log10(Chi_class.epsibin(8)),log10(Chi_class.chibin(3))),'fontsize',20)
set(ax(1),'fontsize',15)
set(ax(2),'fontsize',15)
xlabel(ax(2),'k /cpm')
ylabel(ax(2),'$(\phi^T_k)^2$ / $degC^2$ $cpm^{-2}$','fontsize',20,'interpreter','latex')
ylabel(ax(1),'$(\phi^T_k)^2$ / $degC^2$ $cpm^{-2}$','fontsize',20,'interpreter','latex')

linkaxes(ax)
legend([l1,l2],{'$\overline{\phi^T_k}$','batchelor'},'location','northwest','interpreter','latex')

print([WWpath 'mar3_example_noise_FPO7_batchelor.png'],'-dpng')
