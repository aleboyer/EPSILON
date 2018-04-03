
load('/Users/aleboyer/ARNAUD/SCRIPPS/SHEAR_PROBE/PROCESS_EPSI/PLOT_EPSILOMETER/Drop21503_MMP_EPSI.mat')
A=detrend(EPSI_a1,'constant');
A=A./max(abs(A));
B=detrend(MMP_a1,'constant');
B=B./max(abs(B));

ax(1)=subplot(311);
plot(ax(1),timeaxis,A,'r')
hold on
plot(ax(1),timeaxis,B,'b')
hold off
legend('EPSI','MMP','location','southeast')
title(ax(1),'21503-Normalized acceleration channel','fontsize',20)
ax(1).FontSize=15;


ax(2)=subplot(312);
plot(ax(2),timeaxis,EPSI_v1,'r')
hold off
ylabel(ax(2),'Volt','fontsize',20)
title(ax(2),'EPSI Shear channel','fontsize',20)
ax(2).FontSize=15;

ax(3)=subplot(313);
plot(ax(3),timeaxis,MMP_v1,'b')
title(ax(3),'MMP Shear channel','fontsize',20)
ylabel(ax(3),'Volt','fontsize',20)
ax(3).FontSize=15;
xlabel(ax(3),'minutes','fontsize',20)

fig=gcf;
fig.PaperPosition = [0 0 20 10];
print('-dpng','../../../MMP_EPSI_accel_shear_21503.png')
