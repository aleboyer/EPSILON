
loglog(MS{1}.f,squeeze(MS{1}.Pf(3,10:30,:)))
xlabel('Hz','fontsize',15)
ylabel('s^{-2} /Hz','fontsize',15)
title('SP1810 Shear 1','fontsize',15)
print([WWpath deployement '_SP1810_shear1.png'],'-dpng2')

loglog(MS{1}.f,squeeze(MS{1}.Pf(4,10:30,:)))
xlabel('Hz','fontsize',15)
ylabel('s^{-2} /Hz','fontsize',15)
title('SP1810 Shear 2','fontsize',15)
print([WWpath deployement '_SP1810_shear2.png'],'-dpng2')


loglog(MS{1}.f,squeeze(MS{1}.Pf(1,10:30,:)))
xlabel('Hz','fontsize',15)
ylabel('C^{-2} /Hz','fontsize',15)
title('SP1810 temp 1','fontsize',15)
print([WWpath deployement '_SP1810_temp1.png'],'-dpng2')

loglog(MS{1}.f,squeeze(MS{1}.Pf(2,10:30,:)))
xlabel('Hz','fontsize',15)
ylabel('C^{-2} /Hz','fontsize',15)
title('SP1810 temp 2','fontsize',15)
print([WWpath deployement '_SP1810_temp2.png'],'-dpng2')


