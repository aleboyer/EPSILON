
load('/Users/aleboyer/ARNAUD/SCRIPPS/SHEAR_PROBE/PROCESS_EPSI/data/SPROULnovember/EPSI/d1/EPSITEST_sproul_sbe2.mat')
timesbe=time;
load('/Users/aleboyer/ARNAUD/SCRIPPS/SHEAR_PROBE/PROCESS_EPSI/data/SPROULnovember/EPSI/d1/EPSITEST_sproul_epsi2.mat')

testdate=datetime(timesbe,'ConvertFrom','posixtime');

dataT=interp1(time,Sensor2,timesbe);

CALFPO7=polyfit(dataT,T,3);
test=polyval(CALFPO7,dataT);

matlab_timesbe=datenum(testdate);

subplot('Position',[.1 .8 .8 .15]);
plot(matlab_timesbe,P,'linewidth',2)
set(gca,'XTick',datenum('09-Nov-2017 04:40:00'):5/1440: ...
                matlab_timesbe(end))
set(gca,'XTickLabel','')
axis ij
ylabel('Depth /m','Fontsize',15)
ylim([-10,300])
set(gca,'Fontsize',12)
legend('SBE49 Pressure','location','southwest')

subplot('Position',[.1 .35 .8 .45]);
plot(matlab_timesbe,T,'r','linewidth',2)
hold on
plot(matlab_timesbe,test,'b','linewidth',2)
set(gca,'XTick',datenum('09-Nov-2017 04:40:00'):5/1440: ...
                matlab_timesbe(end))
            set(gca,'XTickLabelRotation',45)
ylim([9,21])            
ylabel('Temperature / C\circ','Fontsize',15)
set(gca,'Fontsize',12)
legend('SBE49-T','FP07')
set(gca,'XTickLabel','')

dataS=interp1(time,Sensor3,timesbe);
subplot('Position',[.1 .2 .8 .15]);
plot(matlab_timesbe,dataS,'b','linewidth',2)
set(gca,'XTick',datenum('09-Nov-2017 04:40:00'):5/1440: ...
                matlab_timesbe(end))
            set(gca,'XTickLabelRotation',45)
set(gca,'XTickLabel',datestr(datenum('09-Nov-2017 04:40:00'):5/1440: ...
                matlab_timesbe(end),'HH:MM:SS'))
ylabel('Shear / Volt','Fontsize',15)
xlabel(datestr(matlab_timesbe(1)),'Fontsize',15)
set(gca,'Fontsize',12)
ylim([0,2.5])


fig=gcf;
fig.PaperPosition = [0 0 10 8];
print('-dpng','../FIGURE/SBE49T_FPO7_EPSISPROUL.png')





