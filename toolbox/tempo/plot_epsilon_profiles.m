%load profile epsi from EPSI SPROUL

figure
for i=1:5
    subplot(1,5,i)
    hold on
    semilogx(Epsilon{i}.epsilon1,Epsilon{i}.pr,'r')
    semilogx(Epsilon{i}.epsilon2,Epsilon{i}.pr,'b')
    hold off
    axis ij
    xlim([0 1e-6])
    ylim([0 275])
    title(['Profile ' num2str(i)])
    set(gca,'Xscale','log')
    xlabel('\epsilon','fontsize',15)
    if i==1
        ylabel('Pressure (db)','fontsize',15)
        legend('Probe1','Probe2','location','south')
    end
    disp(num2str(i))
    set(gca,'fontsize',15)
end


fig=gcf;
fig.PaperPosition = [0 0 15 8];
print('-dpng2','Epsi_SPROUL_epsilon_Profile_d3.png')