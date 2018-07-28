    fig=figure('Position',[100,100,1600,2000]);
    v = VideoWriter('toto.avi');
    v.FrameRate=10;
    open(v)

    for i=1:10
        scatter(i,i,50,'k','filled')
        xlim([1,10])
        ylim([1,10])
        frame=getframe(fig);
        writeVideo(v,frame)
        clf;
    end
    
    close(v)
    close all

