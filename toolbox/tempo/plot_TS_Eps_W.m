%% plot T-S-rho, epsilon, and wall rate
root_data='/Users/aleboyer/ARNAUD/SCRIPPS/SHEAR_PROBE/';
root_script='/Users/aleboyer/ARNAUD/SCRIPPS/WireWalker/WW_process_profile/';
Cruise_name='EPSIWW'; % 
WW_name='EPSI'; % 
deployement='d1';
WWpath=sprintf('%s/%s/WW/%s/%s/L1/',root_data,Cruise_name,WW_name,deployement);
epsipath=sprintf('%s/%s/WW/%s/%s/epsi/',root_data,Cruise_name,WW_name,deployement);
name_rbr=[WW_name '_rbr_' deployement];



%% addpath 
addpath toolbox/generalplotting/
addpath toolbox/seawater2/

%% load data 
load([epsipath 'EpsiProfile.mat'],'Epsilon')
load([epsipath 'EpsiProfile.mat'],'Profiles')

%
Hfmain=figure;
orient landscape
wysiwyg
ig=get(Hfmain,'position');
set(Hfmain,'position',ig*.85); % to fit on screen better
%
% Define positions for the three panels
pos_tsd=[.1 .2 .4 .75];
pos_eps=[.52 .2 .2 .75];
pos_chi=[.74 .2 .2 .75];
%
% Define position for text under panels 2 & 3
text_pos=[.52 .02 .4 .16];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot temp, salinity, sigma_theta on the left
j=10;
for j=1:length(Profiles)
    salinity   = Profiles{j}.S(:);
    temp       = Profiles{j}.T(:);
    pr_thetasd = Profiles{j}.P(:);
    epsilon    = Epsilon{j}.epsilon(:);
    pr_eps     = Epsilon{j}.pr(:);
    w          = Epsilon{j}.w(:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %MHA 5/2015: use salinity in PSU!!
    theta=sw_ptmp(salinity,temp,pr_thetasd,0);
    sgth=sw_dens(salinity,theta,zeros(size(salinity)))-1000;
    
    
    % Determine plotting limits
    pmin=0;
    ig=find(~isnan(pr_thetasd)); pmax=max(pr_thetasd(ig));
    pmax=1.03*pmax;
    igt=find(~isnan(theta)); igs=find(~isnan(salinity)); igth=find(~isnan(sgth));
    tmin=min(theta(igt)); tmax=max(theta(igt)); tspan=tmax-tmin;
    tmin=tmin-.03*tspan; tmax=tmax+.03*tspan;
    %S=1000*salinity;
    S=salinity;
    smin=min(S(igs)); smax=max(S(igs)); span=smax-smin;
    if span==0; span=1; end
    smin=smin-0.1*span; smax=smax+0.1*span;
    sgmin=min(sgth(igth)); sgmax=max(sgth(igth)); sgspan=sgmax-sgmin;
    sgmin=sgmin-.03*sgspan; sgmax=sgmax+.03*sgspan;
    Tlimits=[tmin tmax pmin pmax];
    Slimits=[smin smax pmin pmax];
    SGlimits=[sgmin sgmax pmin pmax];
    limits=[Tlimits; Slimits; SGlimits];
    %
    ipl_pr  = 1:length(pr_thetasd); % plot only every 5th scan
    ieps_pr = 1:length(pr_eps); % plot only every 5th scan
    X=[theta(ipl_pr) S(ipl_pr) sgth(ipl_pr)];
    Y=pr_thetasd(ipl_pr)*ones(1,3);
    xlabeltext=['       \theta / {}^o C     ';
        '              s            ';
        '\sigma_{\theta} / kg m^{-3}'];
    titletext='';
    ylabeltext='p  / MPa  ';
    linetype=[' r';' g';' b'];
    dy= 0.02;
    axestobelinked = multixaxis3(X,Y,xlabeltext,ylabeltext,titletext,[],[],linetype,[],[],...
        pos_tsd,dy,limits,1,1);
    
    axestobelinked(1).XColor='r';
    axestobelinked(3).XColor='b';
    linkaxes(axestobelinked,'y');
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot epsilon in the middle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Set plotting limits
    h = axes('position',pos_eps,'box','on','yticklabel','','ydir','reverse', ...
        'xscale','log','linewidth',1,'fontsize',12);
    axestobelinked = [axestobelinked h];
    hold on
    Hl_eps1=plot(epsilon,pr_eps,'r');
    set(Hl_eps1,'linewidth',2);
    % Hl_eps2=plot(epsilon(:,2),pr_eps,'g');
    % set(Hl_eps2,'linewidth',[1]);
    igeps=find(pr_eps>.05);
    maxeps=max(max(epsilon(igeps,:)));
    if maxeps>1e-3; maxeps=1e-3; end
    mineps=1e-10; maxeps=10^(ceil(log10(maxeps)));
    if maxeps<=mineps
        maxeps=10*mineps;
    end
    xtick=(floor(log10(mineps)):ceil(log10(maxeps)));
    xtick=10.^xtick;
    
    %
    xlabel('\epsilon / W kg^{-1} (v1-r, v2-g)')
    axis([mineps maxeps pmin pmax])
    set(gca,'xtick',xtick)
    set(gca,'yminortick','on')
    linkaxes(axestobelinked,'y');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot chi on the right
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if exist('chi','var') ~=0
        if ~isempty(tlch)
            h=axes('position',pos_chi,'box','on','yticklabel','','ydir','reverse', ...
                'xscale','log');
            hold on
            Hl_chi1=plot(chi(:,1),pr_chi,'r');
            if length(tlch)==2
                Hl_chi2=plot(chi(:,2),pr_chi,'g');
                % following line causing problems 7/22/96  replace w/below
                %igchi=find(~isnan(chi(:,1)) & ~isnan(chi(:,2)) & pr_chi>.05);
                igchi=find(~isnan(chi(:,1))& pr_chi>.05);
            elseif length(tlch)==1
                igchi=find(~isnan(chi(:,1))& pr_chi>.05);
            end
            %minchi=min(chi(igchi));minchi=max(minchi,1e-11); minchi=10^(floor(log10(minchi)));
            %maxchi=max(chi(igchi));maxchi=min(maxchi,1e-5); maxchi=10^(ceil(log10(maxchi)));
            minchi=1e-10; maxchi=1e-2;
            if minchi>maxchi/100; maxchi=100*minchi; end
            xtick=(floor(log10(minchi)):ceil(log10(maxchi)));
            xtick=10.^xtick;
            %
            set(Hl_chi1,'linewidth',[1.5]);
            xlabel('\chi / K^2 s^{-1}')
            axis([minchi maxchi pmin pmax]);
            axestobelinked = [axestobelinked h];
            set(gca,'xtick',xtick)
            set(gca,'yminortick','on')
        end
        linkaxes(axestobelinked,'y');
        % if we're not calculating chi, plot w in it's own plot
    else
        axes('position',pos_chi,'box','on','yticklabel','','ydir','reverse', ...
            'xscale','linear','linewidth',1,'fontsize',12)
        set(gca,'yminortick','on')
        fpr=0.02; % for plot limits, ignore shallower w's
        if pmax>0.15
            fpr=0.1;
        end
        ig=find(~isnan(w) & pr_eps>fpr);
        if ~isempty(ig)
            wmin=min(w(ig)); wmax=max(w(ig));
            psp = max(pmax-pmin, 0.01); % try to exclude drop-ends from range calc
            %ig = find(pr_thetasd>(pmin+psp/20) & pr_thetasd<(pmax-psp/20) & ~isnan(w));
            ig = find(pr_eps>(pmin+psp/20) & pr_eps<(pmax-psp/20) & ~isnan(w));
            if ~isempty(ig), wmin = max(wmin,min(w(ig))); end
            wspan=max(wmax-wmin,0.01);
            wmin=wmin-.03*wspan; wmax=wmax+.03*wspan;
            hold on
            Hl_w=plot(w(ieps_pr),pr_eps(ieps_pr),'b','linewidth',1);
            xlabel('w / m s^{-1}')
            axis([wmin wmax pmin pmax])
            plotw='y';
        else
            plotw='n';
        end
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Label plot under panels 2 & 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    axes('position',text_pos,'box','off')
    hold on
    axis off
    %
    line1_str=['cruise=' Cruise_name ',  drop=' int2str(j) ',  EPSIWW'];
    Ht1=text(.05,.65, line1_str);
    set(Ht1,'fontsize',12,'fontweight','bold')
    %
    xd=fix(clock);
    pd=yearday(xd(3),xd(2),xd(1),xd(4),xd(5),xd(6));
    %
    %
    year=datestr(Epsilon{1}.timeaxis,'yyyy');
    yday=datenum2yday(Epsilon{1}.timeaxis);
    line2_str=['data acquired: ' year ',  ' num2str(yday,6) ...
        ',  processed:' int2str(xd(1)) ';  ' num2str(pd,6)];
    text(.05,.45,line2_str)
    %
    if exist('last_pbot','var')
        lpb=num2str(last_pbot,5);
    else
        lpb=[];
    end
    line3_str=['max(Pr thetasd)=' num2str(max(pr_thetasd),5) ',  pbot=' lpb];
    text(.05,.30,line3_str)
    %
    
    fig=gcf;
    fig.PaperPosition = [0 0 20 30];
    fig.PaperOrientation='Portrait';
    print(sprintf('../FIGURE/PROFILE/%s_TSEW_Profile%i.png',name_rbr,j),'-dpng2')
    close all
end    
