function [MS]=calc_turbulence_epsi_FastCTD(Profile,tscan,f,Sv,nb_profile)
% calc_eps_mmp.m
%   
%	REQUIRES the following, set up in batchprocess4G_mmp.m, setup_epschi3_mmp.m:
%		drop,procdata,cruise,mmpid, hfperscan,dt_hf,eps_step,FS_hf, f,df;
%		neps,cntr_scan,pr_eps,w_eps, t,s,kvis;  save_eps_spec,displ_shear_spec;
%	CREATES file eps<drop>.mat, saving:  
%       MS structure for Micro Structure. Inside MS you ll find
%       temperature spectra in degC Hz^-1
%       Horizontal  velocity spectra in m^2/s^-2 Hz^-1
%       Acceleration/speed spectra in s^-1 Hz^-1 
%
%
%
% REVISED  June-2017 by Arnaud Le Boyer
% REVISED  June-2000 by Dave Winkel,
%		from 02sep96 M.Gregg version calc_eps3_mmp.m

%% Gravity  ... of the situation :)
G       = 9.81;

%% Length of the Profile
T       = length(Profile.nbsample);
df      = f(1);
%% define number of scan in the profile
Lscan   = tscan*2*f(end);
nbscan  = floor(T/Lscan);

%% we compute spectra on scan with 50% overlap
nbscan=2*nbscan-1;

%% define the fall rate of the Profile. 
%  Add Profile.w with w the vertical vel. 
%  We are using the pressure from other sensors (CTD);
Profile = compute_fallrate_downcast(Profile);
Profile.w=Profile.w/100;

% on line calibration of the FPO7s
CALFPO7_1=polyfit(Profile.Sensor1,Profile.T,3);
CALFPO7_2=polyfit(Profile.Sensor2,Profile.T,3);


%% define the index in the profile for each scan
total_indscan = arrayfun(@(x) (1+floor(Lscan/2)*(x-1):1+floor(Lscan/2)*(x-1)+Lscan-1),1:nbscan,'un',0);
total_w       = cellfun(@(x) nanmean(Profile.w(x)),total_indscan); 

ind_downcast = find((total_w)>.20);

MS.indscan = total_indscan(ind_downcast);
MS.w       = cellfun(@(x) nanmean(Profile.w(x)),total_indscan(ind_downcast)); 
MS.t       = cellfun(@(x) nanmean(Profile.T(x)),total_indscan(ind_downcast)); % needed to compute Kvis 
MS.s       = cellfun(@(x) nanmean(Profile.S(x)),total_indscan(ind_downcast)); % needed to compute Kvis 
MS.pr      = cellfun(@(x) nanmean(Profile.P(x)),total_indscan(ind_downcast)); % needed to compute Kvis 

fmax=45; %defined by get_coherence_epsi.m -> the WW is vibrating around 30 hz 

% loop to estimate eps
% get and split the data

data(1,:,:) = cell2mat(cellfun(@(x) polyval(CALFPO7_1,Profile.Sensor1(x)),MS.indscan,'un',0).');
data(2,:,:) = cell2mat(cellfun(@(x) polyval(CALFPO7_2,Profile.Sensor2(x)),MS.indscan,'un',0).');
data(3,:,:) = cell2mat(cellfun(@(x) Profile.Sensor3(x),MS.indscan,'un',0).');
data(4,:,:) = cell2mat(cellfun(@(x) Profile.Sensor4(x),MS.indscan,'un',0).');
% acceleration divided by speed to get the same units (s-1)
data(5,:,:) = cell2mat(cellfun(@(x,y) Profile.Sensor6(x)./y,MS.indscan,num2cell(MS.w),'un',0).');
data(6,:,:) = cell2mat(cellfun(@(x,y) Profile.Sensor7(x)./y,MS.indscan,num2cell(MS.w),'un',0).');
data(7,:,:) = cell2mat(cellfun(@(x,y) Profile.Sensor8(x)./y,MS.indscan,num2cell(MS.w),'un',0).');

% compute kinematic viscosity
MS.kvis=nu(MS.s,MS.t,MS.pr);
% remove the mean to compute the fft

% Profile Power and Co spectrum and Coherence. (Coherence still needs to be averaged over few scans afterwork)
[f1,P1,P11,Co12]=get_profile_spectrum(data,f);
%TODO comment on the Co12 sturcutre and think about reducing the size of
%the Coherence spectra (doublon)

indf1=find(f1>=0);
indf1=indf1(1:end-1);

P11=2*P11(:,:,indf1);
%% get MADRE filters
h_freq=get_filters_MADRE('MADRE2.1',f1(indf1));

%% correct transfert functions for accel spectra
P11(5:7,:,:)=P11(5:7,:,:)./...
    shiftdim(repmat(ones(nbscan,1)*h_freq.electAccel,[1,1,3]),2);

%% correct transfert functions for shear spectra
TF1 =@(x) (Sv.'.*x/(2*G)).^2 .* h_freq.shear .* haf_oakey(f1(indf1),x);     % should add epsi filter
TFshear=cell2mat(cellfun(@(x) TF1(x),num2cell(MS.w),'un',0).');
TFshear=reshape(TFshear,[2,nbscan,length(indf1)]);
P11(3:4,:,:) = P11(3:4,:,:) ./ TFshear;      % vel frequency spectra m^2/s^-2 Hz^-1

%% correct transfert functions for temperature spectra
TFtemp=cell2mat(cellfun(@(x) h_freq.FPO7(x),num2cell(MS.w),'un',0).');
TFtemp=shiftdim(repmat(TFtemp,[1,1,2]),2);
P11(1:2,:,:) = Emp_Corr_fac * P11(1:2,:,:)./TFtemp;


% convert frequency to wavenumber
k=cellfun(@(x) f/x, num2cell(MS.w),'un',0);
dk=cellfun(@(x) df/x, num2cell(MS.w),'un',0);


Pvelk1  = (P1*speed) ./ htotal(:,1);                                  % vel spec as function of k
Pvelk2  = (P2*speed) ./ htotal(:,2);                                  % vel spec as function of k
MS.Pshear1{j} = (2*pi*k).^2 .* Pvelk1;                           % shear spec  as function of k
MS.Pshear2{j} = (2*pi*k).^2 .* Pvelk2;                           % shear spec  as function of k
MS.k{j}       = k;                                             % k
MS.kvis{j}    = kvis;
        % calc epsilon by integrating to k with 90% variance of Panchev spec
        % unless spectrum is noisy at lower k.
        %
        % Set kmax for integration to highest bin below pump spike,
        % which is between 49 and 52 Hz in a 1024-pt spectrum
        MS.kmax(j)=fmax/speed; % Lowest estimate below pump spike in 1024-pt record
        % Check that data window > 0.5 m, as needed for initial estimate
        if tscan*speed>0.5
            [MS.epsilon1(j),MS.kc1(j)]=eps1_mmp(k,MS.Pshear1{j},kvis,speed,dk,MS.kmax(j));
            [MS.epsilon2(j),MS.kc2(j)]=eps1_mmp(k,MS.Pshear2{j},kvis,speed,dk,MS.kmax(j));
            [MS.kpan1{j},MS.Ppan1{j}] = panchev(MS.epsilon1(j),kvis);
            [MS.kpan2{j},MS.Ppan2{j}] = panchev(MS.epsilon2(j),kvis);
            %
            if mod(nb_profile,1)==0
                % display each spectrum if displ_shear_spec='yes'
                
                ref1(1,1)=max(MS.Pshear1{j}); ref1(2,1)=min(MS.Pshear1{j});
                ref2(1,1)=max(MS.Pshear2{j}); ref2(2,1)=min(MS.Pshear2{j});
                kcref1(1,1)=MS.kc1(j); kcref1(2,1)=MS.kc1(j);
                kcref2(1,1)=MS.kc2(j); kcref2(2,1)=MS.kc2(j);
                kmaxref(1,1)=MS.kmax(j); kmaxref(2,1)=MS.kmax(j);
                loglog(MS.k{j},MS.Pshear1{j},'b');
                hold on
                loglog(MS.k{j},MS.Pshear2{j},'k');
                loglog(kcref1,ref1,'b'); %loglog(kcref,ref,'x')
                loglog(kcref2,ref2,'k'); %loglog(kcref,ref,'x')
                loglog(kmaxref,ref1,'g'); %loglog(kmaxref,ref,'+')
                title(['Profile ' int2str(nb_profile) 'Scan=' int2str(j) ', pr=' ...
                    num2str(MS.pr(j)) ', epsilon1=' num2str(MS.epsilon1(j)) ...
                    ', epsilon2=' num2str(MS.epsilon2(j))]);
                xlabel('k (cpm)','fontsize',15)
                ylabel('\phi(k)','fontsize',15)
                set(gca,'fontsize',15)
                loglog(MS.kpan1{j},MS.Ppan1{j},'bo');
                loglog(MS.kpan2{j},MS.Ppan2{j},'ko');
                %ylim([0.1*min(MS.Ppan1{j}) 1e3*max(MS.Ppan1{j})])
                xlim([1e-1 1e3])
                ylim([1e-9 1e-1])
                hold off
                frame=getframe(fig);
                writeVideo(v,frame)
                clf;
            end
        end
    end
end% end nbscan, end profile processing
if mod(nb_profile,10)==0
    close(v)
    close all
end



