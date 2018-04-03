function [f1raw,P11raw]=calc_raw_epsi(Profile,tscan,f)
% calc_eps_mmp.m
%   
%	REQUIRES the following, set up in batchprocess4G_mmp.m, setup_epschi3_mmp.m:
%		drop,procdata,cruise,mmpid, hfperscan,dt_hf,eps_step,FS_hf, f,df;
%		neps,cntr_scan,pr_eps,w_eps, t,s,kvis;  save_eps_spec,displ_shear_spec;
%	CREATES file eps<drop>.mat, saving:  
%       Praw structure for Micro Structure. Inside Praw you ll find
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
% Profile = compute_fallrate_downcast(Profile);
% Profile.w=Profile.w/100;
 Profile.w=Profile.Sensor1*0+.5; % trick for bench data
% 
% %% define the index in the profile for each scan
 total_indscan = arrayfun(@(x) (1+floor(Lscan/2)*(x-1):1+floor(Lscan/2)*(x-1)+Lscan-1),1:nbscan,'un',0);
 total_w       = cellfun(@(x) nanmean(Profile.w(x)),total_indscan); 
% 
 ind_downcast = find((total_w)>.20);
% 
 Praw.indscan   = total_indscan(ind_downcast);
 Praw.nbscan    = nbscan;
% Praw.fmax      = 45; % arbitrary cut off frequency usually extract from coherence spectra shear/accel 
% Praw.nbchannel = nbscan;
% 
% Praw.w       = cellfun(@(x) nanmean(Profile.w(x)),total_indscan(ind_downcast)); 
% Praw.t       = cellfun(@(x) nanmean(Profile.T(x)),total_indscan(ind_downcast)); % needed to compute Kvis 
% Praw.s       = cellfun(@(x) nanmean(Profile.S(x)),total_indscan(ind_downcast)); % needed to compute Kvis 
% Praw.pr      = cellfun(@(x) nanmean(Profile.P(x)),total_indscan(ind_downcast)); % needed to compute Kvis 
% 

% loop to estimate eps
% get and split the data
data(1,:,:) = cell2mat(cellfun(@(x) Profile.Sensor1(x),Praw.indscan,'un',0).');
data(2,:,:) = cell2mat(cellfun(@(x) Profile.Sensor2(x),Praw.indscan,'un',0).');
data(3,:,:) = cell2mat(cellfun(@(x) Profile.Sensor3(x),Praw.indscan,'un',0).');
data(4,:,:) = cell2mat(cellfun(@(x) Profile.Sensor4(x),Praw.indscan,'un',0).');
% acceleration divided by speed to get the same units (s-1)
data(5,:,:) = cell2mat(cellfun(@(x) Profile.Sensor6(x),Praw.indscan,'un',0).');
data(6,:,:) = cell2mat(cellfun(@(x) Profile.Sensor7(x),Praw.indscan,'un',0).');
data(7,:,:) = cell2mat(cellfun(@(x) Profile.Sensor8(x),Praw.indscan,'un',0).');

% Profile Power and Co spectrum and Coherence. (Coherence still needs to be averaged over few scans afterwork)
[f1raw,P1raw,P11raw,Co12raw]=get_profile_spectrum(data,f);
%TODO comment on the Co12 sturcutre and think about reducing the size of
%the Coherence spectra (doublon)

indf1=find(f1raw>=0);
indf1=indf1(1:end-1);
f1raw=f1raw(indf1);
Lf1raw=length(indf1);

P11raw=2*P11raw(:,:,indf1);
