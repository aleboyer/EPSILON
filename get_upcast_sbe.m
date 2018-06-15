function [up,down,dataup,datadown] = get_upcast_sbe(data,crit)


% rename variable to make it easy
if isfield(data,'info')
    info=data.info;
    data=rmfield(data,'info');
end
pdata=-data.P;
tdata=data.time;
L=length(tdata);

% buid a filter 
dt=median(diff(tdata)); % sampling period
T=tdata(end)-tdata(1);  % length of the record
disp('check if time series is shorter than 3 hours')
if T<3/24  
    warning('time serie is less than 3 hours, very short for data processing, watch out the results')
end

disp('smooth the pressure to define up and down cast')
Nb  = 3; % filter order
fnb = 1/(2*dt); % Nyquist frequency
fc  = 1/600/dt; % 100 dt (give "large scale patern") 
[b,a]= butter(Nb,fc/fnb,'low');
filt_pdata=filtfilt(b,a,pdata);
shortdata=filt_pdata(floor(L/3):floor(2*L/3));
ddata=diff(shortdata);
prime_data=filt_pdata;
%crit=1;

Start_ind    =  find(prime_data>=-crit,1,'first');
nb_profile   =  1;
do_it        =  0;

while (do_it==0)
    End_ind=Start_ind+find(prime_data(Start_ind+1:end)>=-crit,1,'first');
    if ~isempty(End_ind)
        zoom_data=prime_data(Start_ind:End_ind);
        enddown_ind=find(zoom_data==min(zoom_data));
        down{nb_profile}=Start_ind:Start_ind+enddown_ind;
        up{nb_profile}=Start_ind+enddown_ind+1:End_ind;
        nb_profile=nb_profile+1;
        Start_ind=End_ind+find(prime_data(End_ind+1:end)<=-crit,1,'first');
    else
        do_it=1;
    end
    if mod(nb_profile,10)==0
        disp(nb_profile)
    end
end


% 
% 
% while (do_it==0)
%     End_ind=Start_ind+find(d_primedata(Start_ind+1:end)>=-crit,1,'first');
%     ind_down(nb_profile)=Start_ind:End_ind+find(d_primedata(Start_ind:End_ind) == ...
%                               min(d_primedata(Start_ind:End_ind)));
%     Start_ind=find(prime_data(End_ind+1:end)>=crit,1,'first');
%     End_ind=Start_ind+find(prime_data(Start_ind+1:end)<=crit,1,'first');
%     if ~isempty(End_ind)
%         ind_up(nb_profile)=Start_ind+find(prime_data(Start_ind:End_ind) == ...
%                                 max(prime_data(Start_ind:End_ind)));
%         nb_profile=nb_profile+1;
%         Start_ind=find(prime_data(End_ind+1:end)<=-crit,1,'first');
%     else
%         do_it=1;
%     end
%     if mod(nb_profile,10)==0
%         %disp(nb_profile)
%     end
% end
% 
% 
% up=arrayfun(@(x,y) (x:y),ind_down(1:end-1),ind_up,'un',0);    
% down=arrayfun(@(x,y) (x:y),ind_up(1:end),ind_down(2:end),'un',0);    
% 
% 
% % up{end+1}=ind_down(end):L;
% % down{end+1}=1:ind_down(end);
% 
 dataup=cellfun(@(x) structfun(@(y) y(x),data,'un',0),up,'un',0);
 datadown=cellfun(@(x) structfun(@(y) y(x),data,'un',0),down,'un',0);

end
