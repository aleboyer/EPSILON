function [up,down,dataup,datadown,info] = get_upcast_rbr(data)


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
meandata=nanmean(shortdata(abs(ddata)>=.5*nanmax(abs(ddata))));
prime_data=filt_pdata-meandata;


Start_ind    =  find(prime_data<=0,1,'first');
nb_profile   =  1;
sign_profile = -1;
do_it        =  0;
while (do_it==0)
    End_ind=Start_ind+find(sign_profile*prime_data(Start_ind+1:end)<=0,1,'first');
    ind_down(nb_profile)=Start_ind+find(prime_data(Start_ind:End_ind) == ...
                              min(prime_data(Start_ind:End_ind)));
    sign_profile=sign_profile*-1;
    Start_ind=End_ind;
    End_ind=Start_ind+find(sign_profile*prime_data(Start_ind+1:end)<=0,1,'first');
    if ~isempty(End_ind)
        ind_up(nb_profile)=Start_ind+find(prime_data(Start_ind:End_ind) == ...
                                max(prime_data(Start_ind:End_ind)));
        nb_profile=nb_profile+1;
        sign_profile=sign_profile*-1;
        Start_ind=End_ind;
    else
        do_it=1;
    end
    if mod(nb_profile,10)==0
        disp(nb_profile)
    end
end
    
up=arrayfun(@(x,y) (x:y),ind_down(1:end-1),ind_up,'un',0);    
down=arrayfun(@(x,y) (x:y),ind_up(1:end-1),ind_down(2:end-1),'un',0);    


% up{end+1}=ind_down(end):L;
% down{end+1}=1:ind_down(end);

dataup=cellfun(@(x) structfun(@(y) y(x),data,'un',0),up,'un',0);
datadown=cellfun(@(x) structfun(@(y) y(x),data,'un',0),down,'un',0);

end
