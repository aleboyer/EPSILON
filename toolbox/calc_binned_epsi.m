function Epsilon_class=calc_binned_epsi(MS,epsi_bin)
    %% create epsi_bin or include min and max epsi value to epsi_bin
    if nargin<2
        epsi_bin=10.^(-10:.5:-4.5);
    end
    
    % trick for conditional anonymous function
    iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();

    
    % change struct to cell TODO check if it would be worth it to do it for
    %  the whole process
    if iscell(MS)
        S_MS=[MS{:}];
    end
        
    %% create common k(cpm) axis
    maxk = max([S_MS.k]);
    mink = min([S_MS.k]);
    if mink==0
        mink=sort([S_MS.k]);
        mink=mink(find(mink>0,1,'first'));
    end
    dk   = mink;
    
    k    = mink:dk:maxk;
    
    % project Epsilon fields onto the common k and omega axis 
    Psheark=cellfun(@(x) shiftdim(interp1(x.k,shiftdim(x.Pshear_k,1),k),1) ,MS,'un',0);
    Psheark=[Psheark{:}];
    
    epsilon=cat(1,S_MS.epsilon);
    kvis=real(cat(1,S_MS.kvis));
    
    index1=arrayfun(@(x) find(epsilon(:,1)>=x-.5*x & epsilon(:,1)<x+.5*x),epsi_bin,'un',0);
    index2=arrayfun(@(x) find(epsilon(:,2)>=x-.5*x & epsilon(:,2)<x+.5*x),epsi_bin,'un',0);
    %% if thehre is only 1 epsilon value in the bin, empty the bin.
    index1=cellfun(@(x) iif(length(x)<=1,[],length(x)>1,x),index1,'un',0);
    index2=cellfun(@(x) iif(length(x)<=1,[],length(x)>1,x),index2,'un',0);
    
    Epsilon_class=struct();
    Epsilon_class.k=k;
    Epsilon_class.bin=epsi_bin;
    Epsilon_class.nbin1=cellfun(@length,index1);
    Epsilon_class.nbin2=cellfun(@length,index2);
    Epsilon_class.Pshear1=cellfun(@(x) squeeze(Psheark(1,x,:)),index1,'un',0);
    Epsilon_class.Pshear2=cellfun(@(x) squeeze(Psheark(2,x,:)),index2,'un',0);
    Epsilon_class.mPshear1=cell2mat(cellfun(@(x) nanmean(x,1),Epsilon_class.Pshear1,'un',0).');
    Epsilon_class.mPshear2=cell2mat(cellfun(@(x) nanmean(x,1),Epsilon_class.Pshear2,'un',0).');
    Epsilon_class.kvis=cellfun(@(x) nanmean(kvis(x)),index1,'un',0);
    Epsilon_class.kvis=cellfun(@(x) iif(isnan(x),nanmean([Epsilon_class.kvis{:}]),~isnan(x),x), ...
                   Epsilon_class.kvis,'un',0);
   
    [kpan,Ppan] = cellfun(@(x,y) panchev(x,y),num2cell(epsi_bin),Epsilon_class.kvis,'un',0);
    Epsilon_class.kpan = cell2mat(kpan).';
    Epsilon_class.Ppan = cell2mat(Ppan).';
    
end

