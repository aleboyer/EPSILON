
addpath toolbox/
addpath toolbox/FILTER/



fid=fopen('/Users/aleboyer/ARNAUD/SCRIPPS/NISKINE/epsifish2/d6/raw/Meta_NISKINE_d6.dat','r');

count=0;
while(count<10)
    l=fgetl(fid);
    spl=strsplit(l,':');
    Meta_Data.(strtrim(spl{1}))=strtrim(spl{2});
    count=count+1;
end
%% MADRE
fgetl(fid); % empty line
l    = fgetl(fid); % MADRE
spl  = strsplit(l,':');
count=0;
while(count<2)
    l    = fgetl(fid); % MADRE fields
    spl1 = strsplit(l,':');
    Meta_Data.(strtrim(spl{1})).(strtrim(spl1{1}))=strtrim(spl1{2})
    count=count+1;
end

%% MAP
fgetl(fid); % empty line
l    = fgetl(fid); % MAP
spl  = strsplit(l,':');
count=0;
while(count<2)
    l    = fgetl(fid); % MAP fields
    spl1 = strsplit(l,':');
    Meta_Data.(strtrim(spl{1})).(strtrim(spl1{1}))=strtrim(spl1{2})
    count=count+1;
end
%% firmware
fgetl(fid); % empty line
l    = fgetl(fid); % firmware
spl  = strsplit(l,':');
count=0;
while(count<5)
    l    = fgetl(fid); % firmware fields
    spl1 = strsplit(l,':');
    Meta_Data.(strtrim(spl{1})).(strtrim(spl1{1}))=strtrim(spl1{2});
    count=count+1;
end

%% aux1
fgetl(fid); % empty line
l    = fgetl(fid); % aux1
spl  = strsplit(l,':');
count=0;
while(count<3)
    l    = fgetl(fid); % aux1 fields 
    spl1 = strsplit(l,':');
    Meta_Data.(strtrim(spl{1})).(strtrim(spl1{1}))=strtrim(spl1{2});
    count=count+1;
end

%% epsi
fgetl(fid); % empty line
l    = fgetl(fid); % epsi
spl  = strsplit(l,':');
count=0;

%s1 
l    = fgetl(fid); % epsi field s1
spl1 = strsplit(l,':');
l    =fgetl(fid);  %  s1 field SN
spl2 = strsplit(l,':');
Meta_Data.(strtrim(spl{1})).(strtrim(spl1{1})).(strtrim(spl2{1}))=strtrim(spl2{2});
l    =fgetl(fid);  %  s1 field Sv
spl2 = strsplit(l,':');
Meta_Data.(strtrim(spl{1})).(strtrim(spl1{1})).(strtrim(spl2{1}))=strtrim(spl2{2});
l    =fgetl(fid);  %  s1 field ADCfilter
spl2 = strsplit(l,':');
Meta_Data.(strtrim(spl{1})).(strtrim(spl1{1})).(strtrim(spl2{1}))=strtrim(spl2{2});


%s2 
l    = fgetl(fid); % epsi field s2
spl1 = strsplit(l,':');
l    =fgetl(fid);  %  s2 field SN
spl2 = strsplit(l,':');
Meta_Data.(strtrim(spl{1})).(strtrim(spl1{1})).(strtrim(spl2{1}))=strtrim(spl2{2});
l    =fgetl(fid);  %  s2 field Sv
spl2 = strsplit(l,':');
Meta_Data.(strtrim(spl{1})).(strtrim(spl1{1})).(strtrim(spl2{1}))=strtrim(spl2{2});
l    =fgetl(fid);  %  s2 field ADCfilter
spl2 = strsplit(l,':');
Meta_Data.(strtrim(spl{1})).(strtrim(spl1{1})).(strtrim(spl2{1}))=strtrim(spl2{2});


%t1 
l    = fgetl(fid); % epsi field t1
spl1 = strsplit(l,':');
l    =fgetl(fid);  %  t1 field SN
spl2 = strsplit(l,':');
Meta_Data.(strtrim(spl{1})).(strtrim(spl1{1})).(strtrim(spl2{1}))=strtrim(spl2{2});
l    =fgetl(fid);  %  t1 field ADCfilter
spl2 = strsplit(l,':');
Meta_Data.(strtrim(spl{1})).(strtrim(spl1{1})).(strtrim(spl2{1}))=strtrim(spl2{2});


%t2 
l    = fgetl(fid); % epsi field t2
spl1 = strsplit(l,':');
l    =fgetl(fid);  %  t2 field SN
spl2 = strsplit(l,':');
Meta_Data.(strtrim(spl{1})).(strtrim(spl1{1})).(strtrim(spl2{1}))=strtrim(spl2{2});
l    =fgetl(fid);  %  t2 field ADCfilter
spl2 = strsplit(l,':');
Meta_Data.(strtrim(spl{1})).(strtrim(spl1{1})).(strtrim(spl2{1}))=strtrim(spl2{2});

l    =fgetl(fid);  %  t2 field ADCfilter
spl1 = strsplit(l,':');
Meta_Data.(strtrim(spl{1})).(strtrim(spl1{1}))=strtrim(spl1{2});


fclose(fid);

save([Meta_Data.RAWpath ...
    'Meta_' Meta_Data.mission ...
    '_' Meta_Data.deployement '.mat'],'Meta_Data')
