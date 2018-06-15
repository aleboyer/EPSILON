fid=fopen('mpafp07filter6a.txt');
C=textscan(fid,'%s %s %s','Headerlines',1);
freq=char(C{1});
freq=str2num(freq);

valueIn=C{2}(:);
valueIn1=zeros(1,length(valueIn));
for n=1:length(valueIn)
    if valueIn{n}(2)=='-'
        
        valueIn1(n)=str2num(valueIn{n}(2:22));
    else
        valueIn1(n)=str2num(valueIn{n}(2:21));
    end
end

valueOut=C{3}(:);
valueOut1=zeros(1,length(valueOut));
for n=1:length(valueOut)
    if valueOut{n}(2)=='-'
        
        valueOut1(n)=str2num(valueOut{n}(2:22));
    else
        valueOut1(n)=str2num(valueOut{n}(2:21));
    end
end


coef_filt=10.^(valueOut1/20); % value of the spice model given by sean are in dB
coef_filt=coef_filt./max(coef_filt);
save('toolbox/FILTER/FPO7_coeffilt.mat','coef_filt','freq');

