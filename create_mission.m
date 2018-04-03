%% Create folder and meta data structure of an Epsilometer mission
%  The schematic of this structure can be found on confluence

addpath toolbox/
addpath toolbox/FILTER/


%% Define main names 
Meta_Data.mission='DEV';
Meta_Data.vehicle_name='SD';
Meta_Data.deployement='d1';
Meta_Data.path_mission='/Volumes/DataDrive/';
Meta_Data.vehicle='bench';   % 'WireWalker' or 'FastCTD'

%% create mission folders
mission_folder_L0=sprintf('%s/%s',...
                        Meta_Data.path_mission,...
                        Meta_Data.mission);
mission_folder_L1=sprintf('%s/%s',...
                        mission_folder_L0,...
                        Meta_Data.vehicle_name);
mission_folder_L2=sprintf('%s/%s',...
                        mission_folder_L1,...
                        Meta_Data.deployement);
                    
L1path   = sprintf('%s/L1/',mission_folder_L2);
Epsipath = sprintf('%s/epsi/',mission_folder_L2);
CTDpath  = sprintf('%s/ctd/',mission_folder_L2);
RAWpath  = sprintf('%s/raw/',mission_folder_L2);

                                            
if ~exist(mission_folder_L0,'dir')
    %% create paths
    eval([ '!mkdir ' mission_folder_L0]);
    eval([ '!mkdir ' mission_folder_L1]);
    eval([ '!mkdir ' mission_folder_L2]);
    eval([ '!mkdir ' L1path]);
    eval([ '!mkdir ' Epsipath]);
    eval([ '!mkdir ' CTDpath]);
    eval([ '!mkdir ' RAWpath]);
end

%% add path fields
Meta_Data.root     = mission_folder_L2;
Meta_Data.L1path   = L1path;
Meta_Data.Epsipath = Epsipath;
Meta_Data.CTDpath  = CTDpath;
Meta_Data.RAWpath  = RAWpath;


%% add Firmware fields
Meta_Data.Firmware.version='MADRE2.1';
Meta_Data.Firmware.sampling_frequency='320Hz';
Meta_Data.Firmware.ADCshear='Unipolar';
Meta_Data.Firmware.ADC_FPO7='Unipolar';
Meta_Data.Firmware.ADC_accellerometer='Unipolar';

%% add auxillary device field
Meta_Data.aux1.name = 'SBE49';
Meta_Data.aux1.SN   = '00000';
Meta_Data.aux1.cal_file='path/to/file';

%% add channels fields
Meta_Data.epsi.s1.SN='102'; % serial number;
Meta_Data.epsi.s2.SN='102'; % serial number;
Meta_Data.epsi.t1.SN='000'; % serial number;
Meta_Data.epsi.t2.SN='000'; % serial number;

Meta_Data.epsi.shearcal_path='/Users/aleboyer/ARNAUD/SCRIPPS/SHEAR_PROBE/CALIBRATION';
Meta_Data.epsi=get_shear_calibration(Meta_Data.epsi);    % Calibration number

Meta_Data=get_filters_name_MADRE(Meta_Data);

save([Meta_Data.RAWpath ...
    'Meta_' Meta_Data.mission ...
    '_' Meta_Data.deployement '.mat'],'Meta_Data')

