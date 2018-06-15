function H=get_filters_MADRE(version,f)
%% Define the electronic filter, ADC filters.
%  The electronics filter Helec = charge amp filter (Hca) + Gain (Hg)
%  The ADC filter is a simple sinc^4 and the gain of the ADC is 2. Other options are available. 
%
%  TODO: enter argument to be able to change the ADC filter without hard
%  coding
%
switch version
    case{'MADRE2.1'}
        
        % shear channels
        %charge amp filter
        
        ca_filter = load('FILTER/charge_coeffilt.mat');
        epsi_ca   = interp1(ca_filter.freq,ca_filter.coef_filt ,f);
        gain_ca      = .5; %TODO check the charge amp gain with sean
        helectronics= epsi_ca*gain_ca;% charge amp from sean spec sheet
        
        sinc4_freq=f/(2*f(end));
        
        H.gainshear=1;
        H.adcshear=(H.gainshear.*(sinc(sinc4_freq)).^4);
        H.electshear=helectronics.^2;
        H.shear=H.electshear .* H.adcshear.^2;
        
        %% FPO7 channels
        H.gainFPO7=1;
        H.electFPO7 = H.gainFPO7.*(sinc(sinc4_freq)).^4;
        
        %speed% convert to m/s
        %tau=0.005 * speed^(-0.32); % thermistor time constant
        H.magsq=@(speed)(1 ./ (1+((2*pi*(0.005 * speed^(-0.32))).*f).^2).^2); % magnitude-squared
        H.phase=@(speed)(-2*atan( 2*pi*f*(0.005 * speed^(-0.32))));
        H.FPO7=@(speed)(H.electFPO7.^2 .* H.magsq(speed));

        %% Accel channels
        H.gainAccel  = 1;
        H.electAccel = (H.gainAccel.*(sinc(sinc4_freq)).^4).^2;

end

