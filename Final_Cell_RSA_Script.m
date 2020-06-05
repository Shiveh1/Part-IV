% This script uses the 2016/2017 Gattoni Cell model with calmodulin modulated
% SERCA
% Modifications by J Musgrave:
% - Does monotonic for 10 mins of sim time and then 5 mins for RSA using
% end monotonic values as ICs
% - Calls a stimulation function (stimfunction) to implement RSA
% - Runs through a range of heart rates and variation amplitudes
% - Saves plots of everything
% - Saves key data
% - Saves parameters at the end of each sim

clear
% Create ALGEBRAIC of correct size
global algebraicVariableCount;  algebraicVariableCount = 88;


varies = linspace(0,1,11); % range of RSA amplitude from 0->1 Hz
global bperiod; 
global stimtimes;
global HR; % mean HR
global vary; % amplitude of RSA

% iterate for different mean HRs
for HR = [5]

% store values for analysis of each HR
outputs=zeros(11,6);

PEAKS = [];

% iterate for different variation amplitudes
for vary=varies
    % Initialise constants and state variables
    if vary == 0 % monotonic sim
    [INIT_STATES, CONSTANTS] = initConsts;
    tlim = 600000; % run a bit longer for monotonic - 10 mins
    else % rsa sim
    INIT_STATES = MONO_STATES;  % use end conditions
    tlim = 300000; % run for 5 mins for RSA
    end

    [~]=stimfunction(0,0,0); % calling stimfunction to get stimtimes and bperiod

    % Set timespan to solve over to be breath period
    tspan = [0, bperiod];

    % Set numerical accuracy options for ODE solver
    options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 1);
  
    i=0;
    TIME=[]; %stores time for whole simulation
    STATES=[]; %states for whole simulation
    iSTATES = []; %states for the current/most recent breath period
    % Solve model with ODE solver

    while i<tlim/bperiod
        %solve odes
        [iTIME, iSTATES] = ode15s(@(VOI, STATES)computeRates(VOI, STATES, CONSTANTS), tspan, INIT_STATES, options);

        len = size(iTIME,1);
        INIT_STATES = iSTATES(len,:);
        TIME=[TIME iTIME'];
        STATES=[STATES iSTATES'];
        tspan = [iTIME(len),iTIME(len)+bperiod];
        i=i+1;
        
        %if vary ~= 0
        %   CONSTANTS(:,55) = (0.01*(1- exp(-0.04*(300/(tlim/bperiod))*i)))  +  0.05; 
        %end
        
    end

    % Compute algebraic variables
    [RATES, ALGEBRAIC] = computeRates(iTIME, iSTATES, CONSTANTS);
    ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, iSTATES, iTIME);

    % use monotonic end state as initial conditions for RSA and
    % save last min of monotonic for graphs
    if vary == 0 
    MONO_STATES=iSTATES(end,:);
    %finding last minute
    [~,idx] = min(abs(TIME+60000-TIME(end)));
    mono_endT = TIME(idx:end)-TIME(idx);
    mono_endS = STATES(:,idx:end);
    end

    %% Analysis
    
    % combining mono and RSA data for graphs (add last min of mono data to
    % beginning
    if vary ~= 0
        TIME = [mono_endT TIME+60000];
        STATES = [mono_endS STATES];
    end
    
    % Free Calcium
    figure
    f=gcf;
    f.WindowState='maximized';
    plot(TIME/1000,STATES(13,:)*1000,'k')
    ylabel('Concentration (uM)')
    xlabel('time (s)')
    title('Free Ca behaviour')
    %saveas(gcf,[num2str(HR),'Hz' num2str(vary) ' Ca.jpg'])

    % Membrane Voltage 
    plot(TIME/1000,STATES(1,:),'k')
    ylabel('Voltage (mV)')
    xlabel('time (s)')
    title('Membrane Voltage behaviour')
    %saveas(gcf,[num2str(HR),'Hz' num2str(vary) ' Voltage.jpg'])

    % Sodium Concentration
    plot(TIME/1000,STATES(2,:),'k')
    ylabel('Concentration (mM)')
    xlabel('time (s)')
    title('Sodium behaviour')
    %saveas(gcf,[num2str(HR),'Hz' num2str(vary) ' Na.jpg'])

    % Potassium Concentration
    plot(TIME/1000,STATES(6,:),'k')
    ylabel('Concentration (mM)')
    xlabel('time (s)')
    title('K behaviour')
    avec = mean(findpeaks(iSTATES(:,6)));
    %saveas(gcf,join([num2str(HR),'Hz' num2str(vary) ' K.jpg']))

    % SERCA Ca Concentration
    plot(TIME/1000,STATES(14,:),'k')
    ylabel('Concentration (mM)')
    xlabel('time (s)')
    title('SR Ca behaviour')
    saveas(gcf,join([num2str(HR),'Hz' num2str(vary) ' SRCa.jpg']))

    % Troponin
    plot(TIME/1000,STATES(18,:),'k')
    ylabel('Concentration (mM)')
    xlabel('time (s)')
    title('Troponin behaviour')
    %saveas(gcf,join([num2str(HR),'Hz' num2str(vary) ' Troponin.jpg']))

    % End time Ca behaviour!
    Ca_SS = iSTATES(:,13)*1000;
    t_SS = iTIME-iTIME(1);
    plot(t_SS,Ca_SS,'k')
    
    % other stats on Ca
    [pks,idx] = findpeaks(Ca_SS);
    avepk = mean(pks);
    PEAKS=[PEAKS,avepk];
    troughs = findpeaks(2-Ca_SS);
    avetrough = -mean(troughs-2);
    aveampl = avepk-avetrough;
    tassel = max(pks)-min(pks);
    auc=trapz(Ca_SS);
    
    for j=1:length(stimtimes)-1
        t = stimtimes(j);
        p = stimtimes(j+1)-t;
        txt = ['\leftarrow ' num2str(p,4) ' ms \rightarrow'];
        text(t+p/2,avepk,txt,'HorizontalAlignment','center')
    end
    txt = ['\leftarrow ' num2str(p,4) ' ms \rightarrow'];
    text(stimtimes(end)+p/2,avepk,txt,'HorizontalAlignment','center')

    ylabel('Concentration (uM)')
    xlabel('time (ms)')
    saveas(gcf,join([num2str(HR),'Hz' num2str(vary) ' Steady State Ca.jpg']))
    close(f)
    
    % index for output
    oi = int8(vary*10+1);
    outputs(oi,1) = oi-1;
    outputs(oi,2) = auc;
    outputs(oi,3) = aveampl;
    outputs(oi,4) = avepk;
    outputs(oi,5) = avetrough;
    outputs(oi,6) = tassel;
    % save all the variables in case we need later
    save([num2str(HR) 'Hz' num2str(vary) '.mat'])
    
    % once we get to last variation, save the outputs matrix for that HR
    if oi == 11
        save([num2str(HR) 'Hz output.mat'],'outputs');
    end
end
    
end
%% Functions below here are part of the cell ML model with very minor edits
% made by JMusgrave
function [STATES, CONSTANTS] = initConsts()
    % This function defines the initial conditions for 6 Hz, used for
    % monotonic simulations
    global HR;
    VOI = 0; CONSTANTS = []; STATES = []; ALGEBRAIC = [];
    STATES(:,1) = -83.5486878017526;
    CONSTANTS(:,1) = 8314;
    CONSTANTS(:,2) = 310;
    CONSTANTS(:,3) = 96487;
    CONSTANTS(:,4) = 0.0001;
    CONSTANTS(:,5) = 1/HR*1000; % heart beat period
    CONSTANTS(:,6) = 3;
    CONSTANTS(:,7) = -0.0012;
    CONSTANTS(:,8) = 25850;
    CONSTANTS(:,9) = 2.585e-5;
    CONSTANTS(:,10) = 2.098e-6;
    CONSTANTS(:,11) = 0.0007;
    STATES(:,2) = 13.6740982817237;
    CONSTANTS(:,12) = 140;
    STATES(:,3) = 0.00264996385357493;
    STATES(:,4) = 0.764915362203661;
    STATES(:,5) = 0.593238313370696;
    CONSTANTS(:,13) = 1.96e-5;
    CONSTANTS(:,14) = 0.883;
    CONSTANTS(:,15) = 0.117;
    CONSTANTS(:,16) = 5.4;
    STATES(:,6) = 136.616455312563;
    STATES(:,7) = 0.00167914037093536;
    STATES(:,8) = 0.942364196302614;
    STATES(:,9) = 0.289751848759626;
    CONSTANTS(:,17) = 1.2e-5;
    STATES(:,10) = 0.00224986598955826;
    STATES(:,11) = 0.253344243516578;
    CONSTANTS(:,18) = 4e-5;
    CONSTANTS(:,19) = 1.45e-6;
    CONSTANTS(:,20) = 0.2;
    STATES(:,12) = 0.00290576866153454;
    CONSTANTS(:,21) = 8.015e-8;
    CONSTANTS(:,22) = 1.38e-7;
    CONSTANTS(:,23) = 0.00138;
    CONSTANTS(:,24) = 3.6;
    CONSTANTS(:,25) = 19;
    CONSTANTS(:,26) = 22;
    CONSTANTS(:,27) = 880;
    CONSTANTS(:,28) = 0.3;
    CONSTANTS(:,29) = 1.8;
    CONSTANTS(:,30) = 1.8;
    CONSTANTS(:,31) = 0.099;
    CONSTANTS(:,32) = 0.02;
    CONSTANTS(:,33) = 0.0007;
    CONSTANTS(:,34) = 50000;
    STATES(:,13) = 0.000247148732859065;
    STATES(:,14) = 1.74941544642849;
    CONSTANTS(:,35) = -9;
    CONSTANTS(:,36) = 7;
    CONSTANTS(:,37) = 11.5;
    CONSTANTS(:,38) = 1;
    CONSTANTS(:,39) = 1450;
    CONSTANTS(:,40) = 1.17;
    CONSTANTS(:,41) = 2.4;
    CONSTANTS(:,42) = 0.05;
    CONSTANTS(:,43) = 0.012;
    CONSTANTS(:,44) = 0.065;
    CONSTANTS(:,45) = 0.0006;
    CONSTANTS(:,46) = 0.0625;
    CONSTANTS(:,47) = 14;
    CONSTANTS(:,48) = 0.01;
    CONSTANTS(:,49) = 100;
    STATES(:,15) = 0.93527246301579;
    STATES(:,16) = 0.00895634833667885;
    STATES(:,17) = 0.0552423715970428;
    CONSTANTS(:,50) = 87.5;
    CONSTANTS(:,51) = 1.38;
    CONSTANTS(:,52) = 0.35;
    CONSTANTS(:,53) = 0.1;
    CONSTANTS(:,54) = 0.0515;
    STATES(:,18) = 0.820228418446148;
    CONSTANTS(:,55) = 0.05;
    CONSTANTS(:,56) = 0.0002;
    CONSTANTS(:,57) = 0.00045;
    CONSTANTS(:,58) = 5e-6;
    CONSTANTS(:,59) = 0.00035;
    CONSTANTS(:,60) = 5.10875e-8;
    CONSTANTS(:,61) = 7.11e-6;
    CONSTANTS(:,62) = 0.04;
    CONSTANTS(:,63) = 40;
    CONSTANTS(:,64) = 0.07;
    STATES(:,19) = 0.0541640622725934;
    CONSTANTS(:,65) = 0.002382;
    CONSTANTS(:,66) = 0.05;
    CONSTANTS(:,67) = 0;
    CONSTANTS(:,68) = 0.00015;
    CONSTANTS(:,69) = 2100.00;
    CONSTANTS(:,70) = 1.00000 - CONSTANTS(:,20);
    CONSTANTS(:,71) = CONSTANTS(:,37)./CONSTANTS(:,38);
    CONSTANTS(:,72) = CONSTANTS(:,42)./CONSTANTS(:,40);
    if (isempty(STATES)), warning('Initial values for states not set'); end
end

function [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS)
    % system of ODEs to solve, includes calculations for lots of algebraic
    % values that are needed to solve
    % RATES(11)=dSTATES(11)/dt
    global algebraicVariableCount;
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        ALGEBRAIC = zeros(1, algebraicVariableCount);
        utilOnes = 1;
    else
        statesRowCount = statesSize(1);
        ALGEBRAIC = zeros(statesRowCount, algebraicVariableCount);
        RATES = zeros(statesRowCount, statesColumnCount);
        utilOnes = ones(statesRowCount, 1);
    end
    ALGEBRAIC(:,9) = 1.00000./(1.00000+exp((STATES(:,1)+87.5000)./10.3000));
    RATES(:,11) = (ALGEBRAIC(:,9) - STATES(:,11))./CONSTANTS(:,69);
    ALGEBRAIC(:,2) = 1.00000./(1.00000+exp((STATES(:,1)+45.0000)./ - 6.50000));
    ALGEBRAIC(:,12) = 1.36000./(( 0.320000.*(STATES(:,1)+47.1300))./(1.00000 - exp(  - 0.100000.*(STATES(:,1)+47.1300)))+ 0.0800000.*exp( - STATES(:,1)./11.0000));
    RATES(:,3) = (ALGEBRAIC(:,2) - STATES(:,3))./ALGEBRAIC(:,12);
    ALGEBRAIC(:,3) = 1.00000./(1.00000+exp((STATES(:,1)+76.1000)./6.07000));
    ALGEBRAIC(:,13) = piecewise({STATES(:,1)>= - 40.0000,  0.453700.*(1.00000+exp( - (STATES(:,1)+10.6600)./11.1000)) }, 3.49000./( 0.135000.*exp( - (STATES(:,1)+80.0000)./6.80000)+ 3.56000.*exp( 0.0790000.*STATES(:,1))+ 310000..*exp( 0.350000.*STATES(:,1))));
    RATES(:,4) = (ALGEBRAIC(:,3) - STATES(:,4))./ALGEBRAIC(:,13);
    ALGEBRAIC(:,4) = 1.00000./(1.00000+exp((STATES(:,1)+76.1000)./6.07000));
    ALGEBRAIC(:,14) = piecewise({STATES(:,1)>= - 40.0000, ( 11.6300.*(1.00000+exp(  - 0.100000.*(STATES(:,1)+32.0000))))./exp(  - 2.53500e-07.*STATES(:,1)) }, 3.49000./( ((STATES(:,1)+37.7800)./(1.00000+exp( 0.311000.*(STATES(:,1)+79.2300)))).*(  - 127140..*exp( 0.244400.*STATES(:,1)) -  3.47400e-05.*exp(  - 0.0439100.*STATES(:,1)))+( 0.121200.*exp(  - 0.0105200.*STATES(:,1)))./(1.00000+exp(  - 0.137800.*(STATES(:,1)+40.1400)))));
    RATES(:,5) = (ALGEBRAIC(:,4) - STATES(:,5))./ALGEBRAIC(:,14);
    ALGEBRAIC(:,15) = 100.000./( 45.1600.*exp( 0.0357700.*(STATES(:,1)+50.0000))+ 98.9000.*exp(  - 0.100000.*(STATES(:,1)+38.0000)));
    ALGEBRAIC(:,5) = 1.00000./(1.00000+exp((STATES(:,1)+10.6000)./ - 11.4200));
    RATES(:,7) = (ALGEBRAIC(:,5) - STATES(:,7))./ALGEBRAIC(:,15);
    ALGEBRAIC(:,16) =  20.0000.*exp( - power((STATES(:,1)+70.0000)./25.0000, 2.00000))+35.0000;
    ALGEBRAIC(:,6) = 1.00000./(1.00000+exp((STATES(:,1)+45.3000)./6.88410));
    RATES(:,8) = (ALGEBRAIC(:,6) - STATES(:,8))./ALGEBRAIC(:,16);
    ALGEBRAIC(:,17) =  1300.00.*exp( - power((STATES(:,1)+70.0000)./30.0000, 2.00000))+35.0000;
    ALGEBRAIC(:,7) = 1.00000./(1.00000+exp((STATES(:,1)+45.3000)./6.88410));
    RATES(:,9) = (ALGEBRAIC(:,7) - STATES(:,9))./ALGEBRAIC(:,17);
    ALGEBRAIC(:,18) = 10000.0./( 45.1600.*exp( 0.0357700.*(STATES(:,1)+50.0000))+ 98.9000.*exp(  - 0.100000.*(STATES(:,1)+38.0000)));
    ALGEBRAIC(:,8) = 1.00000./(1.00000+exp((STATES(:,1)+11.5000)./ - 11.8200));
    RATES(:,10) = (ALGEBRAIC(:,8) - STATES(:,10))./ALGEBRAIC(:,18);
    ALGEBRAIC(:,19) = 1000.00./( 0.118850.*exp((STATES(:,1)+80.0000)./28.3700)+ 0.562300.*exp((STATES(:,1)+80.0000)./ - 14.1900));
    ALGEBRAIC(:,10) = 1.00000./(1.00000+exp((STATES(:,1)+138.600)./10.4800));
    RATES(:,12) = (ALGEBRAIC(:,10) - STATES(:,12))./ALGEBRAIC(:,19);
    ALGEBRAIC(:,24) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./CONSTANTS(:,3)).*log(CONSTANTS(:,16)./STATES(:,6));
    ALGEBRAIC(:,25) =  CONSTANTS(:,13).*STATES(:,7).*( CONSTANTS(:,14).*STATES(:,8)+ CONSTANTS(:,15).*STATES(:,9)).*(STATES(:,1) - ALGEBRAIC(:,24));
    ALGEBRAIC(:,26) =  CONSTANTS(:,17).*STATES(:,10).*STATES(:,11).*(STATES(:,1) - ALGEBRAIC(:,24));
    ALGEBRAIC(:,27) = ( (0.0480000./(exp((STATES(:,1)+37.0000)./25.0000)+exp((STATES(:,1)+37.0000)./ - 25.0000))+0.0100000).*0.00100000)./(1.00000+exp((STATES(:,1) - (ALGEBRAIC(:,24)+76.7700))./ - 17.0000))+( CONSTANTS(:,18).*(STATES(:,1) - (ALGEBRAIC(:,24)+1.73000)))./( (1.00000+exp(( 1.61300.*CONSTANTS(:,3).*(STATES(:,1) - (ALGEBRAIC(:,24)+1.73000)))./( CONSTANTS(:,1).*CONSTANTS(:,2)))).*(1.00000+exp((CONSTANTS(:,16) - 0.998800)./ - 0.124000)));
    ALGEBRAIC(:,32) =  CONSTANTS(:,22).*(STATES(:,1) - ALGEBRAIC(:,24));
    ALGEBRAIC(:,34) =  CONSTANTS(:,24).*power(1.00000+CONSTANTS(:,25)./STATES(:,2), 2.00000).*(1.00000+ (CONSTANTS(:,26)./STATES(:,2)).*exp(((  - CONSTANTS(:,28).*CONSTANTS(:,3).*STATES(:,1))./CONSTANTS(:,1))./CONSTANTS(:,2)))+ power(1.00000+CONSTANTS(:,29)./CONSTANTS(:,16), 2.00000).*(1.00000+ (CONSTANTS(:,12)./CONSTANTS(:,27)).*exp(((  - (1.00000 - CONSTANTS(:,28)).*CONSTANTS(:,3).*STATES(:,1))./CONSTANTS(:,1))./CONSTANTS(:,2)));
    ALGEBRAIC(:,35) = ( CONSTANTS(:,23).*(CONSTANTS(:,24)+1.00000))./ALGEBRAIC(:,34);
    % My (JMusgrave) implementation of the stim function to produce stim
    % current (algebraic 21)
    if (size(VOI)>1)
        for i=1:size(VOI)
            ALGEBRAIC(i,21) = stimfunction(VOI(i),CONSTANTS(:,6),CONSTANTS(:,7));
        end
    else
        ALGEBRAIC(:,21) = stimfunction(VOI,CONSTANTS(:,6),CONSTANTS(:,7));
    end 
    ALGEBRAIC(:,29) =  CONSTANTS(:,19).*STATES(:,12).*CONSTANTS(:,70).*(STATES(:,1) - ALGEBRAIC(:,24));
    RATES(:,6) = (  - ((ALGEBRAIC(:,21)+ALGEBRAIC(:,26)+ALGEBRAIC(:,32)+ALGEBRAIC(:,25)+ALGEBRAIC(:,27)+ALGEBRAIC(:,29)) -  2.00000.*ALGEBRAIC(:,35)).*1.00000)./( CONSTANTS(:,9).*CONSTANTS(:,3));
    ALGEBRAIC(:,1) = ( CONSTANTS(:,3).*STATES(:,1))./( CONSTANTS(:,1).*CONSTANTS(:,2));
    ALGEBRAIC(:,20) =  2.00000.*ALGEBRAIC(:,1);
    ALGEBRAIC(:,45) = piecewise({abs(ALGEBRAIC(:,20))>1.00000e-09, (STATES(:,13)+( (CONSTANTS(:,33)./CONSTANTS(:,31)).*CONSTANTS(:,30).*ALGEBRAIC(:,20).*exp( - ALGEBRAIC(:,20)))./(1.00000 - exp( - ALGEBRAIC(:,20))))./(1.00000+( (CONSTANTS(:,33)./CONSTANTS(:,31)).*ALGEBRAIC(:,20))./(1.00000 - exp( - ALGEBRAIC(:,20)))) }, (STATES(:,13)+ (CONSTANTS(:,33)./CONSTANTS(:,31)).*CONSTANTS(:,30))./(1.00000+CONSTANTS(:,33)./CONSTANTS(:,31)));
    ALGEBRAIC(:,47) = (power(ALGEBRAIC(:,45), 2.00000)+ CONSTANTS(:,48).*power(CONSTANTS(:,44), 2.00000))./( CONSTANTS(:,41).*(power(ALGEBRAIC(:,45), 2.00000)+power(CONSTANTS(:,44), 2.00000)));
    ALGEBRAIC(:,41) = (power(STATES(:,13), 2.00000)+ CONSTANTS(:,48).*power(CONSTANTS(:,44), 2.00000))./( CONSTANTS(:,41).*(power(STATES(:,13), 2.00000)+power(CONSTANTS(:,44), 2.00000)));
    ALGEBRAIC(:,36) = exp((STATES(:,1) - CONSTANTS(:,35))./CONSTANTS(:,36));
    ALGEBRAIC(:,37) = ALGEBRAIC(:,36)./( CONSTANTS(:,38).*(ALGEBRAIC(:,36)+1.00000));
    ALGEBRAIC(:,38) = power(STATES(:,13), 2.00000)./( CONSTANTS(:,40).*(power(STATES(:,13), 2.00000)+power(CONSTANTS(:,44), 2.00000)));
    ALGEBRAIC(:,46) = power(ALGEBRAIC(:,45), 2.00000)./( CONSTANTS(:,40).*(power(ALGEBRAIC(:,45), 2.00000)+power(CONSTANTS(:,44), 2.00000)));
    ALGEBRAIC(:,53) =  (ALGEBRAIC(:,37)+CONSTANTS(:,71)).*( (CONSTANTS(:,71)+CONSTANTS(:,72)+ALGEBRAIC(:,46)).*(CONSTANTS(:,72)+ALGEBRAIC(:,38))+ ALGEBRAIC(:,37).*(CONSTANTS(:,72)+ALGEBRAIC(:,46)));
    ALGEBRAIC(:,54) = ( ALGEBRAIC(:,37).*CONSTANTS(:,72).*(ALGEBRAIC(:,37)+CONSTANTS(:,71)+CONSTANTS(:,72)+ALGEBRAIC(:,38)))./ALGEBRAIC(:,53);
    ALGEBRAIC(:,57) = ( CONSTANTS(:,71).*CONSTANTS(:,72).*(CONSTANTS(:,71)+ALGEBRAIC(:,37)+CONSTANTS(:,72)+ALGEBRAIC(:,46)))./ALGEBRAIC(:,53);
    ALGEBRAIC(:,58) =  ALGEBRAIC(:,54).*ALGEBRAIC(:,47)+ ALGEBRAIC(:,57).*ALGEBRAIC(:,41);
    ALGEBRAIC(:,48) = ( CONSTANTS(:,43).*CONSTANTS(:,49).*(power(ALGEBRAIC(:,45), 2.00000)+ CONSTANTS(:,48).*power(CONSTANTS(:,44), 2.00000)))./( CONSTANTS(:,41).*( CONSTANTS(:,49).*power(ALGEBRAIC(:,45), 2.00000)+ CONSTANTS(:,48).*power(CONSTANTS(:,44), 2.00000)));
    ALGEBRAIC(:,42) = ( CONSTANTS(:,43).*CONSTANTS(:,49).*(power(STATES(:,13), 2.00000)+ CONSTANTS(:,48).*power(CONSTANTS(:,44), 2.00000)))./( CONSTANTS(:,41).*( CONSTANTS(:,49).*power(STATES(:,13), 2.00000)+ CONSTANTS(:,48).*power(CONSTANTS(:,44), 2.00000)));
    ALGEBRAIC(:,60) = ( ALGEBRAIC(:,37).*ALGEBRAIC(:,48)+ CONSTANTS(:,71).*ALGEBRAIC(:,42))./(ALGEBRAIC(:,37)+CONSTANTS(:,71));
    ALGEBRAIC(:,43) = (STATES(:,13)+ (CONSTANTS(:,32)./CONSTANTS(:,31)).*STATES(:,14))./(1.00000+CONSTANTS(:,32)./CONSTANTS(:,31));
    ALGEBRAIC(:,44) = ( ALGEBRAIC(:,43).*(ALGEBRAIC(:,36)+CONSTANTS(:,46)))./( CONSTANTS(:,39).*CONSTANTS(:,45).*(ALGEBRAIC(:,36)+1.00000));
    ALGEBRAIC(:,39) = ( STATES(:,13).*(ALGEBRAIC(:,36)+CONSTANTS(:,46)))./( CONSTANTS(:,39).*CONSTANTS(:,45).*(ALGEBRAIC(:,36)+1.00000));
    ALGEBRAIC(:,55) = ( CONSTANTS(:,71).*( ALGEBRAIC(:,38).*(CONSTANTS(:,71)+CONSTANTS(:,72)+ALGEBRAIC(:,46))+ ALGEBRAIC(:,46).*ALGEBRAIC(:,37)))./ALGEBRAIC(:,53);
    ALGEBRAIC(:,66) =  ALGEBRAIC(:,55).*ALGEBRAIC(:,44)+ ALGEBRAIC(:,57).*ALGEBRAIC(:,39);
    ALGEBRAIC(:,40) = ( CONSTANTS(:,47).*(ALGEBRAIC(:,36)+CONSTANTS(:,46)))./( CONSTANTS(:,39).*( CONSTANTS(:,47).*ALGEBRAIC(:,36)+CONSTANTS(:,46)));
    ALGEBRAIC(:,68) = ALGEBRAIC(:,40);
    RATES(:,15) =   - (ALGEBRAIC(:,58)+ALGEBRAIC(:,66)).*STATES(:,15)+ ALGEBRAIC(:,60).*STATES(:,16)+ ALGEBRAIC(:,68).*STATES(:,17);
    ALGEBRAIC(:,70) = ( CONSTANTS(:,71).*ALGEBRAIC(:,39))./(ALGEBRAIC(:,37)+CONSTANTS(:,71));
    ALGEBRAIC(:,72) = ALGEBRAIC(:,40);
    ALGEBRAIC(:,74) = ((1.00000 - STATES(:,15)) - STATES(:,16)) - STATES(:,17);
    RATES(:,16) = ( ALGEBRAIC(:,58).*STATES(:,15) -  (ALGEBRAIC(:,60)+ALGEBRAIC(:,70)).*STATES(:,16))+ ALGEBRAIC(:,72).*ALGEBRAIC(:,74);
    ALGEBRAIC(:,62) = ( CONSTANTS(:,72).*ALGEBRAIC(:,41))./(CONSTANTS(:,72)+ALGEBRAIC(:,38));
    ALGEBRAIC(:,64) = ALGEBRAIC(:,42);
    RATES(:,17) = ( ALGEBRAIC(:,66).*STATES(:,15) -  (ALGEBRAIC(:,68)+ALGEBRAIC(:,62)).*STATES(:,17))+ ALGEBRAIC(:,64).*ALGEBRAIC(:,74);
    ALGEBRAIC(:,22) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./CONSTANTS(:,3)).*log(CONSTANTS(:,12)./STATES(:,2));
    ALGEBRAIC(:,23) =  CONSTANTS(:,11).*power(STATES(:,3), 3.00000).*STATES(:,4).*STATES(:,5).*(STATES(:,1) - ALGEBRAIC(:,22));
    ALGEBRAIC(:,31) =  CONSTANTS(:,21).*(STATES(:,1) - ALGEBRAIC(:,22));
    ALGEBRAIC(:,75) = ( CONSTANTS(:,54).*( exp( CONSTANTS(:,52).*ALGEBRAIC(:,1)).*power(STATES(:,2), 3.00000).*CONSTANTS(:,30) -  exp( (CONSTANTS(:,52) - 1.00000).*ALGEBRAIC(:,1)).*power(CONSTANTS(:,12), 3.00000).*STATES(:,13)))./( (power(CONSTANTS(:,12), 3.00000)+power(CONSTANTS(:,50), 3.00000)).*(CONSTANTS(:,30)+CONSTANTS(:,51)).*(1.00000+ CONSTANTS(:,53).*exp( (CONSTANTS(:,52) - 1.00000).*ALGEBRAIC(:,1))));
    ALGEBRAIC(:,76) =  ALGEBRAIC(:,75).*CONSTANTS(:,3).*CONSTANTS(:,9);
    ALGEBRAIC(:,28) =  CONSTANTS(:,19).*STATES(:,12).*CONSTANTS(:,20).*(STATES(:,1) - ALGEBRAIC(:,22));
    RATES(:,2) = (  - (ALGEBRAIC(:,23)+ALGEBRAIC(:,31)+ ALGEBRAIC(:,76).*3.00000+ ALGEBRAIC(:,35).*3.00000+ALGEBRAIC(:,28)).*1.00000)./( CONSTANTS(:,9).*CONSTANTS(:,3));
    ALGEBRAIC(:,77) = ( 0.0500000.*(1.00000 - STATES(:,18)).*1.00000)./(1.00000+0.000700000./STATES(:,13));
    RATES(:,18) =  CONSTANTS(:,55).*ALGEBRAIC(:,77).*(ALGEBRAIC(:,77)+STATES(:,18)) -  CONSTANTS(:,56).*STATES(:,18);
    ALGEBRAIC(:,30) = ALGEBRAIC(:,28)+ALGEBRAIC(:,29);
    ALGEBRAIC(:,81) = ( CONSTANTS(:,58).*STATES(:,13))./(CONSTANTS(:,59)+STATES(:,13));
    ALGEBRAIC(:,82) =  ALGEBRAIC(:,81).*2.00000.*CONSTANTS(:,3).*CONSTANTS(:,9);
    ALGEBRAIC(:,52) = piecewise({abs(ALGEBRAIC(:,20))>1.00000e-05, ( (( CONSTANTS(:,33).*ALGEBRAIC(:,20))./(1.00000 - exp( - ALGEBRAIC(:,20)))).*(( CONSTANTS(:,30).*exp( - ALGEBRAIC(:,20)) - STATES(:,13))+ (CONSTANTS(:,32)./CONSTANTS(:,31)).*( CONSTANTS(:,30).*exp( - ALGEBRAIC(:,20)) - STATES(:,14))))./(1.00000+CONSTANTS(:,32)./CONSTANTS(:,31)+( (CONSTANTS(:,33)./CONSTANTS(:,31)).*ALGEBRAIC(:,20))./(1.00000 - exp(ALGEBRAIC(:,20)))) }, ( (( CONSTANTS(:,33).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))).*(( CONSTANTS(:,30).*exp( - 1.00000e-05) - STATES(:,13))+ (CONSTANTS(:,32)./CONSTANTS(:,31)).*( CONSTANTS(:,30).*exp( - 1.00000e-05) - STATES(:,14))))./(1.00000+CONSTANTS(:,32)./CONSTANTS(:,31)+( (CONSTANTS(:,33)./CONSTANTS(:,31)).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))));
    ALGEBRAIC(:,51) = piecewise({abs(ALGEBRAIC(:,20))>1.00000e-05, ( (( CONSTANTS(:,33).*ALGEBRAIC(:,20))./(1.00000 - exp( - ALGEBRAIC(:,20)))).*( CONSTANTS(:,30).*exp( - ALGEBRAIC(:,20)) - STATES(:,13)))./(1.00000+( (CONSTANTS(:,33)./CONSTANTS(:,31)).*ALGEBRAIC(:,20))./(1.00000 - exp( - ALGEBRAIC(:,20)))) }, ( (( CONSTANTS(:,33).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))).*( CONSTANTS(:,30).*exp( - 1.00000e-05) - STATES(:,13)))./(1.00000+( (CONSTANTS(:,33)./CONSTANTS(:,31)).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))));
    ALGEBRAIC(:,56) = ( ALGEBRAIC(:,37).*( ALGEBRAIC(:,46).*(ALGEBRAIC(:,37)+CONSTANTS(:,72)+ALGEBRAIC(:,38))+ ALGEBRAIC(:,38).*CONSTANTS(:,71)))./ALGEBRAIC(:,53);
    ALGEBRAIC(:,67) =  ALGEBRAIC(:,52).*ALGEBRAIC(:,56)+ ALGEBRAIC(:,51).*ALGEBRAIC(:,54);
    ALGEBRAIC(:,69) = ( ALGEBRAIC(:,51).*ALGEBRAIC(:,37))./(ALGEBRAIC(:,37)+CONSTANTS(:,71));
    ALGEBRAIC(:,71) = ( ( STATES(:,15).*ALGEBRAIC(:,67)+ STATES(:,16).*ALGEBRAIC(:,69)).*CONSTANTS(:,34))./CONSTANTS(:,8);
    ALGEBRAIC(:,73) =   - ALGEBRAIC(:,71).*2.00000.*CONSTANTS(:,3).*CONSTANTS(:,9);
    ALGEBRAIC(:,83) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./( 2.00000.*CONSTANTS(:,3))).*log(CONSTANTS(:,30)./STATES(:,13));
    ALGEBRAIC(:,84) =  CONSTANTS(:,60).*(ALGEBRAIC(:,83) - STATES(:,1));
    ALGEBRAIC(:,85) =   - ALGEBRAIC(:,84).*2.00000.*CONSTANTS(:,3).*CONSTANTS(:,9);
    RATES(:,1) =  - (ALGEBRAIC(:,23)+ALGEBRAIC(:,25)+ALGEBRAIC(:,26)+ALGEBRAIC(:,30)+ALGEBRAIC(:,27)+ALGEBRAIC(:,31)+ALGEBRAIC(:,32)+ALGEBRAIC(:,35)+ALGEBRAIC(:,85)+ALGEBRAIC(:,76)+ALGEBRAIC(:,82)+ALGEBRAIC(:,73)+ALGEBRAIC(:,21))./CONSTANTS(:,4);
    ALGEBRAIC(:,49) = ( CONSTANTS(:,32).*(STATES(:,14) - STATES(:,13)))./(1.00000+CONSTANTS(:,32)./CONSTANTS(:,31));
    ALGEBRAIC(:,50) = piecewise({abs(ALGEBRAIC(:,20))>1.00000e-05, ( CONSTANTS(:,32).*((STATES(:,14) - STATES(:,13))+ (( (CONSTANTS(:,33)./CONSTANTS(:,31)).*ALGEBRAIC(:,20))./(1.00000 - exp( - ALGEBRAIC(:,20)))).*(STATES(:,14) -  CONSTANTS(:,30).*exp( - ALGEBRAIC(:,20)))))./(1.00000+CONSTANTS(:,32)./CONSTANTS(:,31)+( (CONSTANTS(:,33)./CONSTANTS(:,31)).*ALGEBRAIC(:,20))./(1.00000 - exp( - ALGEBRAIC(:,20)))) }, ( CONSTANTS(:,32).*((STATES(:,14) - STATES(:,13))+ (( (CONSTANTS(:,33)./CONSTANTS(:,31)).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))).*(STATES(:,14) -  CONSTANTS(:,30).*exp( - 1.00000e-05))))./(1.00000+CONSTANTS(:,32)./CONSTANTS(:,31)+( (CONSTANTS(:,33)./CONSTANTS(:,31)).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))));
    ALGEBRAIC(:,59) =  ALGEBRAIC(:,56).*ALGEBRAIC(:,50)+ ALGEBRAIC(:,49).*ALGEBRAIC(:,55);
    ALGEBRAIC(:,61) = ( ALGEBRAIC(:,49).*ALGEBRAIC(:,38))./(CONSTANTS(:,72)+ALGEBRAIC(:,38));
    ALGEBRAIC(:,63) = ( ( STATES(:,15).*ALGEBRAIC(:,59)+ STATES(:,17).*ALGEBRAIC(:,61)).*CONSTANTS(:,34))./CONSTANTS(:,8);
    ALGEBRAIC(:,65) = ALGEBRAIC(:,63);
    ALGEBRAIC(:,78) = ALGEBRAIC(:,77)+STATES(:,18);
    ALGEBRAIC(:,79) = ( 0.00115800.*ALGEBRAIC(:,78))./(2.22400 - ALGEBRAIC(:,78));
    ALGEBRAIC(:,80) = ( ALGEBRAIC(:,79).*power(STATES(:,13), 2.00000))./(power(CONSTANTS(:,57), 2.00000)+power(STATES(:,13), 2.00000));
    ALGEBRAIC(:,86) =  CONSTANTS(:,61).*(STATES(:,14) - STATES(:,13));
    RATES(:,14) =  (CONSTANTS(:,9)./CONSTANTS(:,10)).*(( - ALGEBRAIC(:,65)+ALGEBRAIC(:,80)) - ALGEBRAIC(:,86));
    ALGEBRAIC(:,87) =  CONSTANTS(:,62).*(CONSTANTS(:,64) - STATES(:,19)) -  CONSTANTS(:,63).*STATES(:,19).*STATES(:,13);
    RATES(:,19) = ALGEBRAIC(:,87);
    ALGEBRAIC(:,88) = power(1.00000+( CONSTANTS(:,65).*CONSTANTS(:,66))./power(CONSTANTS(:,65)+STATES(:,13), 2.00000)+( CONSTANTS(:,67).*CONSTANTS(:,68))./power(CONSTANTS(:,68)+STATES(:,13), 2.00000),  - 1.00000);
    RATES(:,13) =  ALGEBRAIC(:,88).*(((ALGEBRAIC(:,65) - ALGEBRAIC(:,80))+ALGEBRAIC(:,86)+ALGEBRAIC(:,87)) - (  - 2.00000.*ALGEBRAIC(:,76)+ALGEBRAIC(:,73)+ALGEBRAIC(:,82)+ALGEBRAIC(:,85))./( 2.00000.*CONSTANTS(:,9).*CONSTANTS(:,3)));
   RATES = RATES';
end

function ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI)
    % Calculates algebraic variables using calculated state values
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        utilOnes = 1;
    else
        statesRowCount = statesSize(1);
        utilOnes = ones(statesRowCount, 1);
    end
    ALGEBRAIC(:,9) = 1.00000./(1.00000+exp((STATES(:,1)+87.5000)./10.3000));
    ALGEBRAIC(:,2) = 1.00000./(1.00000+exp((STATES(:,1)+45.0000)./ - 6.50000));
    ALGEBRAIC(:,12) = 1.36000./(( 0.320000.*(STATES(:,1)+47.1300))./(1.00000 - exp(  - 0.100000.*(STATES(:,1)+47.1300)))+ 0.0800000.*exp( - STATES(:,1)./11.0000));
    ALGEBRAIC(:,3) = 1.00000./(1.00000+exp((STATES(:,1)+76.1000)./6.07000));
    ALGEBRAIC(:,13) = piecewise({STATES(:,1)>= - 40.0000,  0.453700.*(1.00000+exp( - (STATES(:,1)+10.6600)./11.1000)) }, 3.49000./( 0.135000.*exp( - (STATES(:,1)+80.0000)./6.80000)+ 3.56000.*exp( 0.0790000.*STATES(:,1))+ 310000..*exp( 0.350000.*STATES(:,1))));
    ALGEBRAIC(:,4) = 1.00000./(1.00000+exp((STATES(:,1)+76.1000)./6.07000));
    ALGEBRAIC(:,14) = piecewise({STATES(:,1)>= - 40.0000, ( 11.6300.*(1.00000+exp(  - 0.100000.*(STATES(:,1)+32.0000))))./exp(  - 2.53500e-07.*STATES(:,1)) }, 3.49000./( ((STATES(:,1)+37.7800)./(1.00000+exp( 0.311000.*(STATES(:,1)+79.2300)))).*(  - 127140..*exp( 0.244400.*STATES(:,1)) -  3.47400e-05.*exp(  - 0.0439100.*STATES(:,1)))+( 0.121200.*exp(  - 0.0105200.*STATES(:,1)))./(1.00000+exp(  - 0.137800.*(STATES(:,1)+40.1400)))));
    ALGEBRAIC(:,15) = 100.000./( 45.1600.*exp( 0.0357700.*(STATES(:,1)+50.0000))+ 98.9000.*exp(  - 0.100000.*(STATES(:,1)+38.0000)));
    ALGEBRAIC(:,5) = 1.00000./(1.00000+exp((STATES(:,1)+10.6000)./ - 11.4200));
    ALGEBRAIC(:,16) =  20.0000.*exp( - power((STATES(:,1)+70.0000)./25.0000, 2.00000))+35.0000;
    ALGEBRAIC(:,6) = 1.00000./(1.00000+exp((STATES(:,1)+45.3000)./6.88410));
    ALGEBRAIC(:,17) =  1300.00.*exp( - power((STATES(:,1)+70.0000)./30.0000, 2.00000))+35.0000;
    ALGEBRAIC(:,7) = 1.00000./(1.00000+exp((STATES(:,1)+45.3000)./6.88410));
    ALGEBRAIC(:,18) = 10000.0./( 45.1600.*exp( 0.0357700.*(STATES(:,1)+50.0000))+ 98.9000.*exp(  - 0.100000.*(STATES(:,1)+38.0000)));
    ALGEBRAIC(:,8) = 1.00000./(1.00000+exp((STATES(:,1)+11.5000)./ - 11.8200));
    ALGEBRAIC(:,19) = 1000.00./( 0.118850.*exp((STATES(:,1)+80.0000)./28.3700)+ 0.562300.*exp((STATES(:,1)+80.0000)./ - 14.1900));
    ALGEBRAIC(:,10) = 1.00000./(1.00000+exp((STATES(:,1)+138.600)./10.4800));
    ALGEBRAIC(:,24) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./CONSTANTS(:,3)).*log(CONSTANTS(:,16)./STATES(:,6));
    ALGEBRAIC(:,25) =  CONSTANTS(:,13).*STATES(:,7).*( CONSTANTS(:,14).*STATES(:,8)+ CONSTANTS(:,15).*STATES(:,9)).*(STATES(:,1) - ALGEBRAIC(:,24));
    ALGEBRAIC(:,26) =  CONSTANTS(:,17).*STATES(:,10).*STATES(:,11).*(STATES(:,1) - ALGEBRAIC(:,24));
    ALGEBRAIC(:,27) = ( (0.0480000./(exp((STATES(:,1)+37.0000)./25.0000)+exp((STATES(:,1)+37.0000)./ - 25.0000))+0.0100000).*0.00100000)./(1.00000+exp((STATES(:,1) - (ALGEBRAIC(:,24)+76.7700))./ - 17.0000))+( CONSTANTS(:,18).*(STATES(:,1) - (ALGEBRAIC(:,24)+1.73000)))./( (1.00000+exp(( 1.61300.*CONSTANTS(:,3).*(STATES(:,1) - (ALGEBRAIC(:,24)+1.73000)))./( CONSTANTS(:,1).*CONSTANTS(:,2)))).*(1.00000+exp((CONSTANTS(:,16) - 0.998800)./ - 0.124000)));
    ALGEBRAIC(:,32) =  CONSTANTS(:,22).*(STATES(:,1) - ALGEBRAIC(:,24));
    ALGEBRAIC(:,34) =  CONSTANTS(:,24).*power(1.00000+CONSTANTS(:,25)./STATES(:,2), 2.00000).*(1.00000+ (CONSTANTS(:,26)./STATES(:,2)).*exp(((  - CONSTANTS(:,28).*CONSTANTS(:,3).*STATES(:,1))./CONSTANTS(:,1))./CONSTANTS(:,2)))+ power(1.00000+CONSTANTS(:,29)./CONSTANTS(:,16), 2.00000).*(1.00000+ (CONSTANTS(:,12)./CONSTANTS(:,27)).*exp(((  - (1.00000 - CONSTANTS(:,28)).*CONSTANTS(:,3).*STATES(:,1))./CONSTANTS(:,1))./CONSTANTS(:,2)));
    ALGEBRAIC(:,35) = ( CONSTANTS(:,23).*(CONSTANTS(:,24)+1.00000))./ALGEBRAIC(:,34);
    for i=1:size(VOI)
        ALGEBRAIC(i,21) = stimfunction(VOI(i),CONSTANTS(:,6),CONSTANTS(:,7));
    end
    ALGEBRAIC(:,29) =  CONSTANTS(:,19).*STATES(:,12).*CONSTANTS(:,70).*(STATES(:,1) - ALGEBRAIC(:,24));
    ALGEBRAIC(:,1) = ( CONSTANTS(:,3).*STATES(:,1))./( CONSTANTS(:,1).*CONSTANTS(:,2));
    ALGEBRAIC(:,20) =  2.00000.*ALGEBRAIC(:,1);
    ALGEBRAIC(:,45) = piecewise({abs(ALGEBRAIC(:,20))>1.00000e-09, (STATES(:,13)+( (CONSTANTS(:,33)./CONSTANTS(:,31)).*CONSTANTS(:,30).*ALGEBRAIC(:,20).*exp( - ALGEBRAIC(:,20)))./(1.00000 - exp( - ALGEBRAIC(:,20))))./(1.00000+( (CONSTANTS(:,33)./CONSTANTS(:,31)).*ALGEBRAIC(:,20))./(1.00000 - exp( - ALGEBRAIC(:,20)))) }, (STATES(:,13)+ (CONSTANTS(:,33)./CONSTANTS(:,31)).*CONSTANTS(:,30))./(1.00000+CONSTANTS(:,33)./CONSTANTS(:,31)));
    ALGEBRAIC(:,47) = (power(ALGEBRAIC(:,45), 2.00000)+ CONSTANTS(:,48).*power(CONSTANTS(:,44), 2.00000))./( CONSTANTS(:,41).*(power(ALGEBRAIC(:,45), 2.00000)+power(CONSTANTS(:,44), 2.00000)));
    ALGEBRAIC(:,41) = (power(STATES(:,13), 2.00000)+ CONSTANTS(:,48).*power(CONSTANTS(:,44), 2.00000))./( CONSTANTS(:,41).*(power(STATES(:,13), 2.00000)+power(CONSTANTS(:,44), 2.00000)));
    ALGEBRAIC(:,36) = exp((STATES(:,1) - CONSTANTS(:,35))./CONSTANTS(:,36));
    ALGEBRAIC(:,37) = ALGEBRAIC(:,36)./( CONSTANTS(:,38).*(ALGEBRAIC(:,36)+1.00000));
    ALGEBRAIC(:,38) = power(STATES(:,13), 2.00000)./( CONSTANTS(:,40).*(power(STATES(:,13), 2.00000)+power(CONSTANTS(:,44), 2.00000)));
    ALGEBRAIC(:,46) = power(ALGEBRAIC(:,45), 2.00000)./( CONSTANTS(:,40).*(power(ALGEBRAIC(:,45), 2.00000)+power(CONSTANTS(:,44), 2.00000)));
    ALGEBRAIC(:,53) =  (ALGEBRAIC(:,37)+CONSTANTS(:,71)).*( (CONSTANTS(:,71)+CONSTANTS(:,72)+ALGEBRAIC(:,46)).*(CONSTANTS(:,72)+ALGEBRAIC(:,38))+ ALGEBRAIC(:,37).*(CONSTANTS(:,72)+ALGEBRAIC(:,46)));
    ALGEBRAIC(:,54) = ( ALGEBRAIC(:,37).*CONSTANTS(:,72).*(ALGEBRAIC(:,37)+CONSTANTS(:,71)+CONSTANTS(:,72)+ALGEBRAIC(:,38)))./ALGEBRAIC(:,53);
    ALGEBRAIC(:,57) = ( CONSTANTS(:,71).*CONSTANTS(:,72).*(CONSTANTS(:,71)+ALGEBRAIC(:,37)+CONSTANTS(:,72)+ALGEBRAIC(:,46)))./ALGEBRAIC(:,53);
    ALGEBRAIC(:,58) =  ALGEBRAIC(:,54).*ALGEBRAIC(:,47)+ ALGEBRAIC(:,57).*ALGEBRAIC(:,41);
    ALGEBRAIC(:,48) = ( CONSTANTS(:,43).*CONSTANTS(:,49).*(power(ALGEBRAIC(:,45), 2.00000)+ CONSTANTS(:,48).*power(CONSTANTS(:,44), 2.00000)))./( CONSTANTS(:,41).*( CONSTANTS(:,49).*power(ALGEBRAIC(:,45), 2.00000)+ CONSTANTS(:,48).*power(CONSTANTS(:,44), 2.00000)));
    ALGEBRAIC(:,42) = ( CONSTANTS(:,43).*CONSTANTS(:,49).*(power(STATES(:,13), 2.00000)+ CONSTANTS(:,48).*power(CONSTANTS(:,44), 2.00000)))./( CONSTANTS(:,41).*( CONSTANTS(:,49).*power(STATES(:,13), 2.00000)+ CONSTANTS(:,48).*power(CONSTANTS(:,44), 2.00000)));
    ALGEBRAIC(:,60) = ( ALGEBRAIC(:,37).*ALGEBRAIC(:,48)+ CONSTANTS(:,71).*ALGEBRAIC(:,42))./(ALGEBRAIC(:,37)+CONSTANTS(:,71));
    ALGEBRAIC(:,43) = (STATES(:,13)+ (CONSTANTS(:,32)./CONSTANTS(:,31)).*STATES(:,14))./(1.00000+CONSTANTS(:,32)./CONSTANTS(:,31));
    ALGEBRAIC(:,44) = ( ALGEBRAIC(:,43).*(ALGEBRAIC(:,36)+CONSTANTS(:,46)))./( CONSTANTS(:,39).*CONSTANTS(:,45).*(ALGEBRAIC(:,36)+1.00000));
    ALGEBRAIC(:,39) = ( STATES(:,13).*(ALGEBRAIC(:,36)+CONSTANTS(:,46)))./( CONSTANTS(:,39).*CONSTANTS(:,45).*(ALGEBRAIC(:,36)+1.00000));
    ALGEBRAIC(:,55) = ( CONSTANTS(:,71).*( ALGEBRAIC(:,38).*(CONSTANTS(:,71)+CONSTANTS(:,72)+ALGEBRAIC(:,46))+ ALGEBRAIC(:,46).*ALGEBRAIC(:,37)))./ALGEBRAIC(:,53);
    ALGEBRAIC(:,66) =  ALGEBRAIC(:,55).*ALGEBRAIC(:,44)+ ALGEBRAIC(:,57).*ALGEBRAIC(:,39);
    ALGEBRAIC(:,40) = ( CONSTANTS(:,47).*(ALGEBRAIC(:,36)+CONSTANTS(:,46)))./( CONSTANTS(:,39).*( CONSTANTS(:,47).*ALGEBRAIC(:,36)+CONSTANTS(:,46)));
    ALGEBRAIC(:,68) = ALGEBRAIC(:,40);
    ALGEBRAIC(:,70) = ( CONSTANTS(:,71).*ALGEBRAIC(:,39))./(ALGEBRAIC(:,37)+CONSTANTS(:,71));
    ALGEBRAIC(:,72) = ALGEBRAIC(:,40);
    ALGEBRAIC(:,74) = ((1.00000 - STATES(:,15)) - STATES(:,16)) - STATES(:,17);
    ALGEBRAIC(:,62) = ( CONSTANTS(:,72).*ALGEBRAIC(:,41))./(CONSTANTS(:,72)+ALGEBRAIC(:,38));
    ALGEBRAIC(:,64) = ALGEBRAIC(:,42);
    ALGEBRAIC(:,22) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./CONSTANTS(:,3)).*log(CONSTANTS(:,12)./STATES(:,2));
    ALGEBRAIC(:,23) =  CONSTANTS(:,11).*power(STATES(:,3), 3.00000).*STATES(:,4).*STATES(:,5).*(STATES(:,1) - ALGEBRAIC(:,22));
    ALGEBRAIC(:,31) =  CONSTANTS(:,21).*(STATES(:,1) - ALGEBRAIC(:,22));
    ALGEBRAIC(:,75) = ( CONSTANTS(:,54).*( exp( CONSTANTS(:,52).*ALGEBRAIC(:,1)).*power(STATES(:,2), 3.00000).*CONSTANTS(:,30) -  exp( (CONSTANTS(:,52) - 1.00000).*ALGEBRAIC(:,1)).*power(CONSTANTS(:,12), 3.00000).*STATES(:,13)))./( (power(CONSTANTS(:,12), 3.00000)+power(CONSTANTS(:,50), 3.00000)).*(CONSTANTS(:,30)+CONSTANTS(:,51)).*(1.00000+ CONSTANTS(:,53).*exp( (CONSTANTS(:,52) - 1.00000).*ALGEBRAIC(:,1))));
    ALGEBRAIC(:,76) =  ALGEBRAIC(:,75).*CONSTANTS(:,3).*CONSTANTS(:,9);
    ALGEBRAIC(:,28) =  CONSTANTS(:,19).*STATES(:,12).*CONSTANTS(:,20).*(STATES(:,1) - ALGEBRAIC(:,22));
    ALGEBRAIC(:,77) = ( 0.0500000.*(1.00000 - STATES(:,18)).*1.00000)./(1.00000+0.000700000./STATES(:,13));
    ALGEBRAIC(:,30) = ALGEBRAIC(:,28)+ALGEBRAIC(:,29);
    ALGEBRAIC(:,81) = ( CONSTANTS(:,58).*STATES(:,13))./(CONSTANTS(:,59)+STATES(:,13));
    ALGEBRAIC(:,82) =  ALGEBRAIC(:,81).*2.00000.*CONSTANTS(:,3).*CONSTANTS(:,9);
    ALGEBRAIC(:,52) = piecewise({abs(ALGEBRAIC(:,20))>1.00000e-05, ( (( CONSTANTS(:,33).*ALGEBRAIC(:,20))./(1.00000 - exp( - ALGEBRAIC(:,20)))).*(( CONSTANTS(:,30).*exp( - ALGEBRAIC(:,20)) - STATES(:,13))+ (CONSTANTS(:,32)./CONSTANTS(:,31)).*( CONSTANTS(:,30).*exp( - ALGEBRAIC(:,20)) - STATES(:,14))))./(1.00000+CONSTANTS(:,32)./CONSTANTS(:,31)+( (CONSTANTS(:,33)./CONSTANTS(:,31)).*ALGEBRAIC(:,20))./(1.00000 - exp(ALGEBRAIC(:,20)))) }, ( (( CONSTANTS(:,33).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))).*(( CONSTANTS(:,30).*exp( - 1.00000e-05) - STATES(:,13))+ (CONSTANTS(:,32)./CONSTANTS(:,31)).*( CONSTANTS(:,30).*exp( - 1.00000e-05) - STATES(:,14))))./(1.00000+CONSTANTS(:,32)./CONSTANTS(:,31)+( (CONSTANTS(:,33)./CONSTANTS(:,31)).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))));
    ALGEBRAIC(:,51) = piecewise({abs(ALGEBRAIC(:,20))>1.00000e-05, ( (( CONSTANTS(:,33).*ALGEBRAIC(:,20))./(1.00000 - exp( - ALGEBRAIC(:,20)))).*( CONSTANTS(:,30).*exp( - ALGEBRAIC(:,20)) - STATES(:,13)))./(1.00000+( (CONSTANTS(:,33)./CONSTANTS(:,31)).*ALGEBRAIC(:,20))./(1.00000 - exp( - ALGEBRAIC(:,20)))) }, ( (( CONSTANTS(:,33).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))).*( CONSTANTS(:,30).*exp( - 1.00000e-05) - STATES(:,13)))./(1.00000+( (CONSTANTS(:,33)./CONSTANTS(:,31)).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))));
    ALGEBRAIC(:,56) = ( ALGEBRAIC(:,37).*( ALGEBRAIC(:,46).*(ALGEBRAIC(:,37)+CONSTANTS(:,72)+ALGEBRAIC(:,38))+ ALGEBRAIC(:,38).*CONSTANTS(:,71)))./ALGEBRAIC(:,53);
    ALGEBRAIC(:,67) =  ALGEBRAIC(:,52).*ALGEBRAIC(:,56)+ ALGEBRAIC(:,51).*ALGEBRAIC(:,54);
    ALGEBRAIC(:,69) = ( ALGEBRAIC(:,51).*ALGEBRAIC(:,37))./(ALGEBRAIC(:,37)+CONSTANTS(:,71));
    ALGEBRAIC(:,71) = ( ( STATES(:,15).*ALGEBRAIC(:,67)+ STATES(:,16).*ALGEBRAIC(:,69)).*CONSTANTS(:,34))./CONSTANTS(:,8);
    ALGEBRAIC(:,73) =   - ALGEBRAIC(:,71).*2.00000.*CONSTANTS(:,3).*CONSTANTS(:,9);
    ALGEBRAIC(:,83) =  (( CONSTANTS(:,1).*CONSTANTS(:,2))./( 2.00000.*CONSTANTS(:,3))).*log(CONSTANTS(:,30)./STATES(:,13));
    ALGEBRAIC(:,84) =  CONSTANTS(:,60).*(ALGEBRAIC(:,83) - STATES(:,1));
    ALGEBRAIC(:,85) =   - ALGEBRAIC(:,84).*2.00000.*CONSTANTS(:,3).*CONSTANTS(:,9);
    ALGEBRAIC(:,49) = ( CONSTANTS(:,32).*(STATES(:,14) - STATES(:,13)))./(1.00000+CONSTANTS(:,32)./CONSTANTS(:,31));
    ALGEBRAIC(:,50) = piecewise({abs(ALGEBRAIC(:,20))>1.00000e-05, ( CONSTANTS(:,32).*((STATES(:,14) - STATES(:,13))+ (( (CONSTANTS(:,33)./CONSTANTS(:,31)).*ALGEBRAIC(:,20))./(1.00000 - exp( - ALGEBRAIC(:,20)))).*(STATES(:,14) -  CONSTANTS(:,30).*exp( - ALGEBRAIC(:,20)))))./(1.00000+CONSTANTS(:,32)./CONSTANTS(:,31)+( (CONSTANTS(:,33)./CONSTANTS(:,31)).*ALGEBRAIC(:,20))./(1.00000 - exp( - ALGEBRAIC(:,20)))) }, ( CONSTANTS(:,32).*((STATES(:,14) - STATES(:,13))+ (( (CONSTANTS(:,33)./CONSTANTS(:,31)).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))).*(STATES(:,14) -  CONSTANTS(:,30).*exp( - 1.00000e-05))))./(1.00000+CONSTANTS(:,32)./CONSTANTS(:,31)+( (CONSTANTS(:,33)./CONSTANTS(:,31)).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))));
    ALGEBRAIC(:,59) =  ALGEBRAIC(:,56).*ALGEBRAIC(:,50)+ ALGEBRAIC(:,49).*ALGEBRAIC(:,55);
    ALGEBRAIC(:,61) = ( ALGEBRAIC(:,49).*ALGEBRAIC(:,38))./(CONSTANTS(:,72)+ALGEBRAIC(:,38));
    ALGEBRAIC(:,63) = ( ( STATES(:,15).*ALGEBRAIC(:,59)+ STATES(:,17).*ALGEBRAIC(:,61)).*CONSTANTS(:,34))./CONSTANTS(:,8);
    ALGEBRAIC(:,65) = ALGEBRAIC(:,63);
    ALGEBRAIC(:,78) = ALGEBRAIC(:,77)+STATES(:,18);
    ALGEBRAIC(:,79) = ( 0.00115800.*ALGEBRAIC(:,78))./(2.22400 - ALGEBRAIC(:,78));
    ALGEBRAIC(:,80) = ( ALGEBRAIC(:,79).*power(STATES(:,13), 2.00000))./(power(CONSTANTS(:,57), 2.00000)+power(STATES(:,13), 2.00000));
    ALGEBRAIC(:,86) =  CONSTANTS(:,61).*(STATES(:,14) - STATES(:,13));
    ALGEBRAIC(:,87) =  CONSTANTS(:,62).*(CONSTANTS(:,64) - STATES(:,19)) -  CONSTANTS(:,63).*STATES(:,19).*STATES(:,13);
    ALGEBRAIC(:,88) = power(1.00000+( CONSTANTS(:,65).*CONSTANTS(:,66))./power(CONSTANTS(:,65)+STATES(:,13), 2.00000)+( CONSTANTS(:,67).*CONSTANTS(:,68))./power(CONSTANTS(:,68)+STATES(:,13), 2.00000),  - 1.00000);
    ALGEBRAIC(:,11) = STATES(:,13);
    ALGEBRAIC(:,33) = ALGEBRAIC(:,31)+ALGEBRAIC(:,32);
end

% Compute result of a piecewise function
function x = piecewise(cases, default)
    set = [0];
    for i = 1:2:length(cases)
        if (length(cases{i+1}) == 1)
            x(cases{i} & ~set,:) = cases{i+1};
        else
            x(cases{i} & ~set,:) = cases{i+1}(cases{i} & ~set);
        end
        set = set | cases{i};
        if(set), break, end
    end
    if (length(default) == 1)
        x(~set,:) = default;
    else
        x(~set,:) = default(~set);
    end
end

