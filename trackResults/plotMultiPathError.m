%% This function Generates Multipath Error Envelopes between clean and 
%% Multipath Signals
% Created by Matthew Alcock and Dr Paul Blunt

%% Clean the environment ==================================================

close all;clear all;

%% Load Tracking Results= =================================================

multipath = load('trackResultsMultipath_GPS_1_Amp_0p5_0_Inc4_55s_EL_0p16.mat');
noMultipath = load('trackResultsNoMultipath_GPS_1_55s_EL_0p16.mat');

%% Initialize tracking variables ==========================================

secondsTracked = 55;

samplingFreq = 53e6;
codeFreq = 1.023e6;
sampleShiftsPerInc = 5;
incrementTime = 4;
chipsInSeconds = incrementTime*((samplingFreq/codeFreq)/sampleShiftsPerInc);

chipsTracked = secondsTracked/chipsInSeconds;

chipMetres = 299792458/1.023e6;
samples = secondsTracked*10;
MpathShift = floor(samples/40);


%% Initialise Result Structure

Results.meanCodeError = zeros(1, MpathShift);
Results.stdCodeError = zeros(1 , MpathShift);

%% Calculate the Code Phase Difference ====================================

codePhaseDiff = chipMetres.*(noMultipath.trackResults.codePhase - multipath.trackResults.codePhase);
codePhaseDiff(1) = 0;

%% Check for chip roll over ===============================================

for i = 1:1:secondsTracked*10
    if  codePhaseDiff(i) > 150       % Over half a chip
        codePhaseDiff(i) = codePhaseDiff(i) -293;         % remove the chip
    end
    if  codePhaseDiff(i) < -150
        codePhaseDiff(i) = codePhaseDiff(i) + 293;
    end
end

%% Average Error every 4 seconds ==========================================
for n = 1:1:MpathShift
    meanCodeError = mean(codePhaseDiff(((n-1)*40)+25:((n-1)*40)+35));
    stdCodeError = std(codePhaseDiff(((n-1)*40)+25:((n-1)*40)+35))/sqrt(n);
    Results.meanCodeError(n) = meanCodeError;
    Results.stdCodeError(n) = stdCodeError;
    err1 = Results.stdCodeError;
end
%% Plot Results ===========================================================

t = linspace(0,chipsTracked,MpathShift);

figure(59)
subplot(2,2,1)
plot(codePhaseDiff,'r')
title('Raw Multipath Error')
xlabel('Time ms)')
ylabel('Range Error in Metres')
subplot(2,1,2)
plot(t,Results.meanCodeError,'g');
hold on
errorbar(t,Results.meanCodeError,err1,'.g');
title('Filtered Multipath Error')
xlabel('Multipath in Chips')
ylabel('Range Error in Metres')
