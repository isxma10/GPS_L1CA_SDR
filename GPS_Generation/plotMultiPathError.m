addpath Tracking_Results
multipath = load('trackingResultsMpath_Inc_2_Amp_0p5_Phase_86_50s.mat');
noMultipath = load('trackingResultsNoMpath_50s.mat');

chipMetres = 299792458/1.023e6;

codePhaseDiff = chipMetres.*(noMultipath.trackResults.codePhase - multipath.trackResults.codePhase);

figure(69)
plot(codePhaseDiff)

figure(70)
plot(noMultipath.trackResults.CNo,'r')
hold on
plot(multipath.trackResults.CNo,'g')
