multipath = load('trackingResults_PRN19_Inc_-2dB_180Mpath_40s.mat');
noMultipath = load('trackingResults_PRN19_WB_NoMpath_40s.mat');

chipMetres = 299792458/1.023e6;

codePhaseDiff = chipMetres.*(noMultipath.trackResults.codePhase - multipath.trackResults.codePhase);

figure(69)
plot(codePhaseDiff)

figure(70)
plot(noMultipath.trackResults.CNo)
hold on
plot(multipath.trackResults.CNo)
