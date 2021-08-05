function [multiSettings]  = initSettingsMultipath(signalLength)
%Functions initializes and saves multiSettings. multiSettings can be edited inside of
%the function, updated from the command line or updated using a dedicated
%GUI - "setmultiSettings".  
%
%All multiSettings are described inside function code.
%
%multiSettings = initmultiSettingsNSL_26MHz()
%
%   Inputs: none
%
%   Outputs:
%       multiSettings     - Fundemental Signal Processing Characteristics 
%       

%% Processing multiSettings ====================================================
% Number of milliseconds to be processed used 36000 + any transients (see
% below - in Nav parameters) to ensure nav subframes are provided
multiSettings.msToProcess        = signalLength*1000;        %[ms]

% Data type used to store one sample
multiSettings.dataType           = 'int8';

% Intermediate, sampling and code frequencies
multiSettings.IF                 = 14.58e6;      %[Hz]
multiSettings.samplingFreq       = 53e6;     %[Hz]
multiSettings.codeFreqBasis      = 1.023e6;      %[Hz]
% Account for any spectrum inversion by the RF front end
multiSettings.spectrumInversion = 1;
% Define number of chips in a code period
multiSettings.codeLength         = 1023;
