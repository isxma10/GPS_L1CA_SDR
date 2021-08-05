%% This function Generates a dataless GPS Signal with/without Multipath Echo
% Developed by Matthew Alcock (MSc) and Dr Paul Blunt

%% Clean the environment ==================================================

close all;clear all;
addpath include             

%% Name the Signal ========================================================

outputFilename = 'GPS_PRN1_Mulitpath_Amp_0p5_Phase_0_60s_Inc4p0.dat';
%outputFilename = 'GPS_PRN1_No_Multipath_60s.dat';

%% Define Signal Parameters ===============================================

signalLength = 60;                  % Signal length in Seconds              
PRN = 1;                            % Determines the spreading code used
Doppler = 270;                      % Shifts the carrier frequency
multiPathAmp = 0.5;                 % Amplitude relative to LOS signal

% Multipath Phase Selection
multiPathPhaseRad = 0;              % In-Phase
%multiPathPhaseRad = pi;            % Out-of-Phase (180 degrees)
r = rand;
%multiPathPhaseRad = r*pi;          % Random phase between 0 and 180 degs

noiseAmp = 10;                      % Set the noise level

%% Initialise Settings ====================================================

[multiSettings]  = initSettingsMultipath(signalLength);

%% Generate the Signal ====================================================

generateCAcodeMultipath(outputFilename, PRN, Doppler, multiPathAmp, multiPathPhaseRad,noiseAmp, multiSettings);
