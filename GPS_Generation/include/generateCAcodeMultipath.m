function generateCAcodeMultipath(outputFilename, PRN, Doppler, multiAmp, multiPhase, noiseAmp, settings)
% Performs generates a statice CA code signal with multipath
%% Initialize result structure ============================================
%===============

codePeriods = settings.msToProcess;     % For GPS one C/A code is one ms

samplesPerChip = settings.samplingFreq/settings.codeFreqBasis;
% initialise waitbar
hwb = waitbar(0,'generating...');

% open the output file
fidOut = fopen(outputFilename, 'wb');


%% Start processing ==============================================
    
    % Get a vector with the C/A code sampled 1x/chip
    caCode = generateCAcode(PRN);
    % Then make it possible to do early and late versions
    caCode = [caCode(1023) caCode caCode(1)];

    %--- Perform various initializations ------------------------------

    % define initial code frequency basis of NCO
    % define residual code phase (in chips)
    remCodePhase  = 0.0;
    % define carrier frequency which is used over whole tracking period
    carrFreq      = settings.IF + Doppler;
    codeFreq      = settings.codeFreqBasis + Doppler/1540;
    % define residual carrier phase
    remCarrPhase  = 0.0;

    navBit = 1;
    sampleShift = 0;
    
    %=== Process the number of specified code periods =================
    for loopCnt =  1:codePeriods

%% GUI update -------------------------------------------------------------
        % The GUI is updated every 50ms. This way Matlab GUI is still
        % responsive enough. At the same time Matlab is not occupied
        % all the time with GUI task.
        if (rem(loopCnt, 50) == 0)
            try
                waitbar(loopCnt/codePeriods, ...
                        hwb, ...
                        ['generating PRN#', int2str(PRN), ...
                        '; Completed ',int2str(loopCnt), ...
                        ' of ', int2str(codePeriods), ' ms']);                       
            catch
                % The progress bar was closed. It is used as a signal
                % to stop, "cancel" processing. Exit.
                disp('Progress bar closed, exiting...');
                return
            end
        end

%% Read next block of data ------------------------------------------------            
        % Find the size of a "block" or code period in whole samples

        % Update the phasestep based on code freq (variable) and
        % sampling frequency (fixed)
        codePhaseStep = codeFreq / settings.samplingFreq;

        blksize = ceil((settings.codeLength-remCodePhase) / codePhaseStep);

%% Set up all the code phase tracking information -------------------------

        % Define index into prompt code vector
        tcode       = remCodePhase : ...
                      codePhaseStep : ...
                      ((blksize-1)*codePhaseStep+remCodePhase);
        tcode2      = ceil(tcode) + 1;
        promptCode  = caCode(tcode2);

        remCodePhase = (tcode(blksize) + codePhaseStep) - 1023.0;

%% Generate the carrier frequency to mix the signal to baseband -----------
        time    = (0:blksize) ./ settings.samplingFreq;

        % Get the argument to sin/cos functions
        trigarg = ((carrFreq * 2.0 * pi) .* time) + remCarrPhase;
        remCarrPhase = rem(trigarg(blksize+1), (2 * pi));
       
        % Start incrementing Mpath at 9 seconds
        if loopCnt > 9000
            % increment multipath delay every 4s
            if (rem(loopCnt,4000)== 0)
                sampleShift = sampleShift + 5;
            end
        end
        if (rem(loopCnt,20) ==0)
            navBit = navBit.*-1;
        end

        % Finally compute the signal to mix the collected data to bandband
        modulatedSig = cos(trigarg(1:blksize)).*promptCode.*navBit + multiAmp.*cos(trigarg(1:blksize)+multiPhase).*circshift(promptCode, sampleShift).*navBit;
        %noMpathSig = cos(trigarg(1:blksize)).*promptCode.*navBit;
        % add noise
        noisySig = modulatedSig + randn(1, blksize).*noiseAmp;
        %noisySig = noMpathSig + randn(1, blksize).*noiseAmp;
        %% scale to 8-bits (3 sigma in the range so, std = 42)        
        scale = 42/std(noisySig);
        dataOut = noisySig*scale;
        
        % write the output data to a file
        fwrite(fidOut, dataOut, 'int8');
end % if a PRN is assigned

% Close the waitbar
close(hwb);
fclose(fidOut);
