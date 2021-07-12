% This example shows how a simple transmitter-channel-receiver simulation may be created using functions from the Matlab LTE Toolbox. The example generates a frame worth of data on one antenna port. As no transport channel is created in this example the data is random bits, QPSK modulated and mapped to every symbol in a subframe. A cell specific reference signal and primary and secondary synchronization signals are created and mapped to the subframe. 10 subframes are individually generated to create a frame. The frame is OFDM modulated, passed through an Extended Vehicular A Model (EVA5) fading channel, additive white Gaussian noise added and demodulated. MMSE equalization using channel and noise estimation is applied and finally the received and equalized resource grids are plotted.

  
% Cell-Wide Settings
% specified in enb, only one transmit antenna is used, referred to TS 36.211 R13 
enb.NDLRB = 15;                 % Number of resource blocks
enb.CellRefP = 1;               % One transmit antenna port
enb.NCellID = 10;               % Cell ID
enb.CyclicPrefix = 'Normal';    % Normal cyclic prefix
enb.DuplexMode = 'FDD';         % FDD


% SNR config
% in decibels
SNRdB = 22;             % Desired SNR in dB
SNR = 10^(SNRdB/20);    % Linear SNR  
rng('default');         % Configure random number generators


% channel model configuration
% specified in cfg structure, fading channel with an EVA delay profile and 120Hz Doppler frequency, paprameters along with MIMO correlation and other channel model specific parameters are set
cfg.Seed = 1;                  % Channel seed
cfg.NRxAnts = 1;               % 1 receive antenna
cfg.DelayProfile = 'EVA';      % EVA delay spread
cfg.DopplerFreq = 120;         % 120Hz Doppler frequency
cfg.MIMOCorrelation = 'Low';   % Low (no) MIMO correlation
cfg.InitTime = 0;              % Initialize at time zero
cfg.NTerms = 16;               % Oscillators used in fading model
cfg.ModelType = 'GMEDS';       % Rayleigh fading model type
cfg.InitPhase = 'Random';      % Random initial phases
cfg.NormalizePathGains = 'On'; % Normalize delay profile power 
cfg.NormalizeTxAnts = 'On';    % Normalize for transmit antennas


% channel estimator configuration
% A user defined window is used to average pilot sysmbols to reduce the effect of noise. The averaging window size is configured in terms of resource elements (REs), in time and frequency. A conservative 9-by-9 window is used in this example as an EVA delay profile and 120Hz Doppler frequency cause the channel changes quickly over time and frequency. A 9-by-9 window includes the 4 pilots immediately surrounding the pilot of interest when averaging. Selecting an averaging window is discussed in Channel Estimation.
cec.PilotAverage = 'UserDefined'; % Pilot averaging method
cec.FreqWindow = 9;               % Frequency averaging window in REs
cec.TimeWindow = 9;               % Time averaging window in REs

% Interpolation is performed by the channel estimator between pilot estimates to create a channel estimate for all REs. To improve the estimate multiple subframes can be used when interpolating. An interpolation window of 3 subframes with a centered interpolation window uses pilot estimates from 3 consecutive subframes to estimate the center subframe.
cec.InterpType = 'Cubic';         % Cubic interpolation
cec.InterpWinSize = 3;            % Interpolate up to 3 subframes 
                                  % simultaneously
cec.InterpWindow = 'Centred';     % Interpolation windowing method


% Subframe Resource grid size
% grid dimensions are determined using lteDLResourceGridSize
gridsize = lteDLResourceGridSize(enb);
K = gridsize(1);    % Number of subcarriers
L = gridsize(2);    % Number of OFDM symbols in one subframe
P = gridsize(3);    % Number of transmit antenna ports


% Transmit Resource Grid
% will be populated with subframes
txGrid = [];


% Payload data generation
% No transport channels is used, data sent are random QPSK modulated symbols
% Number of bits needed is size of resource grid (K*L*P) * number of bits
% per symbol (2 for QPSK)
numberOfBits = K*L*P*2; 
% Create random bit stream
inputBits = randi([0 1], numberOfBits, 1); 
% Modulate input bits
inputSym = lteSymbolModulate(inputBits,'QPSK');


% Frame generation
% Frame will be created by generating individual subframes within a loop and appending each created subframes to the previous subframes. The collection of appended subframes are contained within txGrid. This appending is repeated ten times to create a frame. When the OFDM modulated time domain waveform is passed through a channel, the waveform will experience a delay.
% To avoid any samples being missed due to this delay an extra subframe is generated, therefore 11 subframes are generated in total. For each subframe the Cell-Specific Reference Signal (Cell RS) is added. The Primary Synchronization Signal (PSS) and Secondary Synchronization Signal (SSS) are also added. Note that these synchronization signals only occur in subframes 0 and 5, but the LTE Toolbox takes care of generating empty signals and indices in the other subframes so that the calling syntax here can be completely uniform across the subframes.
% For all subframes within the frame
for sf = 0:10
    
    % Set subframe number
    enb.NSubframe = mod(sf,10);
    
    % Generate empty subframe
    subframe = lteDLResourceGrid(enb);
    
    % Map input symbols to grid
    subframe(:) = inputSym;

    % Generate synchronizing signals
    pssSym = ltePSS(enb); % ltePSS(enb) returns a complex column vector containing the primary synchronization signal (PSS) values for cell-wide settings in the enb structure.
    sssSym = lteSSS(enb); % lteSSS(enb) returns a complex column vector containing the secondary synchronization signal (SSS) values for cell-wide settings in structure enb.
    pssInd = ltePSSIndices(enb); % ltePSSIndices(enb) returns a column vector, ind, of resource element (RE) indices, Port 0 oriented, for the Primary Synchronization Signal (PSS) for the given cell-wide settings structure. By default, the indices are returned in one-based linear indexing form that can directly index elements of a 3-D array representing the resource array. These indices are ordered as the PSS modulation symbols should be mapped. Alternative indexing formats can also be generated
    sssInd = lteSSSIndices(enb); 

    % Map synchronizing signals to the grid
    subframe(pssInd) = pssSym;
    subframe(sssInd) = sssSym;

    % Generate cell specific reference signal symbols and indices
    cellRsSym = lteCellRS(enb); % lteCellRS(enb) returns cell-specific reference signal symbols for cell-wide settings in the enb structure. sym is a complex-valued column vector containing cell-specific reference signal symbols. Unlike other physical channels and signals, the symbols for multiple antennas are concatenated into a single column rather than returned in a matrix with a column for each antenna. The reason for this behavior is that the number of symbols varies across the antenna ports.
    cellRsInd = lteCellRSIndices(enb);

    % Map cell specific reference signal to grid
    subframe(cellRsInd) = cellRsSym;
    
    % Append subframe to grid to be transmitted
    txGrid = [txGrid subframe]; %#ok
      
end



% OFDM Modulation
% Use lteOFDMModulate to transform the frequency domain OFDM symbols into the time domain. Returns two values, a matrix txWaveform and a structure info containing the sampling rate.
[txWaveform,info] = lteOFDMModulate(enb,txGrid); 
txGrid = txGrid(:,1:140);     


% Fading Channel
cfg.SamplingRate = info.SamplingRate;
% Pass data through the fading channel model
rxWaveform = lteFadingChannel(cfg,txWaveform);



% Additive Noise
% The SNR is given by SNR = Es/N0 where Es is the energy of the signal of interest and N0 is the noise power. The noise added before OFDM demodulation will be amplified by the FFT. Therefore to normalize the SNR at the receiver (after OFDM demodulation) the noise must be scaled. The amplification is the square root of the size of the FFT. The size of the FFT can be determined from the sampling rate of the time domain waveform (info.SamplingRate) and the subcarrier spacing (15 kHz). The power of the noise to be added can be scaled so that Es and N0 are normalized after the OFDM demodulation to achieve the desired SNR (SNRdB).
% Calculate noise gain
N0 = 1/(sqrt(2.0*enb.CellRefP*double(info.Nfft))*SNR);
% Create additive white Gaussian noise
noise = N0*complex(randn(size(rxWaveform)),randn(size(rxWaveform)));   
% Add noise to the received time domain waveform
rxWaveform = rxWaveform + noise;


%Synchronization
% function lteDLFrameOffset returns a value offset which indicates how many samples the waveform has been delayed. The offset is considered identical for waveforms received on all antennas. The received time domain waveform can then be manipulated to remove the delay using offset.
offset = lteDLFrameOffset(enb,rxWaveform);
rxWaveform = rxWaveform(1+offset:end,:);

%OFDM Demodulation
% Function lteOFDMDemodulation results 3-dimensional matrix. The number of rows represents the number of subcarriers. The number of columns equals the number of OFDM symbols in a subframe. The number of subcarriers and symbols is the same for the returned grid from OFDM demodulation as the grid passed into lteOFDMModulate. The number of planes (3rd dimension) in the grid corresponds to the number of receive antennas.
rxGrid = lteOFDMDemodulate(enb,rxWaveform);


% Channel estimation
% Function lteDLChannelEstimate is used, configured by the structure cec. 
% lteDLChannelEstimate assumes the first subframe within the resource grid is subframe number enb.NSubframe and therefore the subframe number must be set prior to calling the function
enb.NSubframe = 0;
[estChannel, noiseEst] = lteDLChannelEstimate(enb,cec,rxGrid); % The function returns a 4-D array of complex weights which the channel applies to each resource element in the transmitted grid for each possible transmit and receive antenna combination.
  

  
% MMSE Equalization
% This function uses the estimate of the channel estChannel and noise noiseEst to equalize the received resource grid rxGrid. The function returns eqGrid which is the equalized grid. The dimensions of the equalized grid are the same as the original transmitted grid (txGrid) before OFDM modulation.
eqGrid = lteEqualizeMMSE(rxGrid, estChannel, noiseEst);
  
  
% Anaysis
% Compare the received resource grid with the equalized resource grid. The error between the transmitted and equalized grid and transmitted and received grids are calculated.
% To allow easy inspection the received and equalized grids are plotted on a logarithmic scale using surf within hDownlinkEstimationEqualizationResults.m. These diagrams show that performing channel equalization drastically reduces the error in the received resource grid.  
% Calculate error between transmitted and equalized grid
eqError = txGrid - eqGrid;
rxError = txGrid - rxGrid;

% Compute EVM across all input values
% EVM of pre-equalized receive signal
EVM = comm.EVM;
EVM.AveragingDimensions = [1 2];
preEqualisedEVM = EVM(txGrid,rxGrid);
fprintf('Percentage RMS EVM of Pre-Equalized signal: %0.3f%%\n', ...
        preEqualisedEVM); 
% EVM of post-equalized receive signal
postEqualisedEVM = EVM(txGrid,eqGrid);
fprintf('Percentage RMS EVM of Post-Equalized signal: %0.3f%%\n', ...
        postEqualisedEVM); 
% Plot the received and equalized resource grids 
hDownlinkEstimationEqualizationResults(rxGrid, eqGrid);  
  
  
  
  
