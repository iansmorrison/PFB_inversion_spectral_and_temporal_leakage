
%% GLOBAL PARAMETERS

% PFB parameters
N = 256;  % number of output channels - power of 2, min OS_Nu, max 256
OS_Nu = 4;  % OS numerator - should be a sub-multiple of N
OS_De = 3;  % OS denominator
taps_per_chan = 12;
quantisePFBcoeff = 16;  % choose between 0 (no quantisation), 8, 12, 13, 14, 15 or 16

updateFilterDesign = 1;   % 1 to update filter design, 0 otherwise
displayFilterDesign = 1;  % 1 to display filter design plot, 0 otherwise

% Length of forward FFT to process fine channels
fftLength = 2^8;

% Length of test vector block
blockLength = fftLength*N*OS_De/OS_Nu;

quantisePFBoutput = 0;   % 1 to quantise, 0 to not

equaliseRipple = 1;  % 1 to equalise PFB ripple, 0 to not

% Number of PFB filter arms
PFBarms = (N*OS_De)/OS_Nu;

% Metrics
amplitudeTF = zeros(PFBarms,1);
maxSpuriousPower = zeros(PFBarms,1);
totalSpuriousPower = zeros(PFBarms,1);
averageSpuriousPower = zeros(blockLength,1);


%% DESIGN PFB PROTOTYPE FILTER

Ntaps = N*taps_per_chan + 1;  % must be odd; one more than N*taps_per_chan
    
if (updateFilterDesign)
    disp('updating PFB filter design');

    fprintf('designing PFB prototype filter\n');
    if (displayFilterDesign)
        fprintf('\nPress any key to continue...\n');
    end;
    designPFB(N,OS_Nu,OS_De,Ntaps-1,fftLength,quantisePFBcoeff,displayFilterDesign);  % subtract 1 from num taps because design functions adds 1
end;


%% MAIN LOOP OVER ALL PFB ARMS

for impulseOffset = 1 : PFBarms
    
    % Print loop number
    fprintf('\nPFB Arm %d of %d\n',impulseOffset,PFBarms);

    % GENERATE TEST VECTOR (input to PFB)
    % impulse of amplitude 1.0
    % generate double block length: will use middle half, after PFB primed
    disp('generating test vector');
    inputSig = complex(zeros(2*blockLength,1,'single'));
    impulse_offset = blockLength + impulseOffset;  % location of impulse (first is at blockLength + 1)
    impulse_width = 1;  % number of samples width of impusle
    inputSig(impulse_offset:impulse_offset+impulse_width-1,1) = 1.0;

    % figure;
    % subplot(211); plot((1:2*blockLength),real(inputSig(1:2*blockLength,1))); box on; grid on; title('sig Real'); 
    % subplot(212); plot((1:2*blockLength),imag(inputSig(1:2*blockLength,1))); box on; grid on; title('sig Imag'); xlabel('time');

    spectrum = fft(inputSig)./(2*blockLength);

    % figure;
    % subplot(211); plot((1:2*blockLength),abs(spectrum(1:2*blockLength,1))); box on; grid on; title('sig Mag'); 
    % subplot(212); plot((1:2*blockLength),phase(spectrum(1:2*blockLength,1))); box on; grid on; title('sig Phase'); xlabel('frequency');


    %% PFB Channelize - one block

    disp('channelising test vector');

    % derive PFB output fine channel vector length
    NFineOut = floor(blockLength*OS_Nu/OS_De/N);

    output = complex(zeros(2*NFineOut,N));

    output = PFBchannelizerCSIRO(transpose(inputSig),N,OS_Nu,OS_De);

    % scale up by number of channels
    output = output.*N;

    % optionally quantise to 8 or fewer bits
    if (quantisePFBoutput)
        intOutput = int8(output.*(127.0/2.0));
        output = double(intOutput).*(0.00787401574803*2.0);
    end

    % offset to start of plot/fft
    inputOffset = fftLength/2;

    % plot first 8 channels
    % figure;
    % subplot(8,2,1); plot((inputOffset+1:inputOffset+NFineOut),real(output(inputOffset+1:inputOffset+NFineOut,1))); box on; grid on; title('Channel 1 Real'); xlabel('time');
    % subplot(8,2,2); plot((inputOffset+1:inputOffset+NFineOut),imag(output(inputOffset+1:inputOffset+NFineOut,1))); box on; grid on; title('Channel 1 Imag'); xlabel('time');
    % subplot(8,2,3); plot((inputOffset+1:inputOffset+NFineOut),real(output(inputOffset+1:inputOffset+NFineOut,2))); box on; grid on; title('Channel 2 Real'); xlabel('time');
    % subplot(8,2,4); plot((inputOffset+1:inputOffset+NFineOut),imag(output(inputOffset+1:inputOffset+NFineOut,2))); box on; grid on; title('Channel 2 Imag'); xlabel('time');
    % subplot(8,2,5); plot((inputOffset+1:inputOffset+NFineOut),real(output(inputOffset+1:inputOffset+NFineOut,3))); box on; grid on; title('Channel 3 Real'); xlabel('time');
    % subplot(8,2,6); plot((inputOffset+1:inputOffset+NFineOut),imag(output(inputOffset+1:inputOffset+NFineOut,3))); box on; grid on; title('Channel 3 Imag'); xlabel('time');
    % subplot(8,2,7); plot((inputOffset+1:inputOffset+NFineOut),real(output(inputOffset+1:inputOffset+NFineOut,4))); box on; grid on; title('Channel 4 Real'); xlabel('time');
    % subplot(8,2,8); plot((inputOffset+1:inputOffset+NFineOut),imag(output(inputOffset+1:inputOffset+NFineOut,4))); box on; grid on; title('Channel 4 Imag'); xlabel('time');
    % subplot(8,2,9); plot((inputOffset+1:inputOffset+NFineOut),real(output(inputOffset+1:inputOffset+NFineOut,5))); box on; grid on; title('Channel 5 Real'); xlabel('time');
    % subplot(8,2,10); plot((inputOffset+1:inputOffset+NFineOut),imag(output(inputOffset+1:inputOffset+NFineOut,5))); box on; grid on; title('Channel 5 Imag'); xlabel('time');
    % subplot(8,2,11); plot((inputOffset+1:inputOffset+NFineOut),real(output(inputOffset+1:inputOffset+NFineOut,6))); box on; grid on; title('Channel 6 Real'); xlabel('time');
    % subplot(8,2,12); plot((inputOffset+1:inputOffset+NFineOut),imag(output(inputOffset+1:inputOffset+NFineOut,6))); box on; grid on; title('Channel 6 Imag'); xlabel('time');
    % subplot(8,2,13); plot((inputOffset+1:inputOffset+NFineOut),real(output(inputOffset+1:inputOffset+NFineOut,7))); box on; grid on; title('Channel 7 Real'); xlabel('time');
    % subplot(8,2,14); plot((inputOffset+1:inputOffset+NFineOut),imag(output(inputOffset+1:inputOffset+NFineOut,7))); box on; grid on; title('Channel 7 Imag'); xlabel('time');
    % subplot(8,2,15); plot((inputOffset+1:inputOffset+NFineOut),real(output(inputOffset+1:inputOffset+NFineOut,8))); box on; grid on; title('Channel 8 Real'); xlabel('time');
    % subplot(8,2,16); plot((inputOffset+1:inputOffset+NFineOut),imag(output(inputOffset+1:inputOffset+NFineOut,8))); box on; grid on; title('Channel 8 Imag'); xlabel('time');

    disp('calculating fine channel spectra');

    % analyse middle half: length = fftLength
    samples = output(inputOffset+1:inputOffset+fftLength,:);
    spectra = fft(samples);
    spectra = fftshift(spectra,1);
    spectra = fftshift(spectra,2);

    % figure;
    % subplot(8,1,1); plot((1:fftLength),abs(spectra(1:fftLength,1))); box on; grid on; title('Channel 1 Mag'); xlabel('frequency'); xlim([1 fftLength]);
    % subplot(8,1,2); plot((1:fftLength),abs(spectra(1:fftLength,2))); box on; grid on; title('Channel 2 Mag'); xlabel('frequency'); xlim([1 fftLength]);
    % subplot(8,1,3); plot((1:fftLength),abs(spectra(1:fftLength,3))); box on; grid on; title('Channel 3 Mag'); xlabel('frequency'); xlim([1 fftLength]);
    % subplot(8,1,4); plot((1:fftLength),abs(spectra(1:fftLength,4))); box on; grid on; title('Channel 4 Mag'); xlabel('frequency'); xlim([1 fftLength]);
    % subplot(8,1,5); plot((1:fftLength),abs(spectra(1:fftLength,5))); box on; grid on; title('Channel 5 Mag'); xlabel('frequency'); xlim([1 fftLength]);
    % subplot(8,1,6); plot((1:fftLength),abs(spectra(1:fftLength,6))); box on; grid on; title('Channel 6 Mag'); xlabel('frequency'); xlim([1 fftLength]);
    % subplot(8,1,7); plot((1:fftLength),abs(spectra(1:fftLength,7))); box on; grid on; title('Channel 7 Mag'); xlabel('frequency'); xlim([1 fftLength]);
    % subplot(8,1,8); plot((1:fftLength),abs(spectra(1:fftLength,8))); box on; grid on; title('Channel 8 Mag'); xlabel('frequency'); xlim([1 fftLength]);


    %% DISCARD OVERSAMPLED PORTIONS & OPTIONALLY EQUALISE RIPPLE

    fprintf('discarding oversampled portions');
    if (equaliseRipple)
        fprintf(' and equalising ripple');
    end
    fprintf('\n');
    FN = complex(zeros(fftLength*OS_De/OS_Nu,N));
    for chan = 1:N
        % Keep only the pass-band
        discard = (1.0 - (OS_De/OS_Nu))/2.0;
        FN(:,chan) = spectra(round(discard*fftLength)+1:round((1.0-discard)*fftLength),chan);

        if(equaliseRipple)
            % load PFB prototype filter transfer function
            load('TF_points.mat');

            % use just the baseband passband section of transfer function
            % - apply to both halves of channel
            passband_len = (fftLength/2)*OS_De/OS_Nu;
            for ii = 1:passband_len
                FN(ii,chan) = FN(ii,chan)/abs(H0(passband_len-ii+2));
                FN(passband_len+ii,chan) = FN(passband_len+ii,chan)/abs(H0(ii));
            end;
        end;
    end;


    %% Combine chunks & back-transform

    fprintf('combining channels and back transforming\n');
    FN_width = fftLength*OS_De/OS_Nu;
    FFFF = FN(FN_width/2+1:FN_width,1); % upper half of chan 1 is first part of FFFF
    for chan = 2 : N
        FFFF = [FFFF; FN(:,chan)];
    end;
    FFFF = [FFFF; FN(1:FN_width/2,1)]; % lower half of chan 1 is last part of FFFF

    len = length(FFFF);
    %save('N_channels','FFFF');

    % back transform
    z1 = ifft(fftshift(FFFF))./(OS_Nu/OS_De);  % re-scale by OS factor
    
    % shift centre to common location for all impulseOffset cases
    z1 = circshift(z1, -impulseOffset);

%     figure;
%     subplot(211); plot((1:len),abs(FFFF)); box on; grid on; title('FFFF Mag'); 
%     subplot(212); plot((1:len),angle(FFFF)); box on; grid on; title('FFFF Phase'); xlabel('time');

%     figure;
%     subplot(211); plot((1:len),real(z1(1:len))); box on; grid on; title('z1 Real'); 
%     subplot(212); plot((1:len),imag(z1(1:len))); box on; grid on; title('z1 Imag'); xlabel('time');

%     figure;
%     plot((1:len),20.0*log10(abs(z1(1:len))+0.00000000000001)); box on; grid on; title('Output Power By Sample'); xlabel('sample'); ylabel('Power (dB)');
%     pause

    output_offset =  (OS_De*N/OS_Nu)*(fftLength-inputOffset) - (Ntaps-1)/2;

    amplitudeTF(impulseOffset,1) = abs(z1(output_offset,1));
    %fprintf('\nAmplitude TF = %f (%3.2f dB)\n',amplitudeTF(impulseOffset,1),10.0*log10(amplitudeTF(impulseOffset,1)));
    
    z1(output_offset,1) = 0.0;
    maxSpuriousPower(impulseOffset,1) = (max(abs(z1)))^2;
%     fprintf('\nMax Spurious Power = %e (%3.2f dB)\n',maxSpuriousPower,10.0*log10(maxSpuriousPower));
      
    for ii = 1:blockLength
      totalSpuriousPower(impulseOffset,1) = totalSpuriousPower(impulseOffset,1) + (abs(z1(ii,1)))^2;
    end;

    averageSpuriousPower = averageSpuriousPower + abs(z1).^2;
        

end;

averageSpuriousPower = averageSpuriousPower./PFBarms;  % averaged over number of impulses tested

fprintf('\nMax Average Spurious Power = %3.2f dB\n',10.0*log10(max(averageSpuriousPower)));
fprintf('\nSNR = %3.2f dB\n',10.0*log10(1.0/sum(averageSpuriousPower)));  % total signal power = 1 (per impulse)

figure;
plot((1:PFBarms),amplitudeTF(1:PFBarms,1)); box on; grid on; title('Amplitude vs Impulse Location'); xlabel('Impulse Location'); ylabel('Output Amplitude');

powerTFdB = 20.0*log10(amplitudeTF);  % abs value before this
figure;
plot((1:PFBarms),powerTFdB(1:PFBarms,1)); box on; grid on; title('Normalised Power vs Impulse Location'); xlabel('Impulse Location'); ylabel('Output Power (dB)');

maxSpuriousPower = 10.0*log10(maxSpuriousPower);  % power value before this
figure;
plot((1:PFBarms),maxSpuriousPower(1:PFBarms,1)); box on; grid on; title('Max Spurious Power'); xlabel('Impulse Location'); ylabel('Max Spurious Power Component (dB)');

totalSpuriousPower = 10.0*log10(totalSpuriousPower);  % power value before this
figure;
plot((1:PFBarms),totalSpuriousPower(1:PFBarms,1)); box on; grid on; title('Total Spurious Power'); xlabel('Impulse Location'); ylabel('Total Spurious Power (dB)');

averageSpuriousPower = 10.0*log10(averageSpuriousPower);  % power value before this
figure;
plot((1:blockLength),averageSpuriousPower(1:blockLength,1)); box on; grid on; title('Average Spurious Power'); xlabel('Time Sample'); ylabel('Average Spurious Power (dB)');

fprintf('\n\nDone! Press any key to close plots and exit...\n\n');
pause;
close all;
clear all;
