%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1kW FC Auto%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: Import 'Output Type' as Numeric Matrix
t1=autodiffofftransientS1(:,1); %time
x_t1=autodiffofftransientS1(:,2); %datapoints
% figure(6)
% plot(t1,x_t1)
% title('1kW FC Auto Raw Data vs. Time')

L1=length(x_t1); %number of samples (samples)
Fs1=1/mean(diff(t1)); %sampling frequencyc(samples/sec)
f1=(0:L1-1) * Fs1/L1; %frequency (Hz, 1/s) [samples/sec/samples]
disp(Fs1)
disp(f1)
xlabel('relativetime')
ylabel('current')

%FFT
X1=fft(x_t1);
X_norm1=1/L1 * X1;
% figure(7)
% plot(real(X1))
% title('1kW FC Auto FFT')


%Amplitude and Phase Spectrum
figure(8)
stem(f1,abs(X_norm1),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel('Amplitude')
title('1kW FC Auto FFT Amp')

% figure(9)
% stem(f1,angle(X_norm1),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel ('Phase Angle')
% title('1kW FC Auto FFT Phase')

% %table of values (Coeffs, Freq, Amplitude, Angle)
% table2 = table(X_norm1,f1',abs(X_norm1),angle(X_norm1));
% table2.Properties.VariableNames = {'FFTCoeffs','Frequency','Amplitude','Phase'};
% disp(table2)
% 
%Reconstruct Using k Dominant frequencies
X_recon1=X1;
k=30; %keeps first k and last k + the mean
for i=(k+10):(L1-k) %starting with the 40ith term and everything in between =0
    X_recon1(i)=0;
end 

%IFFT
x_recon1= ifft(X_recon1);
figure(10)
plot(t1,x_t1,'k-*','LineWidth',2) %raw data
hold on, grid on
plot(t1,x_recon1,'b--','LineWidth',2) %reconstructed
xlabel('Time(s)'), ylabel('Current (Amps)')
title('1kW FC Auto Current Height vs. Time')
12 - append(num2str(k),'FFT Dominant Terms');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%30kW FC Auto%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: Import 'Output Type' as Numeric Matrix
t5=autodiffofftransientS2(:,1); %time
x_t5=autodiffofftransientS2(:,2); %datapoints
% figure(26)
% plot(t5,x_t5)
% title('30kW FC Auto Raw Data vs. Time')

L5=length(x_t5); %number of samples (samples)
Fs5=1/mean(diff(t5)); %sampling frequencyc(samples/sec)
f5=(0:L5-1) * Fs5/L5; %frequency (Hz, 1/s) [samples/sec/samples]
disp(Fs5)
disp(f5)
xlabel('relativetime')
ylabel('current')

%FFT
X5=fft(x_t5);
X_norm5=1/L5 * X5;
% figure(27)
% plot(real(X5))
% title('30kW FC Auto FFT')


%Amplitude and Phase Spectrum
figure(28)
stem(f5,abs(X_norm5),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel('Amplitude')
title('30kW FC Auto FFT Amp')

% figure(29)
% stem(f5,angle(X_norm5),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel ('Phase Angle')
% title('30kW FC Auto FFT Phase')

% %table of values (Coeffs, Freq, Amplitude, Angle)
% table6 = table(X_norm5,f5',abs(X_norm5),angle(X_norm5));
% table6.Properties.VariableNames = {'FFTCoeffs','Frequency','Amplitude','Phase'};
% disp(table6)
% 
%Reconstruct Using k Dominant frequencies
X_recon5=X5;
k=30; %keeps first k and last k + the mean
for i=(k+10):(L5-k) %starting with the 40ith term and everything in between =0
    X_recon5(i)=0;
end 

%IFFT
x_recon5= ifft(X_recon5);
figure(30)
plot(t5,x_t5,'k-*','LineWidth',2) %raw data
hold on, grid on
plot(t5,x_recon5,'b--','LineWidth',2) %reconstructed
xlabel('Time(s)'), ylabel('Current (Amps)')
title('30kW FC Auto Current Height vs. Time')
12 - append(num2str(k),'FFT Dominant Terms');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%500kW FC Auto%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Note: Import 'Output Type' as Numeric Matrix
t9=autodiffofftransientS3(:,1); %time
x_t9=autodiffofftransientS3(:,2); %datapoints
% figure(46)
% plot(t9,x_t9)
% title('500kW FC Auto Raw Data vs. Time')

L9=length(x_t9); %number of samples (samples)
Fs9=1/mean(diff(t9)); %sampling frequencyc(samples/sec)
f9=(0:L9-1) * Fs9/L9; %frequency (Hz, 1/s) [samples/sec/samples]
disp(Fs9)
disp(f9)
xlabel('relativetime')
ylabel('current')

%FFT
X9=fft(x_t9);
X_norm9=1/L9 * X9;
% figure(47)
% plot(real(X9))
% title('500kW FC Auto FFT')


%Amplitude and Phase Spectrum
figure(48)
stem(f9,abs(X_norm9),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel('Amplitude')
title('500kW FC Auto FFT Amp')

% figure(49)
% stem(f9,angle(X_norm9),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel ('Phase Angle')
% title('500kW FC Auto FFT Phase')

% %table of values (Coeffs, Freq, Amplitude, Angle)
% table10 = table(X_norm9,f9',abs(X_norm9),angle(X_norm9));
% table10.Properties.VariableNames = {'FFTCoeffs','Frequency','Amplitude','Phase'};
% disp(table10)
% 
%Reconstruct Using k Dominant frequencies
X_recon9=X9;
k=30; %keeps first k and last k + the mean
for i=(k+10):(L9-k) %starting with the 40ith term and everything in between =0
    X_recon9(i)=0;
end 

%IFFT
x_recon9= ifft(X_recon9);
figure(50)
plot(t9,x_t9,'k-*','LineWidth',2) %raw data
hold on, grid on
plot(t9,x_recon9,'b--','LineWidth',2) %reconstructed
xlabel('Time(s)'), ylabel('Current (Amps)')
title('500kW FC Auto Current Height vs. Time')
12 - append(num2str(k),'FFT Dominant Terms');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%900kW FC Auto%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: Import 'Output Type' as Numeric Matrix
t13=autodiffofftransientS4(:,1); %time
x_t13=autodiffofftransientS4(:,2); %datapoints
% figure(66)
% plot(t13,x_t13)
% title('900kW FC Auto Raw Data vs. Time')

L13=length(x_t13); %number of samples (samples)
Fs13=1/mean(diff(t13)); %sampling frequencyc(samples/sec)
f13=(0:L13-1) * Fs13/L13; %frequency (Hz, 1/s) [samples/sec/samples]
disp(Fs13)
disp(f13)
xlabel('relativetime')
ylabel('current')

%FFT
X13=fft(x_t13);
X_norm13=1/L13 * X13;
% figure(67)
% plot(real(X13))
% title('900kW FC Auto FFT')


%Amplitude and Phase Spectrum
figure(68)
stem(f13,abs(X_norm13),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel('Amplitude')
title('900kW FC Auto FFT Amp')

% figure(69) %Nice!
% stem(f13,angle(X_norm13),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel ('Phase Angle')
% title('900kW FC Auto FFT Phase')

% %table of values (Coeffs, Freq, Amplitude, Angle)
% table14 = table(X_norm13,f13',abs(X_norm13),angle(X_norm13));
% table14.Properties.VariableNames = {'FFTCoeffs','Frequency','Amplitude','Phase'};
% disp(table14)
% 
%Reconstruct Using k Dominant frequencies
X_recon13=X13;
k=30; %keeps first k and last k + the mean
for i=(k+10):(L13-k) %starting with the 40ith term and everything in between =0
    X_recon13(i)=0;
end 

%IFFT
x_recon13= ifft(X_recon13);
figure(70)
plot(t13,x_t13,'k-*','LineWidth',2) %raw data
hold on, grid on
plot(t13,x_recon13,'b--','LineWidth',2) %reconstructed
xlabel('Time(s)'), ylabel('Current (Amps)')
title('900kW FC Auto Current Height vs. Time')
12 - append(num2str(k),'FFT Dominant Terms');

% 