%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1kW FC Manu%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: Import 'Output Type' as Numeric Matrix
t=aprEXPMATLAB(:,1); %time
x_t=aprEXPMATLAB(:,2); %datapoints
% figure(1)
% plot(t,x_t)
% title('1kW FC Manu Raw Data vs. Time')

L=length(x_t); %number of samples (samples)
Fs=1/mean(diff(t)); %sampling frequencyc(samples/sec)
f=(0:L-1) * Fs/L; %frequency (Hz, 1/s) [samples/sec/samples]
disp(Fs)
disp(f)
xlabel('relativetime')
ylabel('current')


%FFT
X=fft(x_t);
X_norm=1/L * X;
% figure(2)
% plot(real(X))
% title('1kW FC Manu FFT')

%Amplitude and Phase Spectrum
figure(3)
stem(f,abs(X_norm),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel('Amplitude')
title('1kW FC Manu FFT Amp')

% figure(4)
% stem(f,angle(X_norm),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel ('Phase Angle')
% title('1kW FC Manu FFT Phase')
% 
% %table of values (Coeffs, Freq, Amplitude, Angle)
% table1 = table(X_norm,f',abs(X_norm),angle(X_norm));
% table1.Properties.VariableNames = {'FFTCoeffs','Frequency','Amplitude','Phase'};
% disp(table1)
% 
%Reconstruct Using k Dominant frequencies
X_recon=X;
k=30; %keeps first k and last k + the mean
for i=(k+10):(L-k) %starting with the 40ith term and everything in between =0
    X_recon(i)=0;
end 

%IFFT
x_recon= ifft(X_recon);
figure(5)
plot(t,x_t,'k-*','LineWidth',2) %raw data
hold on, grid on
plot(t,x_recon,'b--','LineWidth',2) %reconstructed
xlabel('Time(s)'), ylabel('Current (Amps)')
title('1kW FC Manu Current Height vs. Time')
12 - append(num2str(k),'FFT Dominant Terms');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%30kW FC Manu%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: Import 'Output Type' as Numeric Matrix
t4=aprEXPMATLABS4(:,1); %time
x_t4=aprEXPMATLABS4(:,2); %datapoints
% figure(21)
% plot(t4,x_t4)
% title('30kW FC Manu Raw Data vs. Time')

L4=length(x_t4); %number of samples (samples)
Fs4=1/mean(diff(t4)); %sampling frequencyc(samples/sec)
f4=(0:L4-1) * Fs4/L4; %frequency (Hz, 1/s) [samples/sec/samples]
disp(Fs4)
disp(f4)
xlabel('relativetime')
ylabel('current')


%FFT
X4=fft(x_t4);
X_norm4=1/L4 * X4;
% figure(22)
% plot(real(X4))
% title('30kW FC Manu FFT')

%Amplitude and Phase Spectrum
figure(23)
stem(f4,abs(X_norm4),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel('Amplitude')
title('30kW FC Manu FFT Amp')

% figure(24)
% stem(f4,angle(X_norm4),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel ('Phase Angle')
% title('30kW FC Manu FFT Phase')

% %table of values (Coeffs, Freq, Amplitude, Angle)
% table5 = table(X_norm4,f4',abs(X_norm4),angle(X_norm4));
% table5.Properties.VariableNames = {'FFTCoeffs','Frequency','Amplitude','Phase'};
% disp(table5)
% 
%Reconstruct Using k Dominant frequencies
X_recon4=X4;
k=30; %keeps first k and last k + the mean
for i=(k+10):(L4-k) %starting with the 40ith term and everything in between =0
    X_recon4(i)=0;
end 

%IFFT
x_recon4= ifft(X_recon4);
figure(25)
plot(t4,x_t4,'k-*','LineWidth',2) %raw data
hold on, grid on
plot(t4,x_recon4,'b--','LineWidth',2) %reconstructed
xlabel('Time(s)'), ylabel('Current (Amps)')
title('30kW FC Manu Current Height vs. Time')
12 - append(num2str(k),'FFT Dominant Terms');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%500kW FC Manu%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: Import 'Output Type' as Numeric Matrix
t8=aprEXPMATLABS8(:,1); %time
x_t8=aprEXPMATLABS8(:,2); %datapoints
% figure(41)
% plot(t8,x_t8)
% title('500kW FC Manu Raw Data vs. Time')

L8=length(x_t8); %number of samples (samples)
Fs8=1/mean(diff(t8)); %sampling frequencyc(samples/sec)
f8=(0:L8-1) * Fs8/L8; %frequency (Hz, 1/s) [samples/sec/samples]
disp(Fs8)
disp(f8)
xlabel('relativetime')
ylabel('current')


%FFT
X8=fft(x_t8);
X_norm8=1/L8 * X8;
% figure(42)
% plot(real(X8))
% title('500 kW FC Manu FFT')

%Amplitude and Phase Spectrum
figure(43)
stem(f8,abs(X_norm8),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel('Amplitude')
title('500kW FC Manu FFT Amp')

% figure(44)
% stem(f8,angle(X_norm8),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel ('Phase Angle')
% title('500kW FC Manu FFT Phase')

% %table of values (Coeffs, Freq, Amplitude, Angle)
% table9 = table(X_norm8,f8',abs(X_norm8),angle(X_norm8));
% table9.Properties.VariableNames = {'FFTCoeffs','Frequency','Amplitude','Phase'};
% disp(table9)
% 
%Reconstruct Using k Dominant frequencies
X_recon8=X8;
k=30; %keeps first k and last k + the mean
for i=(k+10):(L8-k) %starting with the 40ith term and everything in between =0
    X_recon8(i)=0;
end 

%IFFT
x_recon8= ifft(X_recon8);
figure(45)
plot(t8,x_t8,'k-*','LineWidth',2) %raw data
hold on, grid on
plot(t8,x_recon8,'b--','LineWidth',2) %reconstructed
xlabel('Time(s)'), ylabel('Current (Amps)')
title('500kW FC Manu Current Height vs. Time')
12 - append(num2str(k),'FFT Dominant Terms');

% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%900kW FC Manu%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: Import 'Output Type' as Numeric Matrix
t12=aprEXPMATLABS12(:,1); %time
x_t12=aprEXPMATLABS12(:,2); %datapoints
% figure(61)
% plot(t12,x_t12)
% title('900kW FC Manu Raw Data vs. Time')

L12=length(x_t12); %number of samples (samples)
Fs12=1/mean(diff(t12)); %sampling frequencyc(samples/sec)
f12=(0:L12-1) * Fs12/L12; %frequency (Hz, 1/s) [samples/sec/samples]
disp(Fs12)
disp(f12)
xlabel('relativetime')
ylabel('current')


%FFT
X12=fft(x_t12);
X_norm12=1/L12 * X12;
% figure(62)
% plot(real(X12))
% title('900kW FC Manu FFT')

%Amplitude and Phase Spectrum
figure(63)
stem(f12,abs(X_norm12),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel('Amplitude')
title('900kW FC Manu FFT Amp')

% figure(64)
% stem(f12,angle(X_norm12),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel ('Phase Angle')
% title('900kW FC Manu FFT Phase')

% %table of values (Coeffs, Freq, Amplitude, Angle)
% table13 = table(X_norm12,f12',abs(X_norm12),angle(X_norm12));
% table13.Properties.VariableNames = {'FFTCoeffs','Frequency','Amplitude','Phase'};
% disp(table13)
% 
%Reconstruct Using k Dominant frequencies
X_recon12=X12;
k=30; %keeps first k and last k + the mean
for i=(k+10):(L12-k) %starting with the 40ith term and everything in between =0
    X_recon12(i)=0;
end 

%IFFT
x_recon12= ifft(X_recon12);
figure(65)
plot(t12,x_t12,'k-*','LineWidth',2) %raw data
hold on, grid on
plot(t12,x_recon12,'b--','LineWidth',2) %reconstructed
xlabel('Time(s)'), ylabel('Current (Amps)')
title('900kW FC Manu Current Height vs. Time')
12 - append(num2str(k),'FFT Dominant Terms');