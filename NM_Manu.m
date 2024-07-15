% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1kW NM Manu%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Note: Import 'Output Type' as Numeric Matrix
% t2=aprEXPMATLABS2(:,1); %time
% x_t2=aprEXPMATLABS2(:,2); %datapoints
% % figure(11)
% % plot(t2,x_t2)
% % title('1kW NM Manu Raw Data vs. Time')
% 
% L2=length(x_t2); %number of samples (samples)
% Fs2=1/mean(diff(t2)); %sampling frequencyc(samples/sec)
% f2=(0:L2-1) * Fs2/L2; %frequency (Hz, 1/s) [samples/sec/samples]
% disp(Fs2)
% disp(f2)
% xlabel('relativetime')
% ylabel('current')
% 
% 
% %FFT
% X2=fft(x_t2);
% X_norm2=1/L2 * X2;
% % figure(12)
% % plot(real(X2))
% % title('1kW NM Manu FFT')
% 
% %Amplitude and Phase Spectrum
% figure(13)
% stem(f2,abs(X_norm2),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel('Amplitude')
% title('1kW NM Manu FFT Amp')
% 
% % figure(14)
% % stem(f2,angle(X_norm2),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel ('Phase Angle')
% % title('1kW NM Manu FFT Phase')
% 
% % %table of values (Coeffs, Freq, Amplitude, Angle)
% % table3 = table(X_norm2,f2',abs(X_norm2),angle(X_norm2));
% % table3.Properties.VariableNames = {'FFTCoeffs','Frequency','Amplitude','Phase'};
% % disp(table3)
% 
% %Reconstruct Using k Dominant frequencies
% X_recon2=X2;
% k=30; %keeps first k and last k + the mean
% for i=(k+10):(L2-k) %starting with the 40ith term and everything in between =0
%     X_recon2(i)=0;
% end 
% 
% %IFFT
% x_recon2= ifft(X_recon2);
% figure(15)
% plot(t2,x_t2,'k-*','LineWidth',2) %raw data
% hold on, grid on
% plot(t2,x_recon2,'b--','LineWidth',2) %reconstructed
% xlabel('Time(s)'), ylabel('Current (Amps)')
% title('1kW NM Manu Current Height vs. Time')
% 12 - append(num2str(k),'FFT Dominant Terms');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%30kW NM Manu%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Note: Import 'Output Type' as Numeric Matrix
% t6=aprEXPMATLABS6(:,1); %time
% x_t6=aprEXPMATLABS6(:,2); %datapoints
% % figure(31)
% % plot(t6,x_t6)
% % title('30kW NM Manu Raw Data vs. Time')
% 
% L6=length(x_t6); %number of samples (samples)
% Fs6=1/mean(diff(t6)); %sampling frequencyc(samples/sec)
% f6=(0:L6-1) * Fs6/L6; %frequency (Hz, 1/s) [samples/sec/samples]
% disp(Fs6)
% disp(f6)
% xlabel('relativetime')
% ylabel('current')
% 
% 
% %FFT
% X6=fft(x_t6);
% X_norm6=1/L6 * X6;
% % figure(32)
% % plot(real(X6))
% % title('30kW NM Manu FFT')
% 
% %Amplitude and Phase Spectrum
% figure(33)
% stem(f6,abs(X_norm6),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel('Amplitude')
% title('30kW NM Manu FFT Amp')
% 
% % figure(34)
% % stem(f6,angle(X_norm6),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel ('Phase Angle')
% % title('30kW NM Manu FFT Phase')
% 
% % %table of values (Coeffs, Freq, Amplitude, Angle)
% % table7 = table(X_norm6,f6',abs(X_norm6),angle(X_norm6));
% % table7.Properties.VariableNames = {'FFTCoeffs','Frequency','Amplitude','Phase'};
% % disp(table7)
% % 
% %Reconstruct Using k Dominant frequencies
% X_recon6=X6;
% k=30; %keeps first k and last k + the mean
% for i=(k+10):(L6-k) %starting with the 40ith term and everything in between =0
%     X_recon6(i)=0;
% end 
% 
% %IFFT
% x_recon6= ifft(X_recon6);
% figure(35)
% plot(t6,x_t6,'k-*','LineWidth',2) %raw data
% hold on, grid on
% plot(t6,x_recon6,'b--','LineWidth',2) %reconstructed
% xlabel('Time(s)'), ylabel('Current (Amps)')
% title('30kW NM Manu Current Height vs. Time')
% 12 - append(num2str(k),'FFT Dominant Terms');
% 
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%500kW NM Manu%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Note: Import 'Output Type' as Numeric Matrix
% t10=aprEXPMATLABS10(:,1); %time
% x_t10=aprEXPMATLABS10(:,2); %datapoints
% % figure(51)
% % plot(t10,x_t10)
% % title('500kW NM Manu Raw Data vs. Time')
% 
% L10=length(x_t10); %number of samples (samples)
% Fs10=1/mean(diff(t10)); %sampling frequencyc(samples/sec)
% f10=(0:L10-1) * Fs10/L10; %frequency (Hz, 1/s) [samples/sec/samples]
% disp(Fs10)
% disp(f10)
% xlabel('relativetime')
% ylabel('current')
% 
% 
% %FFT
% X10=fft(x_t10);
% X_norm10=1/L10 * X10;
% % figure(52)
% % plot(real(X10))
% % title('500kW NM Manu FFT')
% 
% %Amplitude and Phase Spectrum
% figure(53)
% stem(f10,abs(X_norm10),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel('Amplitude')
% title('500kW NM Manu FFT Amp')
% 
% % figure(54)
% % stem(f10,angle(X_norm10),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel ('Phase Angle')
% % title('500kW NM Manu FFT Phase')
% 
% % %table of values (Coeffs, Freq, Amplitude, Angle)
% % table11 = table(X_norm10,f10',abs(X_norm10),angle(X_norm10));
% % table11.Properties.VariableNames = {'FFTCoeffs','Frequency','Amplitude','Phase'};
% % disp(table11)
% % 
% %Reconstruct Using k Dominant frequencies
% X_recon10=X10;
% k=30; %keeps first k and last k + the mean
% for i=(k+10):(L10-k) %starting with the 40ith term and everything in between =0
%     X_recon10(i)=0;
% end 
% 
% %IFFT
% x_recon10= ifft(X_recon10);
% figure(55)
% plot(t10,x_t10,'k-*','LineWidth',2) %raw data
% hold on, grid on
% plot(t10,x_recon10,'b--','LineWidth',2) %reconstructed
% xlabel('Time(s)'), ylabel('Current (Amps)')
% title('500kW NM Manu Current Height vs. Time')
% 12 - append(num2str(k),'FFT Dominant Terms');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%900kW NM Manu%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: Import 'Output Type' as Numeric Matrix
t14=aprEXPMATLABS14(:,1); %time
x_t14=aprEXPMATLABS14(:,2); %datapoints
% figure(71)
% plot(t14,x_t14)
% title('900kW NM Manu Raw Data vs. Time')

L14=length(x_t14); %number of samples (samples)
Fs14=1/mean(diff(t14)); %sampling frequencyc(samples/sec)
f14=(0:L14-1) * Fs14/L14; %frequency (Hz, 1/s) [samples/sec/samples]
disp(Fs14)
disp(f14)
xlabel('relativetime')
ylabel('current')


%FFT
X14=fft(x_t14);
X_norm14=1/L14 * X14;
% figure(72)
% plot(real(X14))
% title('900kW NM Manu FFT')

%Amplitude and Phase Spectrum
figure(73)
stem(f14,abs(X_norm14),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel('Amplitude')
title('900kW NM Manu FFT Amp')

% figure(74)
% stem(f14,angle(X_norm14),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel ('Phase Angle')
% title('900kW NM Manu FFT Phase')

% %table of values (Coeffs, Freq, Amplitude, Angle)
% table15 = table(X_norm14,f14',abs(X_norm14),angle(X_norm14));
% table15.Properties.VariableNames = {'FFTCoeffs','Frequency','Amplitude','Phase'};
% disp(table15)
% 
%Reconstruct Using k Dominant frequencies
X_recon14=X14;
k=30; %keeps first k and last k + the mean
for i=(k+10):(L14-k) %starting with the 40ith term and everything in between =0
    X_recon14(i)=0;
end 

%IFFT
x_recon14= ifft(X_recon14);
figure(75)
plot(t14,x_t14,'k-*','LineWidth',2) %raw data
hold on, grid on
plot(t14,x_recon14,'b--','LineWidth',2) %reconstructed
xlabel('Time(s)'), ylabel('Current (Amps)')
title('900kW NM Manu Current Height vs. Time')
12 - append(num2str(k),'FFT Dominant Terms');
% 
%