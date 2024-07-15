%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1kW NM Auto%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: Import 'Output Type' as Numeric Matrix
t3=nmwdiffoffautoS1(:,1); %time
x_t3=nmwdiffoffautoS1(:,2); %datapoints
% figure(16)
% plot(t3,x_t3)
% title('1kW NM Auto Raw Data vs. Time')

L3=length(x_t3); %number of samples (samples)
Fs3=1/mean(diff(t3)); %sampling frequencyc(samples/sec)
f3=(0:L3-1) * Fs3/L3; %frequency (Hz, 1/s) [samples/sec/samples]
disp(Fs3)
disp(f3)
xlabel('relativetime')
ylabel('current')

%FFT
X3=fft(x_t3);
X_norm3=1/L3 * X3;
% figure(17)
% plot(real(X3))
% title('1kW NM Auto FFT')


%Amplitude and Phase Spectrum
figure(18)
stem(f3,abs(X_norm3),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel('Amplitude')
title('1kW NM Auto FFT Amp')

% figure(19)
% stem(f3,angle(X_norm3),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel ('Phase Angle')
% title('1kW NM Auto FFT Phase')

% %table of values (Coeffs, Freq, Amplitude, Angle)
% table4 = table(X_norm3,f3',abs(X_norm3),angle(X_norm3));
% table4.Properties.VariableNames = {'FFTCoeffs','Frequency','Amplitude','Phase'};
% disp(table4)
% 
%Reconstruct Using k Dominant frequencies
X_recon3=X3;
k=30; %keeps first k and last k + the mean
for i=(k+10):(L3-k) %starting with the 40ith term and everything in between =0
    X_recon3(i)=0;
end 

%IFFT
x_recon3= ifft(X_recon3);
figure(20)
plot(t3,x_t3,'k-*','LineWidth',2) %raw data
hold on, grid on
plot(t3,x_recon3,'b--','LineWidth',2) %reconstructed
xlabel('Time(s)'), ylabel('Current (Amps)')
title('1kW NM Auto Current Height vs. Time')
12 - append(num2str(k),'FFT Dominant Terms');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%30kW NM Auto%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Note: Import 'Output Type' as Numeric Matrix
% t7=nmwdiffoffautoS2(:,1); %time
% x_t7=nmwdiffoffautoS2(:,2); %datapoints
% % figure(36)
% % plot(t7,x_t7)
% % title('30kW NM Auto Raw Data vs. Time')
% 
% L7=length(x_t7); %number of samples (samples)
% Fs7=1/mean(diff(t7)); %sampling frequencyc(samples/sec)
% f7=(0:L7-1) * Fs7/L7; %frequency (Hz, 1/s) [samples/sec/samples]
% disp(Fs7)
% disp(f7)
% xlabel('relativetime')
% ylabel('current')
% 
% %FFT
% X7=fft(x_t7);
% X_norm7=1/L7 * X7;
% % figure(37)
% % plot(real(X7))
% % title('30kW NM Auto FFT')
% 
% 
% %Amplitude and Phase Spectrum
% figure(38)
% stem(f7,abs(X_norm7),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel('Amplitude')
% title('30kW NM Auto FFT Amp')
% 
% % figure(39)
% % stem(f7,angle(X_norm7),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel ('Phase Angle')
% % title('30kW NM Auto FFT Phase')
% 
% % %table of values (Coeffs, Freq, Amplitude, Angle)
% % table8 = table(X_norm7,f7',abs(X_norm7),angle(X_norm7));
% % table8.Properties.VariableNames = {'FFTCoeffs','Frequency','Amplitude','Phase'};
% % disp(table8)
% % 
% %Reconstruct Using k Dominant frequencies
% X_recon7=X7;
% k=30; %keeps first k and last k + the mean
% for i=(k+10):(L7-k) %starting with the 40ith term and everything in between =0
%     X_recon7(i)=0;
% end 
% 
% %IFFT
% x_recon7= ifft(X_recon7);
% figure(40)
% plot(t7,x_t7,'k-*','LineWidth',2) %raw data
% hold on, grid on
% plot(t7,x_recon7,'b--','LineWidth',2) %reconstructed
% xlabel('Time(s)'), ylabel('Current (Amps)')
% title('30kW NM Auto Current Height vs. Time')
% 12 - append(num2str(k),'FFT Dominant Terms');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%500kW NM Auto%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Note: Import 'Output Type' as Numeric Matrix
% t11=nmwdiffoffautoS3(:,1); %time
% x_t11=nmwdiffoffautoS3(:,2); %datapoints
% % figure(56)
% % plot(t11,x_t11)
% % title('500kW NM Auto Raw Data vs. Time')
% 
% L11=length(x_t11); %number of samples (samples)
% Fs11=1/mean(diff(t11)); %sampling frequencyc(samples/sec)
% f11=(0:L11-1) * Fs11/L11; %frequency (Hz, 1/s) [samples/sec/samples]
% disp(Fs11)
% disp(f11)
% xlabel('relativetime')
% ylabel('current')
% 
% %FFT
% X11=fft(x_t11);
% X_norm11=1/L11 * X11;
% % figure(57)
% % plot(real(X11))
% % title('500kW NM Auto FFT')
% 
% 
% %Amplitude and Phase Spectrum
% figure(58)
% stem(f11,abs(X_norm11),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel('Amplitude')
% title('500kW NM Auto FFT Amp')
% 
% % figure(59)
% % stem(f11,angle(X_norm11),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel ('Phase Angle')
% % title('500kW NM Auto FFT Phase')
% 
% % %table of values (Coeffs, Freq, Amplitude, Angle)
% % table12 = table(X_norm11,f11',abs(X_norm11),angle(X_norm11));
% % table12.Properties.VariableNames = {'FFTCoeffs','Frequency','Amplitude','Phase'};
% % disp(table12)
% % 
% %Reconstruct Using k Dominant frequencies
% X_recon11=X11;
% k=30; %keeps first k and last k + the mean
% for i=(k+10):(L11-k) %starting with the 40ith term and everything in between =0
%     X_recon11(i)=0;
% end 
% 
% %IFFT
% x_recon11= ifft(X_recon11);
% figure(60)
% plot(t11,x_t11,'k-*','LineWidth',2) %raw data
% hold on, grid on
% plot(t11,x_recon11,'b--','LineWidth',2) %reconstructed
% xlabel('Time(s)'), ylabel('Current (Amps)')
% title('500kW NM Auto Current Height vs. Time')
% 12 - append(num2str(k),'FFT Dominant Terms');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%900kW NM Auto%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %Note: Import 'Output Type' as Numeric Matrix
% t15=nmwdiffoffautoS4(:,1); %time
% x_t15=nmwdiffoffautoS4(:,2); %datapoints
% figure(76)
% plot(t15,x_t15)
% title('900kW NM Auto Raw Data vs. Time')
% 
% L15=length(x_t15); %number of samples (samples)
% Fs15=1/mean(diff(t15)); %sampling frequencyc(samples/sec)
% f15=(0:L15-1) * Fs15/L15; %frequency (Hz, 1/s) [samples/sec/samples]
% disp(Fs15)
% disp(f15)
% xlabel('relativetime')
% ylabel('current')
% 
% %FFT
% X15=fft(x_t15);
% X_norm15=1/L15 * X15;
% % figure(77)
% % plot(real(X15))
% % title('900kW NM Auto FFT')
% 
% 
% %Amplitude and Phase Spectrum
% figure(78)
% stem(f15,abs(X_norm15),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel('Amplitude')
% title('900kW NM Auto FFT Amp')
% 
% % figure(79)
% % stem(f15,angle(X_norm15),'*','linewidth',2), xlabel('frequency (Hz)'), ylabel ('Phase Angle')
% % title('900kW NM Auto FFT Phase')
% 
% % %table of values (Coeffs, Freq, Amplitude, Angle)
% % table16 = table(X_norm15,f15',abs(X_norm15),angle(X_norm15));
% % table16.Properties.VariableNames = {'FFTCoeffs','Frequency','Amplitude','Phase'};
% % disp(table16)
% 
% %Reconstruct Using k Dominant frequencies
% X_recon15=X15;
% k=30; %keeps first k and last k + the mean
% for i=(k+10):(L15-k) %starting with the 40ith term and everything in between =0
%     X_recon15(i)=0;
% end 
% 
% %IFFT
% x_recon15= ifft(X_recon15);
% figure(80)
% plot(t15,x_t15,'k-*','LineWidth',2) %raw data
% hold on, grid on
% plot(t15,x_recon15,'b--','LineWidth',2) %reconstructed
% xlabel('Time(s)'), ylabel('Current (Amps)')
% title('900kW NM Auto Current Height vs. Time')
12 - append(num2str(k),'FFT Dominant Terms');