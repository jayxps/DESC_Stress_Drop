% read in the sac data
load rh0eacc.mat;
data_mat = ans;
[tt,data,hdr] = fget_sac('rh0eacc.sac');

% you need to plot the data with tt

figure(1);
plot(tt,data,'r-',tt,data_mat,'b--');
xlabel('Time (s)');
ylabel('Amplitude');
legend('SAC data','your .mat data');

% compare the spectra
[MX,Phase,ff]=sacfft(data,tt);

% or load it from SAC's fft program
[ff_sac,data_sac,hdr] = fget_sac('rh0eacc.sac.fft');

figure(2);
plot(ff,MX,'b.-');
hold on;
plot(ff_sac,data_sac,'r--');
hold off;
xlabel('Frequency (Hz)');
ylabel('Amplitude spectra');
legend('sacfft.m','FFT in SAC');
