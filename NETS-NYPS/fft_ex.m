function [Magnitud, Fase]=fft_ex(Fs,signal,Nm)
t=0:1/Fs:(Nm*(1/Fs)-(1/Fs));
size=length(t);

sig=fft(signal);
Magnitud=abs(sig);
angang=angle(sig);
MMag=max(Magnitud);
[TT]=find(Magnitud==MMag);
% f_st = linspace(0,Fs,length(sig)); 
% figure(10)
% subplot(2,1,1)
% plot(f_st,Magnitud/(size/2))
% subplot(2,1,2)
% plot(f_st,angang)
% title('Waveform spectrum')
% xlabel('Frequency (Hz)')
% ylabel('Amplitude')
% ylim([0 0.002])
% xlim([0 1])
% xlim([0 Fs/32])
angulo=angle(sig);
Magnitud=MMag/(size/2);
Fase=angulo(1,TT);
end
