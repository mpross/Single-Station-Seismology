function [A,F]=ampExtraction(signal,sampF,cuts)
% Extracts complex amplitudes of signal
%
% [F,A]=ampExtraction(signal,sampF)

close all

signal=signal-mean(signal);

N=floor(length(signal)/cuts);
F=(2:N/2)/N*sampF;
A=[];

for i=0:cuts-1
    tempSignal=signal(i*N+1:(i+1)*N);
    tempA=fft(hanning(N).*tempSignal)/sqrt(2*N/sampF);
    A=[A tempA(2:floor(N/2))];
    
end

err1=abs(max(A'));
err2=abs(min(A'));
mag=mean(A');
A=mag;

figure(8)
hold on
plot(F,abs(mag));
fill([F, fliplr(F)], [err1, fliplr(err2)],'b','LineStyle','None');
alpha(0.1)
hold off
set(gca,'YScale','log')
set(gca,'XScale','log')

figure(9)
X=[real(A);imag(A)]';
hist3(X,[200 200],'CdataMode','auto','LineStyle','none')
colorbar
view(2)

figure(10)
histogram(real(A),100)
set(gca,'YScale','log')