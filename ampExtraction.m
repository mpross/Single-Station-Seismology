function [A,F]=ampExtraction(signal,sampf)
% Extracts complex amplitudes of signal
%
% [F,A]=ampExtraction(signal,sampF)

startFreq=0.025;
freqStep=.01;
endFreq=0.5;

iter=floor((endFreq-startFreq)/freqStep);

signal=signal-mean(signal);
F=(startFreq:freqStep:endFreq);
A=[];
err=[];
for a=0:iter
    
    temp=[];
    r2=[];
    %% Data Crunching  
    % Bandpass filtering to get data into frequency bins
    freq=(startFreq+a*freqStep);
    [bb,aa]= butter(3,[2*((a-1/2)*freqStep+startFreq)/sampf 2*((a+1/2)*freqStep+startFreq)/sampf],'bandpass');

    filtSignal=filter(bb,aa,signal);

    fitLength=floor(1/(freq/sampf)/4);
    plot(filtSignal)
    for j=1:floor(length(filtSignal)/fitLength)-2

        tim=(j*fitLength:(j+1)*fitLength)'./sampf;
        cut=filtSignal(j*fitLength:(j+1)*fitLength);
        
        x=[sin(2*pi*freq*tim), cos(2*pi*freq*tim)];
        
        w=cut'*x*inv(x'*x);
        
        temp=[temp abs(w(2))+abs(w(1))*i];
        
    end
    A=[A mean(temp)];
    err=[err std(err)];
end


figure(8)
hold on
plot(F,abs(A));
fill([F, fliplr(F)], [err, fliplr(err)],'b','LineStyle','None');
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