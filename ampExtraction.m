function [A, err, F, list]=ampExtraction(signal,sampf)
% Extracts complex amplitudes of signal
%
% [F,A]=ampExtraction(signal,sampF)

startFreq=0.01;
freqStep=.005;
endFreq=1;

iter=floor((endFreq-startFreq)/freqStep);

signal=signal-mean(signal);
F=(startFreq:freqStep:endFreq);

A=[];
err=[];
list=[];

for a=0:iter
    
    temp=[];
    r2=[];
    %% Data Crunching  
    % Bandpass filtering to get data into frequency bins
    freq=(startFreq+a*freqStep);
    [bb,aa]= butter(3,[2*((a-1/2)*freqStep+startFreq)/sampf 2*((a+1/2)*freqStep+startFreq)/sampf],'bandpass');

    filtSignal=filter(bb,aa,signal);

    fitLength=floor(1/(freq/sampf)/4);
    for j=1:floor(length(filtSignal)/fitLength)-2

        tim=(j*fitLength:(j+1)*fitLength)'./sampf;
        cut=filtSignal(j*fitLength:(j+1)*fitLength);
        
        x=[sin(2*pi*freq*tim), cos(2*pi*freq*tim)];
        
        w=cut'*x*inv(x'*x);
        if (sqrt(w(1)^2+w(2)^2) >= max(abs(signal))/100)
            temp=[temp abs(w(2))+abs(w(1))*i];
        end
        
    end
%     h=histogram(abs(temp), 20);
    
    A=[A mean(temp)];
    err=[err std(temp)/sqrt(length(temp))];
%     list=[list; h.Values];
end
