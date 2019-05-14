function [C, err, F]=cohExtraction(signal1, signal2, sampf, freqSpace)
% Extracts coherence between two signals
%
% [C, err, F]=cohExtraction(signal1, signal2, sampf)

iter=length(freqSpace);

signal1=signal1-mean(signal1);
signal2=signal2-mean(signal2);

F=[];
C=[];
err=[];

for a=1:iter-1
    
    %% Data Crunching  
    % Bandpass filtering to get data into frequency bins
    freq=freqSpace(a);
    freqStep=(freqSpace(a+1)-freqSpace(a))/2;
    F=[F; freq];
    
    [bb,aa]= butter(3,[2*(freq-freqStep)/sampf 2*(freq+freqStep)/sampf],'bandpass');

    filtSignal1=filter(bb,aa,signal1);
    filtSignal2=filter(bb,aa,signal2);

    fitLength=floor(2/(freq/sampf));    
        
    temp=zeros(floor(length(filtSignal1)/fitLength)-2,1);
    
    for j=1:floor(length(filtSignal1)/fitLength)-2

        cut1=filtSignal1(j*fitLength:(j+1)*fitLength);
        cut2=filtSignal2(j*fitLength:(j+1)*fitLength);
  
        crossCor=sqrt(max(abs(xcorr(cut1,cut2))).^2./(max(abs(xcorr(cut1))).*max(abs(xcorr(cut2)))));
        
        temp(j)=crossCor;

    end
        
    C=[C mean(temp)];
    err=[err std(temp)/sqrt(length(temp))];
    
end
