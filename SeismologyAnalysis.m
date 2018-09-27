%% Set-up
close all
clear all

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.port', '587');
% props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.starttls.enable','true');

setpref('Internet','E_mail','mprossmatlab@gmail.com');
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username','mprossmatlab');
setpref('Internet','SMTP_Password','matlab123');

earthquakes=["Mexico_5_9" "Oklahoma_4_4" "Indonesia_6_9" "Fiji_8_2" "CostaRica_6_1" ...
    "Fiji_6_8" "Oregon_6_2" "Venezuela_7_3" "Peru_7_1" "Fiji_7_8" "NewZealand_6_9.mat"];
timeStamp=[1214366228 1212587999 1218725806 1218673195 1218583362 ...
    1218688157 1218965525 1218922324 1219136664 1220284172 1220588360];

vel=[];
vFreq=[];
sampF=8;
t0=cputime;

for j=1
    %% Data pull and decimate
   
    earthquakes(j)
    filename=strcat('/home/michael/Google Drive/Seismology/Data/GPS',num2str(timeStamp(j)),'_',earthquakes(j));

    rawData=load(filename);

    rawBRSX=rawData.rawData(1);
    inBRSX=decimate(rawBRSX.data,rawBRSX.rate/8)*1e-9;

    rawBRSY=rawData.rawData(2);
    inBRSY=decimate(rawBRSY.data,rawBRSY.rate/8)*1e-9;

    rawSTSX=rawData.rawData(3); 
    inSTSX=decimate(rawSTSX.data,rawSTSX.rate/8)*1e-9;

    rawSTSY=rawData.rawData(4); 
    inSTSY=decimate(rawSTSY.data,rawSTSY.rate/8)*1e-9;

    rawSTSZ=rawData.rawData(5); 
    inSTSZ=decimate(rawSTSZ.data,rawSTSZ.rate/8)*1e-9;
    
    inBRSX=inBRSX(1:1e5);
    inBRSY=inBRSY(1:1e5);
    inSTSX=inSTSX(1:1e5);
    inSTSY=inSTSY(1:1e5);
    inSTSZ=inSTSZ(1:1e5);

    %% Time cut
    timeThreshold=3e-6;

    timeCut=find(abs(inSTSZ-mean(inSTSZ))>timeThreshold);
    startTime=(timeCut(1)-1000);

    timeCut=find(abs(fliplr(inSTSZ')-mean(inSTSZ))>timeThreshold);
    endTime=(length(inSTSZ)-(timeCut(1)-2000));
    if(startTime<0)
        startTime=1;
    end
    
    if(endTime>length(inSTSZ))
        endTime=length(inSTSZ);
    end

    %% Inversion and filtering

    STSInvertFilt = zpk(-2*pi*[pairQ(8.33e-3,0.707)],-2*pi*[0 0],1);
    STSInvertFilt = 1*STSInvertFilt/abs(freqresp(STSInvertFilt,2*pi*100));

    time=(startTime:endTime)/sampF;

    [b,a]=butter(3,0.01*2/sampF,'high');

    BRSX=inBRSX(startTime:endTime);
    BRSY=inBRSY(startTime:endTime);

    STSX=lsim(STSInvertFilt,inSTSX(startTime:endTime), time-startTime/sampF);
    STSY=lsim(STSInvertFilt,inSTSY(startTime:endTime), time-startTime/sampF);
    STSZ=lsim(STSInvertFilt,inSTSZ(startTime:endTime), time-startTime/sampF);

    STSX=filter(b,a,STSX);
    STSY=filter(b,a,STSY);
    STSZ=filter(b,a,STSZ);
    BRSX=filter(b,a,BRSX);
    BRSY=filter(b,a,BRSY);

    %% Spectra
    avg=9;
    [ABRSX, ~] = asd2(BRSX,1/sampF, avg, 1, @hann);
    [ABRSY, ~] = asd2(BRSY,1/sampF, avg, 1, @hann);
    [ASTSX, ~] = asd2(STSX,1/sampF, avg, 1, @hann);
    [ASTSY, ~] = asd2(STSY,1/sampF, avg, 1, @hann);
    [ASTSZ, F] = asd2(STSZ,1/sampF, avg, 1, @hann);
    
    [COH,~]=coh2(BRSX,BRSY,1/sampF, avg, 1, @hann);
    [COHX,~]=coh2(BRSX,STSZ,1/sampF, avg, 1, @hann);
    [COHY,F2]=coh2(BRSY,STSZ,1/sampF, avg, 1, @hann);

    [T,~]=tfe2(BRSX,BRSY,1/sampF, avg, 1, @hann);
    [TX,~]=tfe2(BRSX,STSZ,1/sampF, avg, 1, @hann);
    [TY,F3]=tfe2(BRSY,STSZ,1/sampF, avg, 1, @hann);

    %% Phase Velocity Calculations

    thresh=0.6;
    Cin=find(and(movmean(sqrt(COHX.^2+COHY.^2),10)>thresh,F2<0.5));
    cohV=ASTSZ(Cin)./sqrt(ABRSY(Cin).^2+ABRSX(Cin).^2);
    cohF=F(Cin);
    vel=[vel; cohV];
    vFreq=[vFreq; cohF];

    %% Plots
    figure(1)
    plot1=plot(time,BRSY,time,BRSX);
    grid on
    set(plot1,'LineWidth',1.5);
    set(gca,'FontSize',16);
    legend('RX','RY')

    figure(2)
    plot1=plot(time,STSX,time,STSY,time,STSZ);
    grid on
    set(plot1,'LineWidth',1.5);
    set(gca,'FontSize',16);

    figure(3)
    plot2=loglog(F,ABRSY,F,ABRSX,F,ASTSX,F,ASTSY,F,ASTSZ);
    grid on
    set(plot2,'LineWidth',1.5);
    set(gca,'FontSize',16);
    ylabel('ASD (m/s or rad /\surd{Hz})')
    xlabel('Frequency (Hz)')
    legend('\theta_x','\theta_y','v_x','v_y','v_z')

end

t=(cputime-t0)/3600
%% Stastics and Plot

[Fstat,I]=sort(vFreq);
Vstat=vel(I);
meanVel=[];
stdVel=[];
fVel=[];
N=10;
for j=0:length(Vstat)/N-1
    meanVel=[meanVel mean(Vstat(j*N+1:j*N+N))];
    stdVel=[stdVel std(Vstat(j*N+1:j*N+N))];
    fVel=[fVel mean(Fstat(j*N+1:j*N+N))];
end

fillF=[fVel, fliplr(fVel)];
fillV=[meanVel+stdVel, fliplr(meanVel-stdVel)];

fig1=figure(4)
plot2=plot(vFreq,vel,'.',fVel,meanVel);
hold on
fill(fillF, fillV,'r','LineStyle','none');
alpha(0.2)
hold off
ylabel('Velocity (m/s)')
xlabel('Frequency (Hz)')
set(plot2,'LineWidth',1.5);
set(gca,'FontSize',16);
set(plot2,'MarkerSize',16);

fig2=figure(5)
plot2=semilogx(vFreq,vel,'.',fVel,meanVel);
hold on
fill(fillF, fillV,'r','LineStyle','none');
alpha(0.2)
hold off
ylabel('Velocity (m/s)')
xlabel('Frequency (Hz)')
set(plot2,'LineWidth',1.5);
set(gca,'FontSize',16);
set(plot2,'MarkerSize',16);

print(fig1,'-dpng','Rayleigh_Dispersion_Linear.png')
print(fig2,'-dpng','Rayleigh_Dispersion_Log.png')
% sendmail('mpross2@uw.edu','Earthquake Analysis Complete',"Completion Time: "+num2str(t)+" hours"+newline+...
%     newline+"Earthquakes: "+strjoin(earthquakes)+newline+"Times: "+num2str(timeStamp),{'Rayleigh_Dispersion_Linear.png','Rayleigh_Dispersion_Log.png'});

