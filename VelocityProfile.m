%% Set-up
% close all
% clear all

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.port', '587');
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.starttls.enable','true');

setpref('Internet','E_mail','mprossmatlab@gmail.com');
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username','mprossmatlab');
setpref('Internet','SMTP_Password',password);

earthquakes=["Mexico_5_9" "Oklahoma_4_4" "Indonesia_6_9" "Fiji_8_2" "CostaRica_6_1" ...
    "Fiji_6_8" "Oregon_6_2" "Venezuela_7_3" "Peru_7_1" "Fiji_7_8" "NewZealand_6_9.mat" "Canada_6_6.mat" "Iceland_6_8.mat"];
timeStamp=[1214366228 1212587999 1218725806 1218673195 1218583362 ...
    1218688157 1218965525 1218922324 1219136664 1220284172 1220588360 1224221998 1225763398];

vel=[];
vFreq=[];
sampF=8;
t0=cputime;

%% Data pull and decimate
% for j=1:length(earthquakes)
for j=1
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
    % 
    % STSX=inSTSX(startTime:endTime);
    % STSY=inSTSY(startTime:endTime);
    % STSZ=inSTSZ(startTime:endTime);

    STSX=filter(b,a,STSX);
    STSY=filter(b,a,STSY);
    STSZ=filter(b,a,STSZ);
    BRSX=filter(b,a,BRSX);
    BRSY=filter(b,a,BRSY);
    

    %% Spectra
    avg=15;
    [ABRSX, ~] = ampExtraction(BRSX, sampF);
    [ABRSY, ~] = ampExtraction(BRSY, sampF);
    [ASTSX, ~] = ampExtraction(STSX, sampF);
    [ASTSY, ~] = ampExtraction(STSY, sampF);
    [ASTSZ, F] = ampExtraction(STSZ, sampF);

    [COH,~]=coh2(BRSX,BRSY,1/sampF, avg, 1, @hann);
    [COHX,~]=coh2(BRSX,STSX,1/sampF, avg, 1, @hann);
    [COHY,F2]=coh2(BRSY,STSY,1/sampF, avg, 1, @hann);

    [T,~]=tfe2(sqrt(BRSX.^2+BRSY.^2),STSZ,1/sampF, avg, 1, @hann);
    [TX,~]=tfe2(BRSX,STSX,1/sampF, avg, 1, @hann);
    [TY,F3]=tfe2(BRSY,STSY,1/sampF, avg, 1, @hann);

    %% Phase Velocity Calculations

    thresh=0.95;
%     Cin=find(and(sqrt(COHX.^2+COHY.^2)>thresh,F2<0.5));
    Cin=find(and(or(abs(ABRSY)>1e-9,abs(ABRSX)>1e-9),abs(ASTSZ)>1e-6));
    % cohV=ASTSZ(Cin)./sqrt(ABRSY(Cin).^2+ABRSX(Cin).^2);
    cohV=abs(ASTSZ(Cin)./sqrt(ABRSY(Cin).^2+ABRSX(Cin).^2));
%     cohV=abs(ASTSZ./sqrt(ABRSY.^2+ABRSX.^2));
    % cohV=abs(real(T(Cin)));
    cohF=F(Cin);
%     cohF=F;
    vel=[vel; cohV'];
    vFreq=[vFreq; cohF'];
    
    obsDispers=movmean(vel,40);
    % depth=obsDispers./vFreq;
end
%% Fit


layers=3;
[bestPar,bestDispers]=dispersionFit(vFreq,vel,layers);
if(layers==3)
    bestDepth=bestPar(1)*heaviside(-(1:5e4)+bestPar(3))+bestPar(4)*heaviside(-(1:5e4)+bestPar(6)).*heaviside((1:5e4)-bestPar(3))+bestPar(7)*heaviside((1:5e4)-bestPar(6));
elseif(layers==2)
    bestDepth=bestPar(1)*heaviside(-(1:5e4)+bestPar(3))+bestPar(4)*heaviside((1:5e4)-bestPar(3));
end
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
plot2=loglog(F,abs(ABRSY),F,abs(ABRSX),F,abs(ASTSX),F,abs(ASTSY),F,abs(ASTSZ));
grid on
set(plot2,'LineWidth',1.5);
set(gca,'FontSize',16);
ylabel('ASD (m/s or rad /\surd{Hz})')
xlabel('Frequency (Hz)')
legend('\theta_x','\theta_y','v_x','v_y','v_z')

t=(cputime-t0)/3600

fig1=figure(4);
plot2=plot(vFreq,vel,'.',vFreq,bestDispers);
ylabel('Velocity (m/s)')
xlabel('Frequency (Hz)')
set(plot2,'LineWidth',1.5);
set(gca,'FontSize',16);
set(plot2,'MarkerSize',16);

fig2=figure(5);
plot2=plot(bestDepth,-(1:5e4)/1e3);
xlabel('Velocity (m/s)')
ylabel('Depth (km)')
set(plot2,'LineWidth',1.5);
set(gca,'FontSize',16);
set(plot2,'MarkerSize',16);
ylim([-40 0])

bestDepth=bestDepth/1e3;
dens=1.6612*bestDepth-0.4721*bestDepth.^2+0.0671*bestDepth.^3-0.0043*bestDepth.^4+0.000106*bestDepth.^5;

fig3=figure(6);
plot2=plot(dens,-(1:5e4)/1e3,1.3+0*(1:5e4),-(1:5e4)/1e3,2.65+0*(1:5e4),-(1:5e4)/1e3);
xlabel('Density (g/cm^3)')
ylabel('Depth (km)')
set(plot2,'LineWidth',1.5);
set(gca,'FontSize',16);
set(plot2,'MarkerSize',16);
legend('Measured','Density of silt loam soil','Density of quartz')
xlim([1 3])
ylim([-40 0])

print(fig1,'-dpng','Rayleigh_Dispersion.png');
print(fig2,'-dpng','Velocity_Depth.png');
print(fig3,'-dpng','Density_Depth.png');
% sendmail('mpross2@uw.edu','Earthquake Analysis Complete',"Completion Time: "+num2str(t)+" hours"+newline+...
%    newline+"Earthquakes: "+strjoin(earthquakes)+newline+"Times: "+num2str(timeStamp),{'Rayleigh_Dispersion.png','Velocity_Depth.png','Density_Depth.png'});

