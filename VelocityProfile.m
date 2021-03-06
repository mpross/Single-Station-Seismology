%% Set-up
close all
% clear all

% try
%     parpool();
% catch
% end

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.port', '587');
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.starttls.enable','true');

setpref('Internet','E_mail','mprossmatlab@gmail.com');
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username','mprossmatlab');
setpref('Internet','SMTP_Password',password);
% 
testBool=false;

earthquakes=["Mexico_5_9" "Oklahoma_4_4" "Indonesia_6_9" "Fiji_8_2" "CostaRica_6_1" ...
    "Fiji_6_8" "Oregon_6_2" "Venezuela_7_3" "Peru_7_1" "Fiji_7_8" "NewZealand_6_9" "Canada_6_6" "Iceland_6_8" ...
    "Peru_7_0" "Peru_7_5" "Papua_New_Guinea_7_5" "Peru_8_0" "El_Salvador_6_6"];
timeStamp=[1214366228 1212587999 1218725806 1218673195 1218583362 ...
    1218688157 1218965525 1218922324 1219136664 1220284172 1220588360 1224221998 1225763398 ...
    1235465459 1234865860 1241873924 1242891692 1243242230];

exclude=["Oklahoma_4_4" "Indonesia_6_9" "CostaRica_6_1" "Fiji_6_8" "Oregon_6_2" "Fiji_7_8"]; 
% exclude=[""];
  
%http://ds.iris.edu/spud/earthmodel/9991844
PREMdepth=[0 12858 25716 38574 51432]/1e3;
PREMdens=[13088.50 13088.46 13088.36 13088.18 13087.92]*0.0001;
PREMvp=[11262.20 11262.17 11262.10 11261.97 11261.79]/1e3;
PREMvs=[3667.80 3667.78 3667.73 3667.64 3667.51]/1e3;

% Love and Rayleigh phase-velocity maps, 5–40 s, of the western and central USA from USArray data
RefFreq=1./[5 10 20 40];
RefVel=[2.93*0.4800 3.19*0.61 3.51*0.96 3.88*0.98]*1e3;

% freqSpace=logspace(-2,0,20);
freqSpace=linspace(1,1e2)*1e-2;

clipPass=zeros(1, length(earthquakes));
vel=[];
vFreq=[];
vErr=[];

sampF=8;
t0=cputime;
% 
% fig8=figure(8);
% polarhistogram([],20,'Normalization','probability')   
% legend('','Mexico','Fiji','Venezula','Peru','NewZealand','Canada','Iceland','Peru','Peru')
% hold on

%% Code tests
if(false)
    disp(' ')

    % Amplitude fitting test
        tim=(1:1e4)/sampF;
        testA = 1e-9;
        ampTest=testA*sin(2*pi*0.1*tim);
        [ATest, ETest, F] = ampExtraction(ampTest', sampF);

        if(and(abs(ATest(find(and(F>0.099, F<0.101))))>8.5e-10,...
                abs(ATest(find(and(F>0.099, F<0.101))))<1.1e-9))
            disp('Fitting Test Passed')
        else
            disp('Fitting Test Failed')
            disp(['Amplitude found: ' num2str(abs(ATest(find(and(F>0.099, F<0.101)))))])
            disp(['Should be: ' num2str(testA)])
            return
        end
    % Sensor clipping check
        clippingCheck
        disp(['Clipping Check ' num2str(sum(clipPass)) ' out of ' num2str(length(clipPass)) ' passed'])
        disp(' ')
end
    
%% Data pull and decimate
for j=1:length(earthquakes)
% for j=length(earthquakes)-1

    if and(or(clipPass(j)==1, not(testBool)),sum(earthquakes(j)==exclude)==0)
        earthquakes(j)
        filename=strcat('/home/michael/Google Drive/Seismology/Data/GPS',num2str(timeStamp(j)),'_',earthquakes(j));

        rawData=load(filename);

        rawRY=rawData.rawData(1);
        inRY=decimate(rawRY.data,rawRY.rate/8)*1e-9;

        rawRX=rawData.rawData(2);
        inRX=decimate(rawRX.data,rawRX.rate/8)*1e-9;

        rawX=rawData.rawData(3);         
        inX=decimate(rawX.data,rawX.rate/8)*1e-9;

        rawY=rawData.rawData(4); 
        inY=decimate(rawY.data,rawY.rate/8)*1e-9;

        rawZ=rawData.rawData(5); 
        inZ=decimate(rawZ.data,rawZ.rate/8)*1e-9;

        inRY=inRY(1:1e5);
        inRX=inRX(1:1e5);
        inX=inX(1:1e5);
        inY=inY(1:1e5);
        inZ=inZ(1:1e5);


        %% Time cut
        timeThreshold=0.5e-6;

        timeCut=find(abs(inZ-mean(inZ))>timeThreshold);
        startTime=(timeCut(1)-1000*sampF);

        timeCut=find(abs(fliplr(inZ')-mean(inZ))>timeThreshold);
        endTime=(length(inZ)-(timeCut(1)-2000*sampF));
        if(startTime<0)
            startTime=1;
        end

        if(endTime>length(inZ))
            endTime=length(inZ);
        end

        %% Inversion and filtering

        STSInvertFilt = zpk(-2*pi*pairQ(8.33e-3,0.707),-2*pi*[0 0],1);
        STSInvertFilt = 1*STSInvertFilt/abs(freqresp(STSInvertFilt,2*pi*100));

        tim=(startTime:endTime)/sampF;

        [b,a]=butter(3,0.01*2/sampF,'high');

        RY=inRY(startTime:endTime)-mean(inRY(startTime:endTime));
        RX=inRX(startTime:endTime)-mean(inRX(startTime:endTime));

        X=lsim(STSInvertFilt,inX(startTime:endTime), tim-startTime/sampF);
        Y=lsim(STSInvertFilt,inY(startTime:endTime), tim-startTime/sampF);
        Z=lsim(STSInvertFilt,inZ(startTime:endTime), tim-startTime/sampF);

        X=filter(b,a,X);
        Y=filter(b,a,Y);
        Z=filter(b,a,Z);
        RY=filter(b,a,RY);
        RX=filter(b,a,RX);
        
        Z=Z(250*sampF:end);
        X=X(250*sampF:end);
        Y=Y(250*sampF:end);        
        RY=RY(250*sampF:end);
        RX=RX(250*sampF:end);
        
        tim=tim(250*sampF:end);

        %% Coherence
        [CX, ERCX, ~]=cohExtraction(RX, Z, sampF, freqSpace);
        [CY, ERCY, F]=cohExtraction(RY, Z, sampF, freqSpace);
        
        in=freqSpace(find(sqrt(CX.^2+CY.^2)/sqrt(2)>0.5));
        %% Spectra
        
        [AV, EV, F] = velExtraction(Z, X, Y, RX, RY, sampF, in, freqSpace);

        %% Phase Velocity Calculations

        vel=[vel; AV'];
        vFreq=[vFreq; F];
        vErr=[vErr; EV'];
        
    else
        disp(earthquakes(j))
        disp('Skipping this event.')
    end
end
%%
vAv=[];
fAv=[];
aAv=[];
errAv=[];
for i=1:length(freqSpace)
    if not(isnan(mean(vel(find(vFreq==freqSpace(i))))))
        vAv=[vAv; mean(vel(find(vFreq==freqSpace(i))))];
        errAv=[errAv; std(vel(find(vFreq==freqSpace(i))))];
        fAv=[fAv; freqSpace(i)];
    end
end

%% Plots
figure(1)
plot1=plot(tim,RX,tim,RY,tim,Z/4000);
grid on
set(plot1,'LineWidth',1.5);
set(gca,'FontSize',16);
legend('RX','RY')

figure(2)
plot1=plot(tim,X,tim,Y,tim,Z);
grid on
set(plot1,'LineWidth',1.5);
set(gca,'FontSize',16);

t=(cputime-t0)/3600

fig1=figure(4);
% plot2=errorbar(vFreq,vel/1e3,abs(vErr(:,1))/1e3,'.');
% plot2=errorbar(F(in),AV(in),EV(in),'.');
hold on
plot11=plot(vFreq,vel/1e3,'.');
plot2=errorbar(fAv,vAv/1e3,abs(errAv/1e3));
plot10=plot(RefFreq,RefVel/1e3,'.');
hold off
ylabel('Velocity (km/s)')
xlabel('Frequency (Hz)')
set(gca,'XScale','log');
ylim([0 5])
set(plot2,'LineWidth',1.5);
set(gca,'FontSize',16);
set(plot2,'MarkerSize',2);
set(plot10,'MarkerSize',30);
% set(plot11,'LineWidth',3);
legend('Single Station','Average','USArray Map')
grid on

figure(77)
plot2=semilogx(fAv,abs(errAv/1e3));

%% Fit

layers=3;
[bestPar,bestDispers]=dispersionFit(fAv,vAv,layers);
if(layers==3)
    bestDepth=bestPar(1)*heaviside(-(1:5e4)+bestPar(3))+bestPar(4)*heaviside(-(1:5e4)+bestPar(6)).*heaviside((1:5e4)-bestPar(3))+bestPar(7)*heaviside((1:5e4)-bestPar(6));
elseif(layers==2)
    bestDepth=bestPar(1)*heaviside(-(1:5e4)+bestPar(3))+bestPar(4)*heaviside((1:5e4)-bestPar(3));
end

hold on
plot2=plot(vAv,bestDispers);
ylabel('Velocity (m/s)')
xlabel('Frequency (Hz)')
set(gca,'XScale','log');
ylim([0 3000])
set(plot2,'LineWidth',1.5);
set(gca,'FontSize',16);
set(plot2,'MarkerSize',16);
hold off

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
legend('Measured','Density of silt loam soil','Density of quartz','Reference earth (PREM)')
xlim([1 3])
ylim([-40 0])

figure(7)
histogram(vErr(find(vErr<1e-6)),100,'Normalization','probability')
xlabel('Sum of Squared Error')
% xlim([0, 1e-6])

print(fig1,'-dpng','Rayleigh_Dispersion.png');
print(fig2,'-dpng','Velocity_Depth.png');
print(fig3,'-dpng','Density_Depth.png');
sendmail('mpross2@uw.edu','Earthquake Analysis Complete',"Completion Time: "+num2str(t)+" hours"+newline+...
   newline+"Earthquakes: "+strjoin(earthquakes)+newline+"Times: "+num2str(timeStamp),{'Rayleigh_Dispersion.png','Velocity_Depth.png', 'Density_Depth.png', 'dispersionSwarm.avi'});
