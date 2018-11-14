
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
t0=cputime;
for j=1:length(earthquakes)
    earthquakes(j)
    filename=strcat('/home/michael/Google Drive/Seismology/Data/GPS',num2str(timeStamp(j)),'_',earthquakes(j));
    if(not(exist(filename+".mat",'file')))
        Dt=3600*10;
        try
            rawData=get_data2([{'L1:ISI-GND_BRS_ITMX_RY_OUT_DQ','L1:ISI-GND_BRS_ITMY_RX_OUT_DQ','L1:ISI-GND_STS_ITMY_X_DQ'...
                'L1:ISI-GND_STS_ITMY_Y_DQ','L1:ISI-GND_STS_ITMY_Z_DQ'}],'raw',timeStamp(j),Dt);
            save(filename,'rawData')
        catch ex
            sendmail('mpross2@uw.edu','Error in Earthquake Data Pull',ex.message);
            ex.message
            break
        end
    end 
end

t=(cputime-t0)/3600

dataSize=0;
s=dir('/home/michael/Google Drive/Seismology/Data')
for k=1:length(s)
    dataSize=dataSize+s(k).bytes;
end
dataSize=dataSize/1e9

sendmail('mpross2@uw.edu','Earthquake Data Pull Complete',"Completion Time: "+num2str(t)+" hours"+newline+...
    "Data Size: "+num2str(dataSize)+" GB"+newline+"Earthquakes: "+strjoin(earthquakes)+newline+"Times: "+num2str(timeStamp));
