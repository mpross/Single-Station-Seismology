function [bestPar,bestDispers]=dispersionFit(obsFreq,obsDispers,n)
%% Fit
bestPar=[];
bestDispers=[];
bestErr=1e1000;
for j=(0:1e3)
    %Assumes two layers and depth is equal to wavelength   
    if(n==2)
        %First layer
        vP1=rand*3e3;
        vS1=rand*3e3; %must be between 0 and /srt(2) Landau
        d1=rand*10e3;

        %Second layer
        vP2=rand*3e3;
        vS2=rand*3e3;
        d2=0;

        %Third layer
        vP3=0;
        vS3=0;
    elseif(n==3)
        %First layer
        vP1=rand*3e3;
        vS1=rand*3e3; %must be between 0 and /srt(2) Landau
        d1=rand*10e3;

        %Second layer
        vP2=rand*3e3;
        vS2=rand*3e3;
        d2=rand*10e3;

        %Third layer
        vP3=rand*3e3;
        vS3=rand*3e3;
    end

    coeffs=[1/vS1^6 0 -8/vS1^4 0 8/vS1^2*(3-2*vS1^2/vP1^2) 0 -16*(1-vS1^2/vP1^2)];
    rots=roots(coeffs);
    fitV1=rots(find(and(and(imag(rots)==0,rots>0),rots<vP1)));

    coeffs=[1/vS2^6 0 -8/vS2^4 0 8/vS2^2*(3-2*vS2^2/vP2^2) 0 -16*(1-vS2^2/vP2^2)];
    rots=roots(coeffs);
    fitV2=rots(find(and(and(imag(rots)==0,rots>0),rots<vP2)));
    if(n==3)
        coeffs=[1/vS3^6 0 -8/vS3^4 0 8/vS3^2*(3-2*vS3^2/vP3^2) 0 -16*(1-vS3^2/vP3^2)];
        rots=roots(coeffs);
        fitV3=rots(find(and(and(imag(rots)==0,rots>0),rots<vP3)));
    end

    fitDispers=[];
    for j=(1:length(obsFreq))
        if(n==2)
            % Two layers square attenuation
            coeffs=[1 -fitV2 (fitV2-fitV1)*d1*obsFreq(j)];
            % Two layers linear attenuation
%             coeffs=[1 -1/2*d1^2*obsFreq(j)^2*(fitV1-fitV2) -obsFreq(j)^3*d1(fitV1-fitV2) obsFreq(j)^4*fitV2/2];
        elseif(n==3)
            % Three layers
            coeffs=[1 -fitV3 (((fitV2-fitV1)*d1+(fitV3-fitV2)*d2)*obsFreq(j))];
        else
            'Invalid layer number.'
            break
        end
        rots=roots(coeffs);        
        if isempty(max(rots(find(and(rots>0,imag(rots)==0)))))
            coeffs=[1 -fitV2 (fitV2-fitV1)*d1*obsFreq(j)];
            rots=roots(coeffs);
            if isempty(max(rots(find(and(rots>0,imag(rots)==0)))))
                fitDispers=[fitDispers; nan];
            else
                fitDispers=[fitDispers; max(rots(find(and(rots>0,imag(rots)==0))))];
            end
        else
            fitDispers=[fitDispers; max(rots(find(and(rots>0,imag(rots)==0))))];
        end
    end
    err=sum((fitDispers-obsDispers).^2);
    if err<bestErr
        bestDispers=fitDispers;
        bestErr=err;
        bestPar=[vP1 vS1 d1 vP2 vS2 d2 vP3 vS3];
    end
end
if(isempty(bestDispers))
        'Can not find solution'
        bestPar=[0 0 0 0 0 0 0 0];
        bestDispers=nan*obsFreq;
end