function [dispers]=dispersionCalc(vP1,vS1,d1,vP2,vS2,d2,vP3,vS3,freq,layers)
% 
% Solves for Ralyeigh wave dispersion curve given a set of P and S wave
% velocities and depths for two or three layer models.
%
% [dispers]=dispersionCalc(vP1,vS1,d1,vP2,vS2,d2,vP3,vS3,freq,layers)

    % Solve for Rayleigh wave phase velocity for first layer
    coeffs=[1/vS1^6 0 -8/vS1^4 0 8/vS1^2*(3-2*vS1^2/vP1^2) 0 -16*(1-vS1^2/vP1^2)];
    rots=roots(coeffs);
    fitV1=max(rots(find(and(and(and(imag(rots)==0,rots>0),rots<vP1),vS3/vP3<1/sqrt(2)))));

    % Solve for Rayleigh wave phase velocity for second layer
    coeffs=[1/vS2^6 0 -8/vS2^4 0 8/vS2^2*(3-2*vS2^2/vP2^2) 0 -16*(1-vS2^2/vP2^2)];
    rots=roots(coeffs);
    fitV2=max(rots(find(and(and(and(imag(rots)==0,rots>0),rots<vP2),vS3/vP3<1/sqrt(2)))));

    % Solve for Rayleigh wave phase velocity for third layer
    if(layers==3)
        coeffs=[1/vS3^6 0 -8/vS3^4 0 8/vS3^2*(3-2*vS3^2/vP3^2) 0 -16*(1-vS3^2/vP3^2)];
        rots=roots(coeffs);
        fitV3=max(rots(find(and(and(and(imag(rots)==0,rots>0),rots<vP3),vS3/vP3<1/sqrt(2)))));
    end

    % Calculate surface dispersion curve for given parameters
    dispers=[];
    for j=(1:length(freq))
        if(layers==2)
            % Two layers
            coeffs=[1 -fitV2 (fitV2-fitV1)*d1*freq(j)];
        elseif(layers==3)
            % Three layers
            coeffs=[1 -fitV3 (((fitV2-fitV1)*d1+(fitV3-fitV2)*d2)*freq(j))];
        else
            'Invalid layer number.'
            break
        end
        rots=roots(coeffs);        
        if isempty(max(rots(find(and(rots>0,imag(rots)==0)))))
            coeffs=[1 -fitV2 (fitV2-fitV1)*d1*freq(j)];
            rots=roots(coeffs);
            if isempty(max(rots(find(and(rots>0,imag(rots)==0)))))
                dispers=[dispers; nan];
            else
                dispers=[dispers; max(rots(find(and(rots>0,imag(rots)==0))))];
            end
        else
            dispers=[dispers; max(rots(find(and(rots>0,imag(rots)==0))))];
        end
    end
end