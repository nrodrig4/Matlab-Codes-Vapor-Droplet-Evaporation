function [noise_IA,noiseEnabled] = getNoise(noise_flag)
%Determines if noise should be added
    
    if noise_flag == 1.0
        noiseEnabled = 'WithNoise'
        noise_IA = 1.35; %a.u.*cm^{-1}
    else
        noiseEnabled = 'WithoutNoise'
        noise_IA = 0.0;
    end

end

