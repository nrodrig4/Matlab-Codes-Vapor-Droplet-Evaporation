function Part2(app,M,C_v,rootdir,resultpos,sm_val,noise_flag)
Pv=15400; % Vapor pressure [Pa]
%C_v = (Pv*M)/(Ru*T); %0.2002
T=296.45; % room temperature [K]
 Ru=8.314472; % universal gas constant [J/(mole K)]
Rad = 6.5; %mm
D=16.71e-6; %m^lol
compound_theory = 4*Pv*M*D*Rad/(1000*Ru*T);
%methanol_theory =
%
%8.6971e-08
rw = 50;
z_oo = 0;
Cf = 2*C_v/pi*atan(Rad/(sqrt(0.5*(rw^2+z_oo^2-Rad^2+ ...
            sqrt((rw^2+z_oo^2-Rad^2)^2+4*z_oo^2*Rad^2)))));
 
  if length(sm_val) == 1     
        part2Norm(app,D,Rad,compound_theory,Cf,C_v,rootdir,resultpos,sm_val,noise_flag)
  else
        rootdir = [rootdir,'/OptimizedSmoothness'];
        if (exist(rootdir)~=1)
            mkdir(rootdir)
        end
        
        part2Optimimized(app,D,Rad,compound_theory,Cf,C_v,rootdir,resultpos,sm_val,noise_flag)
  end
         
end

