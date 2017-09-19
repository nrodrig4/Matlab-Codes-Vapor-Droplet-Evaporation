%IA function 
function  mat = Integrated_Absorbance_minus_Concentration_mod(noise_IA,z_o,mu,M,C_v,R,SavesFolder,x_spacing)
   
    %Setting up General Parameters and pre-allocating memory for speed
    r = 50;
       
    x_loc = x_spacing';  
    mat = zeros(length(x_loc),2);
        
    %creating file name for Integrated Absorbance Data
    filename = strcat(SavesFolder,'/IA_',num2str(z_o),'.txt');
    
    rw = 50;  
    z_oo = 0;
    %Weber's Disc Analytical Solution
    C_f = 2*C_v/pi*atan(R/(sqrt(0.5*(rw^2+z_oo^2-R^2+ ...
    sqrt((rw^2+z_oo^2-R^2)^2+4*z_oo^2*R^2))))) %analytical C at (rw,z_oo)

    %calculating Integrated Absorbance (minus Weber's disc at C(50,0)
    for i = 1:length(x_loc)

        x_o = x_loc(i);
        a = sqrt(r^2-x_o^2);

        fun = @(t) compute_C(t,C_v,x_o,z_o,R,C_f);
        %disp('int1')
        solution = integral(fun,-a,a); %line integral along the y-direction
        %disp('int2')
        solution = solution*1000/(1*10^6*M)*mu-(2.0*rand(size(solution))-1*ones(size(solution)))*noise_IA; %unit?
        
        mat(i,1) = x_o;
        mat(i,2) = solution;

    end    
   
dlmwrite(filename,mat,'delimiter','\t','precision',10)
     
end
    
function fun = compute_C(t,C_v,x_o,z_o,R,C_f)
   fun = 2.*C_v./pi.*atan(R./(sqrt(0.5.*(x_o.^2+t.^2+z_o.^2-R.^2+ ...
        sqrt((x_o.^2+t.^2+z_o.^2-R.^2).^2+4.*z_o.^2.*R.^2)))))-C_f; %t represents y coordinates
   
   if (fun < 0.0)
       fun = 0.0*fun;
   end
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
