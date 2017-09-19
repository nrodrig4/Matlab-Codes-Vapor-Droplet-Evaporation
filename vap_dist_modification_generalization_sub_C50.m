function [ rho, rho_gradr, rho_gradz ] = vap_dist_modification_generalization_sub_C50(r_star, scaler, z_star,scalez, Cv,Cf)

r_star = scaler*r_star;
z_star = scalez*z_star;
% Written by PKZ, Feb. 9, 2011.
% vap_dist computes the vapor phase density distribution above an
% evaporating drop according to the solution of the steady-state Laplace
% equation. The solution to the s-s Laplace equation for "Weber's Disc" is
% given on p. 42-43 in "The Mathematics of Diffusion", 2nd Ed. by J. Crank.
% The boundary conditions used in that book are different than our boundary
% conditions. I transformed the variables to obtain the proper boundary
% conditions. My expressions can be found in my lab notebook on pages dated
% 2/8/11 and 2/9/11 (pp. 200-203).
%
% M = molar mass [g/mole]
% Pv = equilibrium vapor pressure [kPa]
% T = ambient temperature [K]
% Ru = universal gas constant [J/(mole K)]
% rho = density distribution of the vapor above the evaporating disc (drop)
%       in [kg/m^3].
%       The distribution is given in terms of dimensionless coordinates,
%       r_star = r/a and z_star = z/a where a is the radius of the disc
%       (drop). See my lab notebook for the mathematical expression.
%       Density at any location is given as rho(r_star, z_star).
% rho_grad = density gradient along the surface of the disc (drop) in
%       [kg/m^3]. This gradient is with respect to the dimensionless z
%       coordinate, z_star, i.e. d(rho)/dz_star. Since dz_star/dz = 1/a
%       then d(rho)/dz = (1/a)*d(rho)/dz_star. Be careful about units! The
%       disc (drop) radius is generally given in mm whereas the gradient is
%       given in kg/m^3.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Altered 6/8/17-6/10/17 by Nic Rodriguez


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initializing large matrices, rho and rho_grad
rho = zeros(size(r_star,2),size(z_star,2));
rho_gradz = zeros(size(r_star,2),size(z_star,2));
rho_gradr = zeros(size(r_star,2),size(z_star,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:size(r_star,2)
    for j = 1:size(z_star,2)
        alpha = (r_star(i)^2 + z_star(j)^2-1);
        rho(i,j) = (2/pi)*Cv*atan(1/sqrt(0.5*(alpha+sqrt(alpha^2+4*z_star(j)^2))))-Cf;
        if rho(i,j) < 0
            rho(i,j) = 0;
        end
        beta=sqrt(r_star(i)^4+2*r_star(i)^2*z_star(j)^2-2*r_star(i)^2+z_star(j)^4+2*z_star(j)^2+1);
        delta=2*r_star(i)^2+2*z_star(j)^2-2+2*beta;
        rho_gradr(i,j) = -(4*(beta + r_star(i)^2 + z_star(j)^2 - 1.0)*r_star(i)*Cv)/...
            (sqrt(delta)*pi*(r_star(i)^2 + z_star(j)^2 + 1 + beta)*beta);
        rho_gradz(i,j) = -4*z_star(j)*Cv/(beta*pi*sqrt(delta)); 
        
        if (abs(z_star(j))<1e-5)
            if (abs(r_star(i)) < 1+eps)
                rho_gradr(i,j) = 0.0;
                rho_gradz(i,j) = 0.0;
            end 
        end
        
        if(z_star(j) < 0.0)
            rho(i,j) = 0.0;
            rho_gradr(i,j) = 0.0;
            rho_gradz(i,j) = 0.0;
        end 
    end
end
    rho_gradr = scaler*rho_gradr;
    rho_gradz = scalez*rho_gradz;

end

