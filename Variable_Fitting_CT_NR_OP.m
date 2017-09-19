function results=Variable_Fitting_CT_NR_OP(data,z_o,compound,power,DirectoryPath,mu,Fitting)
    h= figure('Visible','off');
    z = z_o;
    X = data(:,1);
    I = data(:,2);
     I2 = data(:,2);
    X2 = X.^2; %squaring x terms before poly fitting gets rid of odd power coefs
    
    
    if isequal(Fitting,'Polyfit')
        [CT,id,order] = polyFitOP(power,z_o,X2,X,I,z,DirectoryPath)
    elseif isequal(Fitting,'CubicSpline')
        r_max = X(length(X));
        [CT,id,order] = cubicSplineOP(z_o,power,r_max,X,X2,I,DirectoryPath)
    elseif isequal(Fitting,'gridfit')
        r_max = X(length(X));
        [CT,id,order] = gridFitFittingOP(z_o,power,r_max,X,X2,I,I2,DirectoryPath)
    end
    
    

    results=CT_functionOP(CT,id,compound,order,z,mu);  %with poly fitting

end
