function results=Variable_Fitting_CT_NR(app,data,z_o,compound,power,DirectoryPath,mu,Fitting)
    h= figure('Visible','off');
    z = z_o;
    X = data(:,1);
    I = data(:,2);
     I2 = data(:,2);
    X2 = X.^2; %squaring x terms before poly fitting gets rid of odd power coefs
    
    
    if isequal(Fitting,'Polyfit')
        [CT,id,order] = polyFit(app,power,z_o,X2,X,I,z,DirectoryPath)
    elseif isequal(Fitting,'CubicSpline')
        r_max = X(length(X));
        [CT,id,order] = cubicSpline(app,z_o,power,r_max,X,X2,I,DirectoryPath)
    elseif isequal(Fitting,'gridfit')
        r_max = X(length(X));
        [CT,id,order] = gridFitFitting(app,z_o,power,r_max,X,X2,I,I2,DirectoryPath)
    end
    
    

    results=CT_function_NR(app,CT,DirectoryPath,id,compound,order,z,mu);  %with poly fitting
        
    ax = app.UIAxes_2;
    fig=figure;           
    scatter(ax,results(:,1),results(:,2));
    xlim(ax,[0 50]);
    xlabel(ax,'r');
    ylabel(ax,'C(r,z)');
    title(ax,['CT output: z = ',num2str(z_o)]);
    legend(ax,'CT Results');
    
   saveCTResults(results,z_o,DirectoryPath)
end
