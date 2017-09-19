function [CT,id,order] = cubicSplineOP(z_o,power,r_max,X,X2,I,DirectoryPath)

    order = power; 
    id = strcat('z_',num2str(z_o));
    CT=zeros(101,power+1);
    CT(:,1)=(0:.5:r_max)';
     
    r=(X(1):.1:X(end));
    IAs=spline(X2,I,r.^2);
    figure
    

   
    r=X.^2;
    IAs=spline(X2,I,(0:.5:50).^2);    
    CT(:,end)=IAs';
    save([strcat(DirectoryPath,'/CT input_',' ',id,'.txt')],'CT','-ascii');
end

