function [CT,id,order] = gridFitFittingOP(z_o,power,r_max,X,X2,I,I2,DirectoryPath)
    order = power;
    
    for i = 1:length(X)/2
        Y(i) = 0;
    end
    for i2 = length(X)/2+1:length(X)
        if mod(length(X),2) == 0 
         Y(i2) = 1;
        else 
         Y(i2+0.5) = 1;
        end
    end
    xxL = size(X)
    yyL = size(Y)
    IL = size(I)
    
    [zgrid,xgrid,ygrid] = gridfit(X,Y',I,0:0.5:50,0:1,'smooth',5);
    zgrid = [zgrid(1,:);zgrid(2,:)];
%     zgrid
%     zgridS = size(zgrid);
%     xgridA = size(xgrid);
% %     pause
  
    
    X = xgrid(1,:);
    I = zgrid(1,:)
    
    for order=power
        id = strcat('z_',num2str(z_o),'_order_',num2str(order));
        

        CTraw = I;
        
        k=1;

        %while CTval>0
        for k=[1:length(CTraw)] %This ensure that the CT input won't have negative values. 
            if CTraw(k)<=0
               CTraw(k)=0;
               CTraw(k:length(CTraw))=0;
            else
               CTraw(k)=CTraw(k);
            end
            k=k+1;
        end

%         for eq=[1:length(poly)] %Enters each polynomial into m
%             m((order),eq)=poly(eq);
%         end

        for i=[1:length(CTraw)]  %input for the CT code, this enters a smooth curve produced by each polynomail fitting. 
            CT(i,(order+1))=CTraw(i);
            CT(i,1)=.5*i-.5;
        end
    end

%    
%     %m
    save([strcat(DirectoryPath,'/CT input_',' ',id,'.txt')],'CT','-ascii');
end

