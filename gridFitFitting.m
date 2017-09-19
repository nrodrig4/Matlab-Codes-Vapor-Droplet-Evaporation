function [CT,id,order] = gridFitFitting(app,z_o,power,r_max,X,X2,I,I2,DirectoryPath)
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
     scatter(xgrid(1,:),zgrid(1,:))
    ax = app.UIAxes;
    plot(ax,xgrid(1,:),zgrid(1,:),X,I,'X')    
    title(app.UIAxes, 'GridFit'), axis([X(1)-1 max(X) -1 max(I)+5]);
    xlabel(app.UIAxes, 'X-Position mm');
    ylabel(app.UIAxes, 'Integrated Absrobance a.u. cm^-1');
    axis(app.UIAxes,[X(1)-1 max(X) -1 max(I)+5]);
    
    X = xgrid(1,:);
    I = zgrid(1,:)
    
    for order=power
        id = strcat('z_',num2str(z_o),'_order_',num2str(order));
        
%         if mod(order,2)==0  %This ensures that the polynomial order is even.
%     %     if mod(order,2)>0   % mod produces the remainder after division. If the remainder is greater than zero, the order must be odd.
%            figure          %Here the number is even and code executes
%             r=(X(1):0.1:X(end));% Specifies the range of X-values (mm) for the graph % If using 0's at the end change the X(end) to 30 %
%     %         poly=polyfit(X,I,order);
%     %                     plot(r,polyval(poly,r),X,I,'X'),xlabel('X-Position mm'),ylabel('Integrated Absrobance a.u. cm^-1'),title(['Z=' z ' ' num2str(order) 'th Order Fit']), axis([X(1)-1 max(X) -1 max(I)+5]);
%             poly=polyfit(X2,I,order);
%             plot(r,polyval(poly,r.^2),X,I,'X'),xlabel('X-Position mm'),ylabel('Integrated Absrobance a.u. cm^-1')
%             title(['z = ',num2str(z),' ',num2str(order),'th order fit']), axis([X(1)-1 max(X) -1 max(I)+5]);
%             saveas(gcf,[DirectoryPath,'/Polyfit_',id],'jpg'); 
%             saveas(gcf,[DirectoryPath,'/Polyfit_',id],'fig'); 
%         else
%             figure          %Here the number is odd and code executes
%             r=(X(1):0.1:X(end));% Specifies the range of X-values (mm) for the graph % If using 0's at the end change the X(end) to 30 %
%             %poly=polyfit(X,I,order);
%             poly=polyfit(X2,I,order);
%             plot(r,polyval(poly,r.^2),X,I,'X'),xlabel('X-Position mm'),ylabel('Integrated Absrobance a.u. cm^-1')
%             title(['z = ',num2str(z),', ',num2str(order) 'th order fit']), axis([X(1)-1 max(X) -1 max(I)+5]);            
%             saveas(gcf,[DirectoryPath,'/Polyfit_',id],'jpg'); 
%             saveas(gcf,[DirectoryPath,'/Polyfit_',id],'fig');
%         end
% 
%         CTraw=polyval(poly,(0:.5:50).^2); 
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
%     id = strcat('z_',num2str(z_o));
%     CT=zeros(length(I),power+1);
%     CT(:,1)=(0:.5:r_max)';
%     CT(:,length(CT(:,1))) = I';
%     r=(X(1):0.5:X(end));
% %     IAs=spline(X2,Zgrid(2,:),r.^2);
%      figure
%      plot(r,zgrid,X,I,'X')
%      xlabel('x-position (mm)'),ylabel('Integrated Absrobance (a.u. cm^{-1})');
%      title('GridFit Model Fit'), axis([X(1)-1 max(X) -1 max(I)+5]);
%      saveas(gcf,[DirectoryPath,'/GridFit_',id],'jpg'); 
%      saveas(gcf,[DirectoryPath,'/GridFit_',id],'fig'); 
%    
%     %m
    save([strcat(DirectoryPath,'/CT input_',' ',id,'.txt')],'CT','-ascii');
end

