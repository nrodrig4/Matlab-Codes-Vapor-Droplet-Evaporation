function results=Variable_Poly_Fit_CT_NR(data,z_o,compound,power,DirectoryPath,mu)
    z = z_o;
    X = data(:,1);
    I = data(:,2);
    X2 = X.^2; %squaring x terms before poly fitting gets rid of odd power coefs
    
    m=zeros(6,15); %This creates a matrix of 0's where rach polynomial is entered. If order is changed, so should M
    CT=zeros(100,7); %This is another matrix of 0's so that CT input will also be recorder for each order of polynomial. 
    warning('off','MATLAB:polyfit:RepeatedPointsOrRescale') %Due to high polyfit value warnings are created
    % for order=power_sym    %power_sym is 1/2 of the inputted power since the fit will be applied to x^2.
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
    save([strcat(DirectoryPath,'/CT input_',' ',id,'.txt')],'CT','-ascii')%saves variable CT

    results=CT_function_NR(CT,DirectoryPath,id,compound,order,z,mu);  
    
    fig=figure
    scatter(results(:,1),results(:,2))
    xlim([0 50])
    xlabel('r')
    ylabel('C(r,z)')
    title(['CT output: z = ',num2str(z_o)])
    saveas(fig,[DirectoryPath,'/CT_output.png'])
end
