function savePolyFitFigure(r,poly,X,I,order,z,DirectoryPath,id)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
            set(0,'DefaultFigureVisible', 'off');
            plot(r,polyval(poly,r.^2),X,I,'X');
            xlabel('X-Position mm');
            ylabel('Integrated Absrobance a.u. cm^-1');
            title(['z = ',num2str(z),', ',num2str(order) 'th order fit']);
            axis([X(1)-1 max(X) -1 max(I)+5]);            
            saveas(gcf,[DirectoryPath,'/Polyfit_',id],'jpg'); 
            saveas(gcf,[DirectoryPath,'/Polyfit_',id],'fig');

end

