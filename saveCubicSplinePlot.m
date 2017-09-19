function saveCubicSplinePlot(r,IAs,X,I,id,DirectoryPath)
    plot(r,IAs,X,I,'X')
    xlabel('x-position (mm)'),ylabel('Integrated Absrobance (a.u. cm^{-1})');
    title('Cubic Spline'), axis([X(1)-1 max(X) -1 max(I)+5]);
    saveas(gcf,[DirectoryPath,'/Cubic_Spline_',id],'jpg'); 
    saveas(gcf,[DirectoryPath,'/Cubic_Spline_',id],'fig'); 
end

