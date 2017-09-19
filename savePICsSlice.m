function savePICsSlice(x,PIC,scaling,x2,PIC2,savepath)
        set(0,'DefaultFigureVisible', 'off');
        figure;
        plot(x,PIC,'x',scaling*x2,mean(PIC2,2),'o');     
        title('Projections');
        legend('Original PIC','Reconstructed PIC');
        xlabel('Beam Offset [mm]');
        ylabel('Integrated Ray [Absorbance]');
        saveas(gcf,[savepath,'/Projections'],'jpg'); %this saves the graphs
        saveas(gcf,[savepath,'/Projections'],'fig');
end

