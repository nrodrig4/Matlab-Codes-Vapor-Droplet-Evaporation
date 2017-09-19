function saveMirrorPlot(Ir_x,conc,savepath,mu)
%Saves Plot of mirrored Line
        set(0,'DefaultFigureVisible', 'off');
        figure;
        plot(Ir_x,conc);
        
        t = sprintf('Slice of Concentration, through origin, mu=%.3e',mu);
        title(t);
	    ylabel('Concentration');
        xlabel('Radial Position [mm]');
        grid on;
        saveas(gcf,[savepath,'/Slice'],'jpg'); %this saves the graphs
        saveas(gcf,[savepath,'/Slice'],'fig');
end

