function saveIAScatterFigure(data,SavesFolder,z_o)
%Saves Figure to SavesFolder

        set(0,'DefaultFigureVisible', 'off');
        fig = figure; 

        scatter(data(:,1),data(:,2));
 
        title(['IA output: z = ',num2str(z_o)])
        xlim([0 50]);
        xlabel('x axis');
        ylabel('IA');
        
        saveas(fig,[SavesFolder,'/IA_output.png'],'png')
end

