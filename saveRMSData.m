function saveRMSData(data,SavesFolder,CV)
%Saves Figure to SavesFolder

        set(0,'DefaultFigureVisible', 'off');
        fig = figure; 

        scatter(data(:,1),data(:,2));
 
        title(['RMS Values for Different Smoothness at ', CV]);
        xlabel('Smoothness');
        ylabel('RMS(%)');
        
        saveas(fig,[SavesFolder,'/RMSData',CV,'.png'],'png')
end