function saveCTResults(results,z_o,DirectoryPath)
    set(0,'DefaultFigureVisible', 'off');
    fig=figure;
    scatter(results(:,1),results(:,2));
    xlim([0 50]);
    xlabel('r');
    ylabel('C(r,z)');
    title(['CT output: z = ',num2str(z_o)]);
    saveas(fig,[DirectoryPath,'/CT_output.png']);
end

