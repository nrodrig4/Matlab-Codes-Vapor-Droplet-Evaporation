function saveCTComparison(results,temp,z_o,t1,SavesFolder,neg_ana_C,error_rel_rms,count,error_rel2_rms,r_max);
        set(0,'DefaultFigureVisible', 'off');
%         fig = figure;
%         scatter(results(:,1),results(:,2))
%         hold on
%         scatter(results(:,1),temp)
%         xlim([0 50])
%         xlabel('r axis')
%         ylabel('C(r,z)')    
%         legend('Fitted Solution', 'Exact Solution')
%         title(['Concentration distribution comparison at z = ',num2str(z_o)])
%         saveas(fig,[SavesFolder,'/C_comparison.png'],'png')
%         saveas(fig,[SavesFolder,'/C_comparison.fig'],'fig')
        
         fig = figure;
        scatter(results(:,1),results(:,2),'o')
        hold on
        scatter(t1,temp,'x')
        if (count>0)
            scatter(neg_ana_C(:,1),neg_ana_C(:,2),'+','MarkerEdgeColor',[0 0 0],'LineWidth',1)
        end
        xlabel('r axis [mm]')
        ylabel('C(r,z) [mole/cm^3]')    
        if (count>0)
            legend('CT Output', 'Mod. Analyt. Sol','Org. Neg. Analyt. Sol')
        else
            legend('CT Output', 'Mod. Analyt. Sol')
        end
        title({['z = ',num2str(z_o),', Rel RMS (org.) = ',num2str(error_rel_rms),'%, Rel RMS (mod.) = ',num2str(error_rel2_rms),'%'];...
                    [' # points with neg analytical conc. = ',num2str(count)]});
        saveas(fig,[SavesFolder,'/C_full_comparison_z_',num2str(z_o),'.png'],'png')
        
        xlim([0 r_max])
        saveas(fig,[SavesFolder,'/C_short_comparison_z_',num2str(z_o),'.png'],'png')
        
        
        
        
        
end

