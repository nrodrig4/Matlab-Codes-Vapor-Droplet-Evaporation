function plot_conc_grad(fsize,RR_org,ZZ_org,R_org,Z_org,rho_org_matrix,R,Z,rho,rho_gradr,rho_gradz,...
                    rgrid,zgrid,rho_gridfit,fit_gradr,fit_gradz,...
                    smooth_gf,hRR1,path0,rmin,rmax,zmin,zmax,bdry_x,bdry_y)   

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%contour plots
        %%%%%%%%%%%%%concentration
        
        %myscale_C = [0.02:0.02:0.18];
        myscale_C = [0.00:0.02:0.2];
        miC = min(myscale_C);maC = max(myscale_C);
                
        %%%%%
        screen_size = get(0, 'ScreenSize');
        fig1=figure('Visible','off');
        set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        set(fig1,'PaperPositionMode','auto')
        axes('FontSize',fsize)
        
        surf(rgrid,zgrid,rho)
        hold on
        plot(bdry_x,bdry_y,'--r','LineWidth',2)

        xlabel('r*','FontSize',fsize)
        ylabel('z*','FontSize',fsize)
        zlabel('Concentration (mol/cm^3)','FontSize',fsize)
        title('Surface Fitted to Analytical Solution (Methanol)','FontSize',18)
        axis([-10,10,0,3.5,0,0.205])
        colorbar
        caxis([miC maC]);
        set(gca,'FontSize',fsize)
        saveas(fig1,[path0,'/1_Surface_Fitted_to_Analytical_Methanol_Data.png'],'png')

        %%%%%
        screen_size = get(0, 'ScreenSize');
        fig1=figure('Visible','off');
        set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        set(fig1,'PaperPositionMode','auto')
        axes('FontSize',fsize)
        
        surf(R_org,Z_org,rho_org_matrix)
        hold on
        plot(bdry_x,bdry_y,'--r','LineWidth',2)

        xlabel('r*','FontSize',fsize)
        ylabel('z*','FontSize',fsize)
        zlabel('Concentration (mol/cm^3)','FontSize',fsize)
        title('Surface Fitted to Analytical Data + Error','FontSize',18)
        axis([-10,10,0,3.5,0,0.205])
        colorbar
        caxis([miC maC]);
        set(gca,'FontSize',fsize)
        saveas(fig1,[path0,'/1_Surface_Fitted_to_Analytical_Methanol_Data_Plus_Error.png'],'png')

        %%%%%%%%
        screen_size = get(0, 'ScreenSize');
        fig1=figure('Visible','off');
        set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        set(fig1,'PaperPositionMode','auto')
        axes('FontSize',fsize)
        
        contourf(rgrid,zgrid,rho,myscale_C)
        hold on
        plot(bdry_x,bdry_y,'--r','LineWidth',2)
        line([-rmin rmin],[zmin zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmin -rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmin rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax rmax],[zmax zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax -rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmax rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        
        xlabel('r*','FontSize',fsize)
        ylabel('z*','FontSize',fsize)
        title('Filled 2D Contour Plot of Analytical Solution (Methanol)','FontSize',18)
        axis equal
        xlim([-5 5])
        ylim([0 3])
        
        caxis([miC maC]);
        colorbar
        set(gca,'FontSize',fsize)
        saveas(fig1,[path0,'/2_Filled_Contour_Fitted_to_Analytical_Methanol_Data.png'],'png')
        
        %%%%
        screen_size = get(0, 'ScreenSize');
        fig1=figure('Visible','off');
        set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        set(fig1,'PaperPositionMode','auto')
        axes('FontSize',fsize)
        
        contourf(R_org,Z_org,rho_org_matrix,myscale_C)
        hold on
        plot(bdry_x,bdry_y,'--r','LineWidth',2)
        line([-rmin rmin],[zmin zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmin -rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmin rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax rmax],[zmax zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax -rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmax rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        
        xlabel('r*','FontSize',fsize)
        ylabel('z*','FontSize',fsize)
        title('Filled 2D Contour Plot of Analytical Methanol Data + Error','FontSize',18)
        axis equal
        xlim([-5 5])
        ylim([0 3])
        
        caxis([miC maC]);
        colorbar
        set(gca,'FontSize',fsize)
        saveas(fig1,[path0,'/2_Filled_Contour_Fitted_to_Analytical_Methanol_Data_Plus_Error.png'],'png')
            
        
        %%%%%
        screen_size = get(0, 'ScreenSize');
        fig1=figure('Visible','off');
        set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        set(fig1,'PaperPositionMode','auto')
        axes('FontSize',fsize)
        
        surf(rgrid,zgrid,rho_gridfit)
        hold on
        plot(bdry_x,bdry_y,'--r','LineWidth',2)

        xlabel('r*','FontSize',fsize)
        ylabel('z*','FontSize',fsize)
        zlabel('Concentration (mol/cm^3)','FontSize',fsize)
        title('Surface Fitted to Computational Methanol Data by Gridfit','FontSize',18)
        axis([-10,10,0,3.5,0,0.205])
        caxis([miC maC]);        
        colorbar
        set(gca,'FontSize',fsize)
        saveas(fig1,[path0,'/3_Surface_Fitted_to_Computational_Methanol_Data_by_Gridfit.png'],'png')
        
        
        
        %%%%%
        screen_size = get(0, 'ScreenSize');
        fig1=figure('Visible','off');
        set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        set(fig1,'PaperPositionMode','auto')
        axes('FontSize',fsize)

        
        contourf(rgrid,zgrid,rho_gridfit,myscale_C)
        hold on
        plot(bdry_x,bdry_y,'--r','LineWidth',2)
        line([-rmin rmin],[zmin zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmin -rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmin rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax rmax],[zmax zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax -rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmax rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        xlabel('r*','FontSize',fsize)
        ylabel('z*','FontSize',fsize)
        title('Filled 2D Contour Plot of Computational Methanol Data by Gridfit','FontSize',18)
        axis equal
        xlim([-5 5])
        ylim([0 3])
        caxis([miC maC]);
        colorbar
        set(gca,'FontSize',fsize)

        saveas(fig1,[path0,'/4_Filled_Contour_Fitted_to_Computational_Methanol_Data_by_Gridfit.png'],'png')
       
        
        %%%%%%%%Error = Approximate - Exact  
        screen_size = get(0, 'ScreenSize');
        fig1=figure('Visible','off');
        set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        set(fig1,'PaperPositionMode','auto')
        axes('FontSize',fsize)
        
        error = (rho_gridfit-rho);
        %error = [-0.035:0.001:0.005];
        min_e = -0.035; %min(min(error));
        max_e = 0.005; %max(max(error));
        contourf(rgrid,zgrid,error)
        hold on
        plot(bdry_x,bdry_y,'--r','LineWidth',2)
        line([-rmin rmin],[zmin zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmin -rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmin rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax rmax],[zmax zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax -rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmax rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')

        xlabel('r*','FontSize',fsize)
        ylabel('z*','FontSize',fsize)
        title('Filled 2D Contour Plot of Error in Concentration (Computational - Analytical)','FontSize',18)
        axis equal
        xlim([-5 5])
        ylim([0 3])
        caxis([min_e max_e]);
        colorbar
        set(gca,'FontSize',fsize)
        saveas(fig1,[path0,'/5_Filled_Contour_Fitted_to_Methanol_Data_Error.png'],'png')
        
%         %%%%%%%%Error = (Approximate - Exact)/Exact  
%         screen_size = get(0, 'ScreenSize');
%         fig1=figure('Visible','off');
%         set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
%         set(fig1,'PaperPositionMode','auto')
%         axes('FontSize',fsize)
%         
%         error = error./rho*100;
%         min_e = min(min(error));
%         max_e = max(max(error));
%         contourf(rgrid,zgrid,error)
%         hold on
%         line([-1 1],[0 0],[0 0],'Linewidth',4,'LineStyle','-','Color','r')
%         xlabel('r*','FontSize',fsize)
%         ylabel('z*','FontSize',fsize)
%         %zlabel('Concentration (mol/cm^3)','FontSize',fsize)
%         title('Filled 2D Contour Plot of % Error in Concentration (Approx - Exact)/Exact*100','FontSize',18)
%         axis equal
%         xlim([-5 5])
%         ylim([0 3])
%         %caxis([min(rho_gridfit), max(rho_gridfit)])
%         caxis([min_e max_e]);
%         colorbar
%         set(gca,'FontSize',fsize)
%         saveas(fig1,[path0,'/Filled_Contour_Fitted_to_Analytical_Methanol_Data_Percent_Error.png'],'png')
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%partial C/partial r*
        myscale_Cr = [-0.5:0.01:0.5];
        %myscale_Cr = [0.0:0.02:0.14];
        miCr = min(myscale_Cr);maCr = max(myscale_Cr);
        
        [s1,s2] = size(fit_gradr);
        fit_gradr_tmp=fit_gradr(1:end,1:end); 
        rho_gradr_tmp=rho_gradr(1:end,1:end);
        
        minfit_gradr=min(min(fit_gradr_tmp));
        maxfit_gradr=max(max(fit_gradr_tmp));
        
        %%%%%
        screen_size = get(0, 'ScreenSize');
        fig1=figure('Visible','off');
        set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        set(fig1,'PaperPositionMode','auto')
        axes('FontSize',fsize)
        
        surf(rgrid(1:end,:),zgrid(1:end,:),rho_gradr_tmp(1:end,:))
        hold on
        plot(bdry_x,bdry_y,'--r','LineWidth',2)

        xlabel('r*','FontSize',fsize)
        ylabel('z*','FontSize',fsize)
        zlabel('\partial C/\partial r*','FontSize',fsize)
        title('Surface Fitted to Analytical \partial C/\partial r*','FontSize',18)
        axis([-10,10,0,3.5,minfit_gradr,maxfit_gradr])
        colorbar
        caxis([miCr maCr]);
        set(gca,'FontSize',fsize)
        saveas(fig1,[path0,'/6_Surface_Fitted_to_Analytical_parCparr.png'],'png')
        
        %%%%%%%%
        screen_size = get(0, 'ScreenSize');
        fig1=figure('Visible','off');
        set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        set(fig1,'PaperPositionMode','auto')
        axes('FontSize',fsize)
        
        
        contourf(rgrid(1:end,:),zgrid(1:end,:),rho_gradr_tmp(1:end,:),myscale_Cr)
        hold on
        plot(bdry_x,bdry_y,'--r','LineWidth',2)
        line([-rmin rmin],[zmin zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmin -rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmin rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax rmax],[zmax zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax -rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmax rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')

        xlabel('r*','FontSize',fsize)
        ylabel('z*','FontSize',fsize)
        title('Filled 2D Contour Plot of Analytical \partial C/\partial r*','FontSize',18)
        axis equal
        xlim([-5 5])
        ylim([0 3])
        colorbar
        caxis([miCr maCr]);
        set(gca,'FontSize',fsize)
        saveas(fig1,[path0,'/7_Filled_Contour_Fitted_to_Analytical_parCparr.png'],'png')

        %%%
        screen_size = get(0, 'ScreenSize');
        fig1=figure ('Visible','off');
        set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        set(fig1,'PaperPositionMode','auto')
        axes('FontSize',fsize)
        
        
        surf(rgrid,zgrid,fit_gradr_tmp)
        hold on
        plot(bdry_x,bdry_y,'--r','LineWidth',2)

        xlabel('r*','FontSize',fsize)
        ylabel('z*','FontSize',fsize)
        zlabel('\partial C/\partial r*','FontSize',fsize)
        title('Surface Fitted to \partial C/\partial r* from Gridfit Data','FontSize',18)
        axis([-10,10,0,3.5,minfit_gradr,maxfit_gradr])
        colorbar
        caxis([miCr maCr]);
        set(gca,'FontSize',fsize)
        saveas(fig1,[path0,'/8_Surface_Fitted_to_parCparr_by_Gridfit_Data.png'],'png')
        
        %%%%%
        screen_size = get(0, 'ScreenSize');
        fig1=figure('Visible','off');
        set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        set(fig1,'PaperPositionMode','auto')
        axes('FontSize',fsize)
        
        contourf(rgrid,zgrid,fit_gradr_tmp,myscale_Cr)
        hold on
        %plot(0,hRR1,'*','markersize',20)
        plot(bdry_x,bdry_y,'--r','LineWidth',2)

        line([-rmin rmin],[zmin zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmin -rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmin rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax rmax],[zmax zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax -rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmax rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')

        xlabel('r*','FontSize',fsize)
        ylabel('z*','FontSize',fsize)
        title('Filled 2D Contour Plot of \partial C/\partial r* from Gridfit Data','FontSize',18)
        axis equal
        xlim([-5 5])
        ylim([0 3])
        colorbar
        caxis([miCr maCr]);
        set(gca,'FontSize',fsize)
        saveas(fig1,[path0,'/9_Filled_Contour_Fitted_to_parCparr_by_Gridfit_Data.png'],'png')
        

        %%%%%%%%Error = Approximate - Exact
        screen_size = get(0, 'ScreenSize');
        fig1=figure('Visible','off');
        set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        set(fig1,'PaperPositionMode','auto')
        axes('FontSize',fsize)
        
        error = (fit_gradr_tmp(1:end,:)-rho_gradr_tmp(1:end,:));
        %error = [-0.35:0.001:0.35];
        min_e = -0.35; %min(min(error));
        max_e = 0.35; %max(max(error));
        contourf(rgrid(1:end,:),zgrid(1:end,:),error)
        hold on
        plot(bdry_x,bdry_y,'--r','LineWidth',2)

        line([-rmin rmin],[zmin zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmin -rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmin rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax rmax],[zmax zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax -rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmax rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')

        xlabel('r*','FontSize',fsize)
        ylabel('z*','FontSize',fsize)
        title('Filled 2D Contour Plot of Error in \partial C/\partial r* (Computational - Analytical)','FontSize',18)
        axis equal
        xlim([-5 5])
        ylim([0 3])
        colorbar
        caxis([min_e max_e]);
        set(gca,'FontSize',fsize)
        saveas(fig1,[path0,'/10_Filled_Contour_Fitted_to_parCparr_error.png'],'png')
        
%         %%%%%%%%Error = (Approximate - Exact)/Exact
%         screen_size = get(0, 'ScreenSize');
%         fig1=figure('Visible','off');
%         set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
%         set(fig1,'PaperPositionMode','auto')
%         axes('FontSize',fsize)
%         
%         error = error./rho_gradr_tmp(1:end,:)*100;
%         contourf(rgrid(1:end,:),zgrid(1:end,:),error)
%         hold on
%         line([-1 1],[0 0],[0 0],'Linewidth',4,'LineStyle','-','Color','r')
%         xlabel('r*','FontSize',fsize)
%         ylabel('z*','FontSize',fsize)
%         title('Filled 2D Contour Plot of % Error in \partial C/\partial r* (Approx - Exact)/Exact*100','FontSize',18)
%         axis equal
%         xlim([-5 5])
%         ylim([0 3])
%         %caxis([min(rho_gridfit), max(rho_gridfit)])
%         colorbar
%         %caxis([miCr maCr]);
%         set(gca,'FontSize',fsize)
%         saveas(fig1,[path0,'/Filled_Contour_Fitted_to_parCparr_Percent_error.png'],'png')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%partial C/partial z*
        myscale_Cz = [-0.5:0.01:0.0];
        miCz = min(myscale_Cz);maCz = max(myscale_Cz);
        
        minfit_gradz=min(min(fit_gradz));
        maxfit_gradz=max(max(fit_gradz));
        %%%%%
        screen_size = get(0, 'ScreenSize');
        fig1=figure('Visible','off');
        set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        set(fig1,'PaperPositionMode','auto')
        axes('FontSize',fsize)
        
        surf(rgrid,zgrid,rho_gradz)
        hold on
        plot(bdry_x,bdry_y,'--r','LineWidth',2)
        
        xlabel('r*','FontSize',fsize)
        ylabel('z*','FontSize',fsize)
        zlabel('\partial C/\partial z*','FontSize',fsize)
        title('Surface Fitted to Analytical \partial C/\partial z*','FontSize',18)
        axis([-10,10,0,3.5,minfit_gradz,maxfit_gradz])
        colorbar
        caxis([miCz maCz]);
        set(gca,'FontSize',fsize)
        saveas(fig1,[path0,'/11_Surface_Fitted_to_Analytical_parCparz.png'],'png')
        
        %%%%%%%%
        screen_size = get(0, 'ScreenSize');
        fig1=figure('Visible','off');
        set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        set(fig1,'PaperPositionMode','auto')
        axes('FontSize',fsize)
        
        contourf(rgrid,zgrid,rho_gradz,myscale_Cz)
        hold on
        plot(bdry_x,bdry_y,'--r','LineWidth',2)

        line([-rmin rmin],[zmin zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmin -rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmin rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax rmax],[zmax zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax -rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmax rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')

        xlabel('r*','FontSize',fsize)
        ylabel('z*','FontSize',fsize)
        title('Filled 2D Contour Plot of Analytical \partial C/\partial z*','FontSize',18)
        axis equal
        xlim([-5 5])
        ylim([0 3])
        colorbar
        caxis([miCz maCz]);
        set(gca,'FontSize',fsize)
        saveas(fig1,[path0,'/12_Filled_Contour_Fitted_to_Analytical_parCparz.png'],'png')
 
                
        %%%
        screen_size = get(0, 'ScreenSize');
        fig1=figure('Visible','off');
        set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        set(fig1,'PaperPositionMode','auto')
        axes('FontSize',fsize)
        
        surf(rgrid,zgrid,fit_gradz)
        hold on
        plot(bdry_x,bdry_y,'--r','LineWidth',2)
        
        xlabel('r*','FontSize',fsize)
        ylabel('z*','FontSize',fsize)
        zlabel('\partial C/\partial z*','FontSize',fsize)
        title('Surface Fitted to \partial C/\partial z* from Gridfit Data','FontSize',18)
        axis([-10,10,0,3.5,minfit_gradz,maxfit_gradz])
        colorbar
        caxis([miCz maCz]);
        set(gca,'FontSize',fsize)
        saveas(fig1,[path0,'/13_Surface_Fitted_to_parCparz_by_Gridfit_Data.png'],'png')
        
        %%%%%
        screen_size = get(0, 'ScreenSize');
        fig1=figure('Visible','off');
        set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        set(fig1,'PaperPositionMode','auto')
        axes('FontSize',fsize)
        
        contourf(rgrid,zgrid,fit_gradz,myscale_Cz)
        hold on
        plot(bdry_x,bdry_y,'--r','LineWidth',2)
        line([-rmin rmin],[zmin zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmin -rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmin rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax rmax],[zmax zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax -rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmax rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')

        xlabel('r*','FontSize',fsize)
        ylabel('z*','FontSize',fsize)
        %zlabel('Concentration (mol/cm^3)','FontSize',fsize)
        title('Filled 2D Contour Plot of \partial C/\partial z* from Gridfit Data','FontSize',18)
        axis equal
        xlim([-5 5])
        ylim([0 3])
        colorbar
        caxis([miCz maCz]);
        set(gca,'FontSize',fsize)
        saveas(fig1,[path0,'/14_Filled_Contour_Fitted_to_parCparz_by_Gridfit_Data.png'],'png')
        

        %%%%%%%%Error = Approximate - Exact
        screen_size = get(0, 'ScreenSize');
        fig1=figure('Visible','off');
        set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        set(fig1,'PaperPositionMode','auto')
        axes('FontSize',fsize)
        
        error = fit_gradz-rho_gradz;
        %error = [-0.15:0.01:0.35];
        min_e = -0.15; %min(min(error));
        max_e = 0.35; %max(max(error));
        contourf(rgrid,zgrid,error)
        hold on
        plot(bdry_x,bdry_y,'--r','LineWidth',2)

        line([-rmin rmin],[zmin zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmin -rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmin rmin],[0 zmin],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax rmax],[zmax zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([-rmax -rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')
        line([rmax rmax],[0 zmax],[0 0],'Linewidth',2,'LineStyle','-','Color','w')

        xlabel('r*','FontSize',fsize)
        ylabel('z*','FontSize',fsize)
        title('Filled 2D Contour Plot of Error in \partial C/\partial z* (Computational - Analytical)','FontSize',18)
        axis equal
        xlim([-5 5])
        ylim([0 3])
        colorbar
        caxis([min_e max_e]);
        set(gca,'FontSize',fsize)
        saveas(fig1,[path0,'/15_Filled_Contour_Fitted_to_parCparz_error.png'],'png')

%         %%%%%%%%Error = (Approximate - Exact)/Exact
%         screen_size = get(0, 'ScreenSize');
%         fig1=figure('Visible','off');
%         set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
%         set(fig1,'PaperPositionMode','auto')
%         axes('FontSize',fsize)
%         
%         error = error/rho_gradz*100;
%         contourf(rgrid,zgrid,error)
%         hold on
%         line([-1 1],[0 0],[0 0],'Linewidth',4,'LineStyle','-','Color','r')
%         xlabel('r*','FontSize',fsize)
%         ylabel('z*','FontSize',fsize)
%         title('Filled 2D Contour Plot of Error in \partial C/\partial z (Approx - Exact)/Exact*100','FontSize',18)
%         axis equal
%         xlim([-5 5])
%         ylim([0 3])
%         %caxis([min(rho_gridfit), max(rho_gridfit)])
%         colorbar
%         %caxis([miCz maCz]);
%         set(gca,'FontSize',fsize)
%         saveas(fig1,[path0,'/Filled_Contour_Fitted_to_parCparz_Percent_error.png'],'png')
        
end