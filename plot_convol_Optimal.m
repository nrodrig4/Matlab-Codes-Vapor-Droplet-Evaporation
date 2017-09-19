function plot_convol_Optimal(OptimalSmoothness,sm,i,MinRMSTot, fsize,path3,Mu,lc,rho,rho_gridfit,rgrid,zgrid,tcoords_top,tcoords_side,...
                        fit_grad_t,fit_grad_ts,exact_grad_t,exact_grad_ts,...
                        minRR,maxRR,minZZ,maxZZ,radius,bdry_x,bdry_y)             

        if Mu*radius == MinRMSTot(i,1) && lc*radius == MinRMSTot(i,2) && sm == OptimalSmoothness(i)
        
            colorvt = [1.0 1.0 1.0];
            colorvs = [1.0 1.0 1.0];

            min_val = min(min([rho_gridfit;rho]));  %0
            max_val = max(max([rho_gridfit;rho]));  %2.0018e-01
            %min_val = 0;%min(min([rho_gridfit;rho]));
            %max_val = 6.3785e-01;%max(max([rho_gridfit;rho]));
            screen_size = get(0, 'ScreenSize');
            fig2=figure('Visible','off');
            set(fig2, 'Position', [0 0 screen_size(3) screen_size(4) ] );
            set(fig2,'PaperPositionMode','auto')
            axes('FontSize',fsize)
            hold all
            contourf(rgrid,zgrid,rho_gridfit)%,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
            caxis([min_val, max_val])
            h=colorbar;
            %alpha(0.5)
            view(2)

            plot3(tcoords_top(:,1),tcoords_top(:,2),zeros(size(tcoords_top(:,1))),'-r','Markersize',5,'Linewidth',4)
            plot3(tcoords_side(:,1),tcoords_side(:,2),zeros(size(tcoords_side(:,1))),'-r','Markersize',5,'Linewidth',4)
            plot3(-tcoords_top(:,1),tcoords_top(:,2),zeros(size(tcoords_top(:,1))),'-r','Markersize',5,'Linewidth',4)
            plot3(-tcoords_side(:,1),tcoords_side(:,2),zeros(size(tcoords_side(:,1))),'-r','Markersize',5,'Linewidth',4)

            quiver3(tcoords_top(:,1),tcoords_top(:,2),zeros(size(tcoords_top(:,1))),...
                -fit_grad_t(:,1),-fit_grad_t(:,2),zeros(size(tcoords_top(:,1))),...
                        'LineWidth',1,'Color',colorvt,'MaxHeadSize',1.0,'AutoScale','off');
            %size(fit_grad_t)
            %fit_grad_t

            quiver3(tcoords_side(:,1),tcoords_side(:,2),zeros(size(tcoords_side(:,1))),...
                -fit_grad_ts(:,1),-fit_grad_ts(:,2),zeros(size(tcoords_side(:,1))),...
                        'LineWidth',1,'Color',colorvs,'MaxHeadSize',1.0,'AutoScale','off');
            %size(fit_grad_ts)
            %fit_grad_ts
            plot3(bdry_x,bdry_y,zeros(size(bdry_x)),'--r','LineWidth',2)

            title(['Control volume of radius ',num2str(Mu*radius),' & height ',num2str(lc*radius),' (Computational)']);
            axis equal
            %axis([minRR maxRR minZZ maxZZ])
            axis([-4 4 0 3])
            xlabel('r*')
            ylabel('z*')
            saveas(fig2,[path3,'control_volume_gridfit_radius_',num2str(Mu*radius),...
             '_height_',num2str(lc*radius),'.png'],'png')


            %------------
            screen_size = get(0, 'ScreenSize');
            fig2=figure('Visible','off');
            set(fig2, 'Position', [0 0 screen_size(3) screen_size(4) ] );
            set(fig2,'PaperPositionMode','auto')
            axes('FontSize',fsize)
            hold all
            contourf(rgrid,zgrid,rho)%,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
            caxis([min_val, max_val])
            h=colorbar;
            %alpha(0.5)
            view(2)

            plot3(tcoords_top(:,1),tcoords_top(:,2),zeros(size(tcoords_top(:,1))),'-r','Markersize',5,'Linewidth',4)
            plot3(tcoords_side(:,1),tcoords_side(:,2),zeros(size(tcoords_side(:,1))),'-r','Markersize',5,'Linewidth',4)
            plot3(-tcoords_top(:,1),tcoords_top(:,2),zeros(size(tcoords_top(:,1))),'-r','Markersize',5,'Linewidth',4)
            plot3(-tcoords_side(:,1),tcoords_side(:,2),zeros(size(tcoords_side(:,1))),'-r','Markersize',5,'Linewidth',4)

            quiver3(tcoords_top(:,1),tcoords_top(:,2),zeros(size(tcoords_top(:,1))),...
                -exact_grad_t(:,1),-exact_grad_t(:,2),zeros(size(tcoords_top(:,1))),...
                        'LineWidth',1,'Color',colorvt,'MaxHeadSize',1.0,'AutoScale','off');
            %tcoords_top
            %size(exact_grad_t)
            %exact_grad_t

            quiver3(tcoords_side(:,1),tcoords_side(:,2),zeros(size(tcoords_side(:,1))),...
                -exact_grad_ts(:,1),-exact_grad_ts(:,2),zeros(size(tcoords_side(:,1))),...
                        'LineWidth',1,'Color',colorvs,'MaxHeadSize',1.0,'AutoScale','off');
            %tcoords_side
            %size(exact_grad_ts)
            %exact_grad_ts
            %pause

            plot3(bdry_x,bdry_y,zeros(size(bdry_x)),'--r','LineWidth',2)        
            title(['Control volume of radius ',num2str(Mu*radius),' & height ',num2str(lc*radius),' (Analytical)']);
            axis equal
            %axis([minRR maxRR minZZ maxZZ])
            axis([-4 4 0 3])
            xlabel('r*')
            ylabel('z*')
            saveas(fig2,[path3,'control_volume_theory_radius_',num2str(Mu*radius),...
             '_height_',num2str(lc*radius),'.png'],'png')


        end
end