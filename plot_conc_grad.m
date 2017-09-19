function plot_conc_grad(fsize,RR_org,ZZ_org,R_org,Z_org,rho_org_matrix,R,Z,rho,rho_gradr,rho_gradz,...
                    rgrid,zgrid,rho_gridfit,fit_gradr,fit_gradz,...
                    smooth_gf,hRR1,path1,path2)   

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%go along r
        %%%%%%%%%%%%%
        for ir = 1: length(RR_org)
            
            %%%%concentration%%%%%%%%%%%%%
            screen_size = get(0, 'ScreenSize');
            fig1=figure('Visible','off');
            set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
            axes('FontSize',fsize)
            set(fig1,'PaperPositionMode','auto')
            hold all

            %size(Z_org(:,ir)),size(rho_org_matrix(:,ir))
            plot(Z_org(:,ir),rho_org_matrix(:,ir),'db',...
                'Linewidth',4,'Markersize',15)
            
            test = R-RR_org(ir);
            [row,col] = find(abs(test) < 1e-5);
            
            %size(zgrid(row,col)),size(rho_gridfit(row,col))
            %temp1 = unique(zgrid(row,col),'rows','stable');
            %temp2 = unique(rho_gridfit(row,col),'rows','stable');
            temp1 = zgrid(row,col);
            temp2 = rho_gridfit(row,col);
            %size(temp1),size(temp2)
            plot(temp1(:,1),temp2(:,1),'or',...
                'Linewidth',2,'Markersize',8,'MarkerFace',[1 0.5 0]) 
            
            %temp1 = unique(Z(row,col),'rows','stable');
            %temp2 = unique(rho(row,col),'rows','stable');
            temp1 = Z(row,col);
            temp2 = rho(row,col);
            %size(temp1),size(temp2)
            plot(temp1(:,1),temp2(:,1),'-k',...
                'Linewidth',2,'Markersize',8) 
            legend('data','gridfit','diffusion-limited')
            xlabel('z')
            ylabel('C')
            title(['Smooth = ',num2str(smooth_gf),...
            ', grid spacing = ',num2str(hRR1),', r = ',num2str(RR_org(ir))]);
            saveas(fig1,[path1,'Conc_Smooth_',num2str(smooth_gf),...
            '_grid_spacing_',num2str(hRR1),'_r_',num2str(RR_org(ir)),'.png'],'png')
                
            %%%%partial C/partial r%%%%%%%%%%%%%
            screen_size = get(0, 'ScreenSize');
            fig2=figure('Visible','off');
            set(fig2, 'Position', [0 0 screen_size(3) screen_size(4) ] );
            axes('FontSize',fsize)
            set(fig2,'PaperPositionMode','auto')
            hold all

            temp1 = zgrid(row,col);
            temp2 = fit_gradr(row,col);
            plot(temp1(:,1),temp2(:,1),'or',...
                'Linewidth',2,'Markersize',8,'MarkerFace',[1 0.5 0]) 
            
            temp1 = Z(row,col);
            temp2 = rho_gradr(row,col);
            %size(temp1),size(temp2)
            plot(temp1(:,1),temp2(:,1),'-k',...
                'Linewidth',2,'Markersize',8) 
            
            legend('gridfit & gradient','diffusion-limited')
            xlabel('z')
            ylabel('\partial C/\partial r')
            title(['Smooth = ',num2str(smooth_gf),...
            ', grid spacing = ',num2str(hRR1),', r = ',num2str(RR_org(ir))]);
            saveas(fig2,[path1,'PartialC_partialr_Smooth_',num2str(smooth_gf),...
            '_grid_spacing_',num2str(hRR1),'_r_',num2str(RR_org(ir)),'.png'],'png')
    
            %%%%partial C/partial z%%%%%%%%%%%%%
            screen_size = get(0, 'ScreenSize');
            fig2=figure('Visible','off');
            set(fig2, 'Position', [0 0 screen_size(3) screen_size(4) ] );
            axes('FontSize',fsize)
            set(fig2,'PaperPositionMode','auto')
            hold all

            temp1 = zgrid(row,col);
            temp2 = fit_gradz(row,col);
            
            plot(temp1(:,1),temp2(:,1),'or',...
                'Linewidth',2,'Markersize',8,'MarkerFace',[1 0.5 0]) 
            
            temp1 = Z(row,col);
            temp2 = rho_gradz(row,col);
            plot(temp1(:,1),temp2(:,1),'-k',...
                'Linewidth',2,'Markersize',8) 
            
            legend('gridfit & gradient','diffusion-limited')
            xlabel('z')
            ylabel('\partial C/\partial z')
            title(['Smooth = ',num2str(smooth_gf),...
            ', grid spacing = ',num2str(hRR1),', r = ',num2str(RR_org(ir))]);
            saveas(fig2,[path1,'PartialC_partialz_Smooth_',num2str(smooth_gf),...
            '_grid_spacing_',num2str(hRR1),'_r_',num2str(RR_org(ir)),'.png'],'png')
        
            %pause
        end
        
        
        %%%%%%%%%%%%%
        %%%go along z 
        %%%%%%%%%%%%%
        for iz = 1: length(ZZ_org)
            screen_size = get(0, 'ScreenSize');
            fig2=figure('Visible','off');
            set(fig2, 'Position', [0 0 screen_size(3) screen_size(4) ] );
            axes('FontSize',fsize)
            set(fig2,'PaperPositionMode','auto')
            hold all

            %size(R_org(iz,:)),size(rho_org_matrix(iz,:))
            plot(R_org(iz,:),rho_org_matrix(iz,:),'db',...
                'Linewidth',4,'Markersize',15)
            
            test = Z-ZZ_org(iz);
            [row,col] = find(abs(test) < 1e-5);
            
            %size(rgrid(row,col)),size(rho_gridfit(row,col))
            temp1 = unique(rgrid(row,col),'rows','stable');
            temp2 = unique(rho_gridfit(row,col),'rows','stable');
            %size(temp1),size(temp2)
            plot(temp1(1,:),temp2(1,:),'or',...
                'Linewidth',2,'Markersize',8,'MarkerFace',[1 0.5 0]) 
 
            
            temp1 = unique(R(row,col),'rows','stable');
            temp2 = unique(rho(row,col),'rows','stable');
            %size(temp1),size(temp2)
            plot(temp1(1,:),temp2(1,:),'-k',...
                'Linewidth',2,'Markersize',8) 
            legend('data','gridfit','diffusion-limited')
            xlabel('r')
            ylabel('C')
            title(['Smooth = ',num2str(smooth_gf),...
            ', grid spacing = ',num2str(hRR1),', z = ',num2str(ZZ_org(iz))]);
            saveas(fig2,[path2,'Conc_Smooth_',num2str(smooth_gf),...
            '_grid_spacing_',num2str(hRR1),'_z_',num2str(ZZ_org(iz)),'.png'],'png')
        
            %%%%partial C/partial r%%%%%%%%%%%%%
            screen_size = get(0, 'ScreenSize');
            fig2=figure('Visible','off');
            set(fig2, 'Position', [0 0 screen_size(3) screen_size(4) ] );
            axes('FontSize',fsize)
            set(fig2,'PaperPositionMode','auto')
            hold all

            temp1 = rgrid(row,col);
            temp2 = fit_gradr(row,col);
            plot(temp1(1,:),temp2(1,:),'or',...
                'Linewidth',2,'Markersize',8,'MarkerFace',[1 0.5 0]) 
            
            temp1 = R(row,col);
            temp2 = rho_gradr(row,col);
            %size(temp1),size(temp2)
            plot(temp1(1,:),temp2(1,:),'-k',...
                'Linewidth',2,'Markersize',8) 
            
            legend('gridfit & gradient','diffusion-limited')
            xlabel('r')
            ylabel('\partial C/\partial r')
            title(['Smooth = ',num2str(smooth_gf),...
            ', grid spacing = ',num2str(hRR1),', z = ',num2str(ZZ_org(iz))]);
            saveas(fig2,[path2,'PartialC_partialr_Smooth_',num2str(smooth_gf),...
            '_grid_spacing_',num2str(hRR1),'_z_',num2str(ZZ_org(iz)),'.png'],'png')
    
            %%%%partial C/partial z%%%%%%%%%%%%%
            screen_size = get(0, 'ScreenSize');
            fig2=figure('Visible','off');
            set(fig2, 'Position', [0 0 screen_size(3) screen_size(4) ] );
            axes('FontSize',fsize)
            set(fig2,'PaperPositionMode','auto')
            hold all

            temp1 = rgrid(row,col);
            temp2 = fit_gradz(row,col);
            
            plot(temp1(1,:),temp2(1,:),'or',...
                'Linewidth',2,'Markersize',8,'MarkerFace',[1 0.5 0]) 
            
            temp1 = R(row,col);
            temp2 = rho_gradz(row,col);
            plot(temp1(1,:),temp2(1,:),'-k',...
                'Linewidth',2,'Markersize',8) 
            
            legend('gridfit & gradient','diffusion-limited')
            xlabel('r')
            ylabel('\partial C/\partial z')
            title(['Smooth = ',num2str(smooth_gf),...
            ', grid spacing = ',num2str(hRR1),', z = ',num2str(ZZ_org(iz))]);
            saveas(fig2,[path2,'PartialC_partialz_Smooth_',num2str(smooth_gf),...
            '_grid_spacing_',num2str(hRR1),'_z_',num2str(ZZ_org(iz)),'.png'],'png')
        
            %pause
        
        end
end