function [trap_approx_int_pointwise, exact_int_pointwise, trap_approx_int, exact_int, tcoords, fit_grad_t, exact_grad_t,...
    trap_approx_int_pointwise_at_0, exact_int_pointwise_at_0, fit_C, exact_C] = ...
    top_integration(hRR1,Mu,lc,R,Z,rho_gridfit,fit_gradr,fit_gradz,...
                            rho,rho_gradr,rho_gradz,radius,fsize,D)
   
            %%%%RMS error will be calculated over all nodes
            %error_rms = 0;
                        
            %%%%integration on the top of the cylinder
            t=(0:hRR1:Mu)';
             
            tcoords=[t lc*ones(size(t))];
%             figure
%             axes('FontSize',fsize)
%             axis equal
%             plot(tcoords(:,1),tcoords(:,2),'or','Markersize',5,'Linewidth',1)
%             axis([0 tcoords(end,1) 0 tcoords(end,2)])
            
            fit_C = zeros(size(tcoords,1),1);
            
            exact_C = zeros(size(tcoords,1),1);
            
            fit_grad_t = zeros(size(tcoords,1),2);
            
            exact_grad_t = zeros(size(tcoords,1),2);
                        
            for i = 1:size(tcoords,1)
                rt=tcoords(i,1);
                zt=tcoords(i,2);
                
                testR = R - rt;
                [row1,col1] = find(abs(testR) < 1e-5);
                testZ = Z - zt;                
                [row2,col2] = find(abs(testZ) < 1e-5);
                %size(row1)
                %size(col1)
                %size(row2)
                %size(col2)
                [row,ia,ib] = intersect(row1,row2);
                [col,ia,ib] = intersect(col1,col2);
                %R(row,col)
                %Z(row,col)
                %pause
                
                 %row
                 %col
                 %size(fit_gradr)
                 %size(fit_gradz)
                fit_grad_t(i,1:2) = [fit_gradr(row,col) fit_gradz(row,col)];
                
                fit_C(i) = rho_gridfit(row,col);
      
                %%%theoretical gradient
                exact_grad_t(i,1:2) = [rho_gradr(row,col) rho_gradz(row,col)];
                
                exact_C(i) = rho(row,col);
                
                %temp = [rt zt rho_gridfit(row,col) rho(row,col)...
                %    fit_grad_t(i,1:2) exact_grad_t(i,1:2)];
                %dlmwrite(evap_file,temp,'delimiter',' ','precision','%e','-append');
                
                % RMS error between exact and fitted node gradients
                %gradr_rms = (exact_grad_t(i,1) - fit_grad_t(i,1))^2/2);
                %gradz_rms = (exact_grad_t(i,2) - fit_grad_t(i,2))^2;
                %error_rms = error_rms + gradz_rms;
            end

            normal_vec=repmat([0 1],length(t),1);

            %size(fit_grad_t)
            %size(normal_vec)
            func = sum(fit_grad_t.*normal_vec,2);
            trap_approx_int=D*2*pi*trapz(t,func.*t)*radius*10^(-3);
            
            trap_approx_int_pointwise=D*2*pi*func.*t*hRR1*radius*10^(-3);
            trap_approx_int_pointwise_at_0=D*2*pi*func(1)*0.5*hRR1*radius*10^(-3); %use this to count for 0 contribution due to "t" in the flux calculation at t = 0
            trap_approx_int_pointwise(1)=0.5*trap_approx_int_pointwise(1);
            trap_approx_int_pointwise(end)=0.5*trap_approx_int_pointwise(end);
            
            %%%%%%%%%%%in theory (using theoretical gradients)
            func = sum(exact_grad_t.*normal_vec,2);
            exact_int=D*2*pi*trapz(t,func.*t)*radius*10^(-3);
            
            exact_int_pointwise=D*2*pi*func.*t*hRR1*radius*10^(-3);
            exact_int_pointwise_at_0=D*2*pi*func(1)*0.5*hRR1*radius*10^(-3);
            exact_int_pointwise(1)=0.5*exact_int_pointwise(1);
            exact_int_pointwise(end)=0.5*exact_int_pointwise(end);
            
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%% concentration
%             screen_size = get(0, 'ScreenSize');
%             fig2=figure('Visible','off');
%             set(fig2, 'Position', [0 0 screen_size(3) screen_size(4) ] );
%             set(fig2,'PaperPositionMode','auto')
%             axes('FontSize',fsize)
%             hold all
%             
%             plot(tcoords(:,1),exact_C,'-k',...
%                 'Linewidth',2,'Markersize',8)      
%             plot(tcoords(:,1),fit_C,'or',...
%                 'Linewidth',2,'Markersize',8,'MarkerFace',[1 0.5 0])
%             legend('exact','gridfit')
%             xlabel('r')
%             ylabel('C')
%             title(['radius ',num2str(Mu),', height ',num2str(lc)])
%             saveas(fig2,[path,'C_radius_',num2str(Mu),...
%              '_height_',num2str(lc),'_top','.png'],'png')
% 
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%% partial derivatives
%             screen_size = get(0, 'ScreenSize');
%             fig2=figure('Visible','off');
%             set(fig2, 'Position', [0 0 screen_size(3) screen_size(4) ] );
%             set(fig2,'PaperPositionMode','auto')
%             axes('FontSize',fsize)
%             hold all
%             
%             plot(tcoords(:,1),exact_grad_t(:,1),'-k',...
%                 'Linewidth',2,'Markersize',8)       
%             plot(tcoords(:,1),fit_grad_t(:,1),'or',...
%                 'Linewidth',2,'Markersize',8,'MarkerFace',[1 0.5 0])
%             xlabel('r')
%             ylabel('\partial C / \partial r')
%             legend('exact gradient','gridfit & gradient')
%             title(['radius ',num2str(Mu),', height ',num2str(lc)])
%             saveas(fig2,[path,'parC_parr_radius_',num2str(Mu),...
%              '_height_',num2str(lc),'_top','.png'],'png')
% 
%             screen_size = get(0, 'ScreenSize');
%             fig2=figure('Visible','off');
%             set(fig2, 'Position', [0 0 screen_size(3) screen_size(4) ] );
%             set(fig2,'PaperPositionMode','auto')
%             axes('FontSize',fsize)
%             hold all
%             
%             plot(tcoords(:,1),exact_grad_t(:,2),'-k',...
%                 'Linewidth',2,'Markersize',8)       
%             plot(tcoords(:,1),fit_grad_t(:,2),'or',...
%                 'Linewidth',2,'Markersize',8,'MarkerFace',[1 0.5 0])
%             xlabel('r')
%             ylabel('\partial C / \partial z')
%             legend('exact gradient','gridfit')
%             title(['radius ',num2str(Mu),', height ',num2str(lc)])
%             saveas(fig2,[path,'parC_parz_radius_',num2str(Mu),...
%              '_height_',num2str(lc),'_top','.png'],'png')         


end
