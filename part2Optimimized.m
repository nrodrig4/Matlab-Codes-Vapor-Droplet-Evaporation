function part2Optimimized(app,D,Rad,compound_theory,Cf,C_v,rootdir,resultpos,sm_val,noise_flag)
        
            RMSVals8x2 = zeros(length(sm_val),1);
            RMSVals8x4 = zeros(length(sm_val),1);
            RMSVals16x2 = zeros(length(sm_val),1);
            RMSVals16x4 = zeros(length(sm_val),1);    
    for smooth = 1:length(sm_val)

            
            sm = sm_val(smooth);
            smoothnessText = ['Smoothness ',num2str(sm), ' of ', num2str(sm_val(length(sm_val)))];
            clc;
            close all;
            clearvars -except RMSVals8x2 RMSVals16x2 RMSVals8x4 RMSVals16x4 RMSVals app sm C_v Cf smooth resultpos sm_val rootdir Rad compound_theory  D noise_flag smoothnessText;
            format shortE;
%             path0=[rootdir,'/Contour_Plots_Scaled_r_',num2str(scaler),'_Scaled_z_',num2str(scalez),'_Noise_',num2str(noise_flag),'_Smooth_',num2str(sm(smm)),'_grid_spacing_',num2str(hRR1*radius),'/'];        
%             path3=[rootdir,'/Evaporation_Rate_Scaled_r_',num2str(scaler),'_Scaled_z_',num2str(scalez),'_Noise_',num2str(noise_flag),'_Smooth_',num2str(sm(smm)),'_grid_spacing_',num2str(hRR1*radius),'/'];

            fsize = 24;

            %%%number of replicates/tests
            num_test = 1;
            %num_test = 5;

            %%%scale the concentration along the r-direction and z-direction
            scaler = 1;  %no scale
            scalez = 1;  %no scale

            %scaler = 0.5;
            %scalez = 2;   %scale > 1: compressed; scale < 1: stretched

            %%%add noise to IA (NOT TO concentration)
            %noise_IA = 1.5168e-3; %kg/m^3 about (1.5168e-3/2.0018e-01)*100 = 0.75772%, i.e., the "error" noise is about 0.7% of the saturated concentration, Cv

            noise = 0;%%%%DO NOT CHANGE THIS even when adding noise to IA

            radius=Rad;   % % this is the radius of the drop (mm)
            org_cell_h=0.5;
            hRR1=org_cell_h;
            org_hZZ1=1.0;  

            RRo=(0:hRR1:35)/radius; % r vector
            ZZ1=1/radius;
            ZZ2=6/radius;
            ZZ3=10/radius;
            ZZ4=20/radius;

            rmin = 8.0; %min r of control volumes
            rmax = 16.0; %max r of control volumes
            zmin = 2.0;   %min z of control volumes
            zmax = 4.0; %max z of control volumes

            %Weber's disc noundary(normalized) 
            bdry_x=[-1,-0.974927912181824,-0.930873748644204,-0.866025403784439,-0.781831482468030,-0.680172737770919,-0.563320058063622,-0.433883739117558,-0.294755174410904,-0.149042266176175,0,0.149042266176175,0.294755174410904,0.433883739117558,0.563320058063622,0.680172737770919,0.781831482468030,0.866025403784439,0.930873748644204,0.974927912181824,1];
            bdry_y=zeros(size(bdry_x));

            h=1;   %3;  %number of grid refinement level: h = 1 means starting out with grid spacing of hRR1


            %% The data is organized and ordered in this section
            RR_orgo=RRo;
            hRR1o=hRR1/radius;
            minRRo=min(RRo);
            maxRRo=max(RRo);

            hZZ1o=org_hZZ1/radius;
            hZZ2o=2*org_hZZ1/radius;
            hZZ3o=5*org_hZZ1/radius;

            ZZ1a=ZZ1:hZZ1o:ZZ2;
            ZZ2a=ZZ2+hZZ2o:hZZ2o:ZZ3;
            ZZo=[ZZ1a ZZ2a ZZ3+hZZ3o:hZZ3o:ZZ4]% z vector
            ZZ_org=ZZo;


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %This rearranges ZZ and RR according to gridfit format
            [Ro,Zo]=meshgrid(RRo,ZZo);
            R_org = Ro;
            Z_org = Zo;
            m=length(RRo);
            n=length(ZZo);

            r_org=reshape(Ro',n*m,1)';% original r vector
            z_org=reshape(Zo',n*m,1)';% original z vector

            for test_indx = 1: num_test
                disp('------------------Test Number---------------')
                test_indx

                RR = RRo;
                RR_org = RR_orgo;
                hRR1 = hRR1o;
                minRR = minRRo;
                maxRR = maxRRo;
                hZZ1 = hZZ1o;
                hZZ2=hZZ2o;
                hZZ3=hZZ3o;
                ZZ = ZZo;

                %% This creates an analytical vapor distribution with the CT output gri

                    %load([rootdir,'/resultpos_with_polyfit.mat']);


                rho_org =resultpos(1:71,:);
                %size(rho_org)
                %size(RR)
                %pause

                rho_org_matrix = rho_org';

                disp('Max and Min values of imposed concentration values BEFORE adding noise')
                max(max(rho_org_matrix))
                min(min(rho_org_matrix))
                disp('Max and Min values of imposed concentration values AFTER adding noise')
                rho_org_matrix = rho_org_matrix+((2.0*rand(size(rho_org_matrix))-1)*noise);
                max(max(rho_org_matrix))
                min(min(rho_org_matrix))

                %%%mirror data over z-axis and r-axis%%%%%%%%%%%%%%
                RR_org=unique([-RR_org RR_org]);
                minRR=min(RR_org);
                maxRR=max(RR_org);
                [R,Z]=meshgrid(RR_org,ZZ_org);
                R_org = R;
                Z_org = Z;
                m=length(RR_org);
                n=length(ZZ_org);
                r_org=reshape(R',n*m,1)';% mirrored r vector
                z_org=reshape(Z',n*m,1)';% original z vector
                size(rho_org_matrix)  %10   149
                rho_org_matrix = [fliplr(rho_org_matrix) rho_org_matrix(:,2:end)];
                rho_org=reshape(rho_org_matrix',n*m,1)';    

                %%%flip about r-axis
                ZZ_org_tmp=[fliplr(-ZZ_org) ZZ_org];
                [R_tmp,Z_tmp]=meshgrid(RR_org,ZZ_org_tmp);
                R_org_tmp = R_tmp;
                Z_org_tmp = Z_tmp;
                m_tmp=length(RR_org);
                n_tmp=length(ZZ_org_tmp);
                r_org_tmp=reshape(R_tmp',n_tmp*m_tmp,1)';% mirrored r vector
                z_org_tmp=reshape(Z_tmp',n_tmp*m_tmp,1)';% original z vector


                %%%%add saturated concentration along droplet boundary
                %r_org_tmp=[r_org_tmp,bdry_x];
                %z_org_tmp=[z_org_tmp,bdry_y];

                size(rho_org_matrix);  %10   297
                rho_org_matrix_tmp = [flipud(rho_org_matrix);rho_org_matrix];
                size(rho_org_matrix_tmp);  %20   297

                rho_org_tmp=reshape(rho_org_matrix_tmp',n_tmp*m_tmp,1)';
                size(rho_org_tmp) ;  %1        5940

                %%%add saturated concentration along the droplet boundary
                %rho_org_tmp = [rho_org_tmp,Cv*ones(size(bdry_x))];

                hRR1_org = hRR1;
                hZZ1_org = hZZ1;
                hZZ2_org = hZZ2;
                hZZ3_org = hZZ3;
                minRR_org = minRR;
                maxRR_org = maxRR;

                min_percent_error_cv_master = [];
                min_rms_error_cv_master = [];
                for smm=1:length(sm)

                    disp('*********Smoothness**********')
                    smooth_gf = sm(smm) 

                    for j=1:h

                        %disp('Normalized grid spacing: ')
                        if j ~= 1
                            hRR1=hRR1/2;% relative cell size in radial direction
                            hZZ1=hZZ1/2;
                            hZZ2=hZZ2/2;
                            hZZ3=hZZ3/2;
                        else
                            hRR1=hRR1_org;
                            hZZ1=hZZ1_org;
                            hZZ2=hZZ2_org;
                            hZZ3=hZZ3_org; 
                            minRR = minRR_org;
                            maxRR = maxRR_org;
                        end
                        disp('Grid size in mm ')
                        hRR1*radius

%                        path0=[rootdir,'/Contour_Plots_Scaled_r_',num2str(scaler),'_Scaled_z_',num2str(scalez),'_Noise_',num2str(noise_flag),'_Smooth_',num2str(sm(smm)),'_grid_spacing_',num2str(hRR1*radius),'/'];        
%                        path3=[rootdir,'/Evaporation_Rate_Scaled_r_',num2str(scaler),'_Scaled_z_',num2str(scalez),'_Noise_',num2str(noise_flag),'_Smooth_',num2str(sm(smm)),'_grid_spacing_',num2str(hRR1*radius),'/'];

%                         if (exist(path0)~=0)
%                             rmdir(path0,'s') 
%                         end
%                         mkdir(path0)
% 
%                         if (exist(path3)~=0)
%                             rmdir(path3,'s') 
%                         end
%                         mkdir(path3)       

                        if (smm > 1) || (j > 1)
                            clear RR ZZ rho rho_gradr rho_gradz fit rgrid zgrid dif fitt rhot rgridt zgridt A max_num max_idx row col...
                                fit_gradr fit_gradz;     
                        end

                        %create grid for gridfit
                        RR=minRR-hRR1:hRR1:maxRR+hRR1; % radial vector
                        ZZ=(-0.0*hRR1:hRR1:ZZ4+hRR1); % z vector 

                        minRR = min(RR);
                        maxRR = max(RR);

                        minZZ = min(ZZ);
                        maxZZ = max(ZZ);        

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %analytic solutions for concentration and derivatives
                        [rho,rho_gradr,rho_gradz] = vap_dist_modification_generalization_sub_C50(RR,scaler,ZZ,scalez,C_v,Cf);
                        rho=rho';  %rows: along z, columns: along r
                        rho_gradz=rho_gradz';
                        rho_gradr=rho_gradr';       

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                        %using gridfit: http://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit
                        [rho_gridfit,rgrid,zgrid]=gridfit(r_org,z_org,rho_org,RR,ZZ,'smooth',smooth_gf); %,'tilesize',100);
                        [rho_gridfit_tmp,rgrid_tmp,zgrid_tmp]=gridfit(r_org_tmp,z_org_tmp,rho_org_tmp,RR,unique([fliplr(-ZZ) ZZ]),'smooth',smooth_gf); %,'tilesize',100);
                        rho_gridfit=rho_gridfit_tmp(round(size(rho_gridfit_tmp,1)/2):end,:);
                        zgrid=zgrid_tmp(round(size(zgrid_tmp,1)/2):end,:);

                        %%%set rho_gridfit to zero if negative
                        [row_sz,col_sz] = size(rho_gridfit);
                        %[row_sz,col_sz,neg_rho_gridfit] = find(rho_gridfit<0.0);
                        count = 0;
                        for ii = 1:row_sz
                            for jj = 1:col_sz
                                if rho_gridfit(ii,jj) < 0
                                   count = count + 1;
                                   rho_gridfit(ii,jj) = 0; 
                                end
                            end
                        end
                        %disp('Number of Negative Concentrations from gridfit')
                        %count


                        [fit_gradr,fit_gradz]=gradient(rho_gridfit,hRR1,hRR1);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        [R,Z]=meshgrid(RR,ZZ);

                       %%%%%%%compare gridfit results with theoretical results by visualization
                       %%%regarding concentration and gradient on the uniform grid
                       %plot_conc_grad2(fsize,RR_org,ZZ_org,R_org,Z_org,rho_org_matrix,R,Z,rho,rho_gradr,rho_gradz,...
                                  % rgrid,zgrid,rho_gridfit,fit_gradr,fit_gradz,...
                                  % smooth_gf,hRR1,path0,rmin/radius,rmax/radius,zmin/radius,zmax/radius,...
                                  % bdry_x,bdry_y);
                        %pause
                        %%%compare concentration and gradient along r and z slices
                %         plot_conc_grad(fsize,RR_org,ZZ_org,R_org,Z_org,rho_org_matrix,R,Z,rho,rho_gradr,rho_gradz,...
                %                     rgrid,zgrid,rho_gridfit,fit_gradr,fit_gradz,...
                %                     smooth_gf,hRR1,path1,path2) 


                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %define control volume & calculate evaporation rate
                        Mu_array = (rmin:8.0:rmax)/radius;
                        %Mu_array = (7.0:0.5:7.5)/radius;
                        Mu_len = length(Mu_array);
                        %pause

                        height_array = (zmin:2.0:zmax)/radius; %assume that the height of the droplet is less than 1 mm.
                        %height_array = (1.0:0.5:1.5)/radius;
                        Height_len = length(height_array);

                        number_of_control_volume = Height_len*Mu_len;
                        %pause

                        test_array = zeros(number_of_control_volume,5);
                        for lc_ind = 1:length(height_array)

                            lc_array = height_array(lc_ind)*ones(size(Mu_array));  

                            trap_approx_int_cv = zeros(size(Mu_array));

                            exact_int_cv = zeros(size(Mu_array));

                            %%%%RMS error will be calculated for each control volume
                            error_rms = zeros(size(Mu_array));
                            avg_per_error = zeros(size(Mu_array));
                            for check=1:length(Mu_array)
                                Mu = Mu_array(check);
                                lc = lc_array(check);

%                                 evap_file1 = [path3,'Top_Evaporation_rate_height_',num2str(lc),'_radius_',num2str(Mu)];
%                                 evap_file2 = [path3,'Side_Evaporation_rate_height_',num2str(lc),'_radius_',num2str(Mu)];

                                %%%top integral
                                %disp('top')
    %                             [trap_approx_int_pointwise,exact_int_pointwise,trap_approx_int,exact_int,tcoords_top,fit_grad_t,exact_grad_t,...
    %                                 trap_approx_int_pointwise_at_0,exact_int_pointwise_at_0,fit_Ct,exact_Ct]
                                [trap_approx_int_pointwise, exact_int_pointwise, trap_approx_int, exact_int, tcoords_top, fit_grad_t, exact_grad_t,...
                                 trap_approx_int_pointwise_at_0, exact_int_pointwise_at_0, fit_C, exact_Ct] =top_integration(hRR1,Mu,lc,R,Z,rho_gridfit,fit_gradr,fit_gradz,...
                                                rho,rho_gradr,rho_gradz,radius,fsize,D);


                                %%%lateral side integral
                                %disp('side')
                                [trap_approx_int_pointwise_cv,exact_int_pointwise_cv,trap_approx_int_cv(check),exact_int_cv(check),tcoords_side,fit_grad_ts,exact_grad_ts,fit_Cs,exact_Cs] =...
                                    side_integration(hRR1,Mu,lc,R,Z,rho_gridfit,fit_gradr,fit_gradz,...
                                                rho,rho_gradr,rho_gradz,radius,fsize,D);


                                %%%%%%%%%%%%%%%%%%%%%
                                trap_approx_int_cv(check) = trap_approx_int_cv(check) + trap_approx_int;

                                exact_int_cv(check) = exact_int_cv(check) + exact_int;            


                               %%%plot control volume
                              % plot_convol(fsize,path0,Mu,lc,rho,rho_gridfit,rgrid,zgrid,tcoords_top,tcoords_side,...
                                       %fit_grad_t,fit_grad_ts,exact_grad_t,exact_grad_ts,...
                                       %minRR,maxRR,minZZ,maxZZ,radius,bdry_x,bdry_y)
                                %pause

                                trap_approx_int_pointwise=trap_approx_int_pointwise(2:end);
                                exact_int_pointwise=exact_int_pointwise(2:end);
                                %abs discretization error
                                top_error_rms = [(trap_approx_int_pointwise_at_0-exact_int_pointwise_at_0)/exact_int_pointwise_at_0;...
                                                 (trap_approx_int_pointwise-exact_int_pointwise)./exact_int_pointwise];
                                indt = find(exact_Ct > -eps);
                                top_error_rms_reduced=top_error_rms(indt); %only count flux components at points with non-negative concentrations

                                side_error_rms = (trap_approx_int_pointwise_cv-exact_int_pointwise_cv)./exact_int_pointwise_cv;
                                inds = find(exact_Cs > -eps);
                                side_error_rms_reduced=side_error_rms(inds); %only count flux components at points with non-negative concentrations

                                cv_error_rms = [top_error_rms_reduced;side_error_rms_reduced];

                                error_rms(check) = 100*sqrt(mean(cv_error_rms.^2));
                                %sqrt(sum(cv_error_rms.^2)/sum([exact_int_pointwise;exact_int_pointwise_cv].^2));

                                avg_per_error(check) = mean(abs(cv_error_rms));
                            end   

                            test = [radius*Mu_array' radius*lc_array' -trap_approx_int_cv' -exact_int_cv'];

                            theory_no_trap = [Mu_array' lc_array' exact_int_cv'];

                            %%%%compute percent error from grid_fit and "gradient"
                            theory = compound_theory;
                            sol = abs(trap_approx_int_cv');
                            percent_error = 100*(sol - theory*ones(size(sol)))./(theory*ones(size(sol)));

                            sol_theo = abs(exact_int_cv');
                            percent_error_theo = 100*(sol - sol_theo)./(sol_theo);

                            test = [test percent_error_theo percent_error avg_per_error' error_rms'];

                           % evap_file = [path3,'Evaporation_Rate_Smooth_',num2str(sm(smm)),'_grid_spacing_',num2str(hRR1*radius)];
                           % dlmwrite(evap_file,test,'delimiter',' ','precision','%e','-append');
                            test_array((lc_ind-1)*Mu_len+1:lc_ind*Mu_len,1:8) = test;

                        end
                        %add control volume with the smallest percent error to the master list
                        %min_percent_error_cv_master = [min_percent_error_cv_master;sm(smm),hRR1*radius,test_sorted(1,:)];

                        %sort the rms errors in ascending orders
                        [~,idx]=sort(test_array(:,end));

                        test_sorted=test_array(idx,:);
                        for i = 1:length(test_sorted(:,1))
           
                        TextAreaText1{i} = [num2str(test_sorted(i,1)),'x',num2str(test_sorted(i,2)),', RMS = ',num2str(test_sorted(i,length(test_sorted(1,:)))),'%'];
                            switch(i)
                                case 1
                                    RMSVals16x4(smooth) = test_sorted(i,length(test_sorted(1,:)));
                                    
                                case 2
                                    RMSVals8x4(smooth) = test_sorted(i,length(test_sorted(1,:)));
                                    
                                case 3
                                    RMSVals8x2(smooth) = test_sorted(i,length(test_sorted(1,:)));
                                    
                                case 4
                                    RMSVals16x2(smooth) = test_sorted(i,length(test_sorted(1,:)));
                                    
                            end
                        end

                        RMSarray(:,1) = test_array(:,1);
                        RMSarray(:,2) = test_array(:,2);

                        RMSVals(:,smooth) = test_array(:,length(test_array));
                        
                        

                        [MinRMS1,Index1] = min(RMSVals(1,:));
                        [MinRMS2,Index2] = min(RMSVals(2,:));
                        [MinRMS3,Index3] = min(RMSVals(3,:));
                        [MinRMS4,Index4] = min(RMSVals(4,:));
                        MinRMS1 = [MinRMS1 sm_val(Index1)];
                        MinRMS2 = [MinRMS2 sm_val(Index2)];
                        MinRMS3 = [MinRMS3 sm_val(Index3)];
                        MinRMS4 = [MinRMS4 sm_val(Index4)];
                        MinRMS = [MinRMS1; MinRMS2; MinRMS3; MinRMS4];
                        MinRMSTot(:,1) = RMSarray(:,1);
                        MinRMSTot(:,2) = RMSarray(:,2);
                        MinRMSTot = [MinRMSTot MinRMS]
                        OptimalSmoothness = MinRMS(:,2);
                        
                        for i = 1:length(RMSarray(:,1))
                            RMSText{i} = [num2str(RMSarray(i,1)),'x',num2str(RMSarray(i,2)),', RMS = ',num2str(MinRMS(i,1)),' at ',num2str(MinRMS(i,2)),' smoothness'];
                        end
                        averageRMS = ['Average RMS = ',num2str(sum(MinRMS(:,1))/length(test_sorted(:,length(test_sorted)))),'%']
                        %CVFluxText = ['Min RMS for each CV';]
                        app.CVFluxRMSTextArea.Value = [smoothnessText;RMSText';averageRMS];

                       % evap_file = [path3,'RMS_Sorted_Evaporation_Rate_Smooth_',num2str(sm(smm)),'_grid_spacing_',num2str(hRR1*radius)];
                       % dlmwrite(evap_file,test_sorted,'delimiter',' ','precision','%e');

                        %%add control volume with the smallest rms error to the master list
                        %min_rms_error_cv_master = [min_rms_error_cv_master;sm(smm),hRR1*radius,test_sorted(1,:)];
                    end
                end
            end
    end
    data8x2(:,1) = sm_val;
    data8x2(:,2) = RMSVals8x2;
    saveRMSData(data8x2,rootdir,'8x2')
    
    data16x2(:,1) = sm_val;
    data16x2(:,2) = RMSVals16x2;
    saveRMSData(data16x2,rootdir,'16x2')
    
    data8x4(:,1) = sm_val;
    data8x4(:,2) = RMSVals8x4;
    saveRMSData(data8x4,rootdir,'8x4')
    
    data16x4(:,1) = sm_val;
    data16x4(:,2) = RMSVals16x4;
    saveRMSData(data16x4,rootdir,'16x4')
    
    for i = 1:length(OptimalSmoothness)
        sm = OptimalSmoothness(i);
        
            clc;
            close all;
            clearvars -except RMSVals app sm C_v Cf smooth resultpos sm_val MinRMSTot rootdir Rad...
                        compound_theory  D noise_flag OptimalSmoothness RMSText averageRMS i;
            format shortE;
            ControlVolumesLog = ['Producing images for control Volumes (',num2str(i),'/',num2str(length(OptimalSmoothness)),')'];
            app.CVFluxRMSTextArea.Value = [RMSText';averageRMS;ControlVolumesLog];

            fsize = 24;

            %%%number of replicates/tests
            num_test = 1;
            %num_test = 5;

            %%%scale the concentration along the r-direction and z-direction
            scaler = 1;  %no scale
            scalez = 1;  %no scale

            %scaler = 0.5;
            %scalez = 2;   %scale > 1: compressed; scale < 1: stretched

            %%%add noise to IA (NOT TO concentration)
            %noise_IA = 1.5168e-3; %kg/m^3 about (1.5168e-3/2.0018e-01)*100 = 0.75772%, i.e., the "error" noise is about 0.7% of the saturated concentration, Cv

            noise = 0;%%%%DO NOT CHANGE THIS even when adding noise to IA

            radius=Rad;   % % this is the radius of the drop (mm)
            org_cell_h=0.5;
            hRR1=org_cell_h;
            org_hZZ1=1.0;  

            RRo=(0:hRR1:35)/radius; % r vector
            ZZ1=1/radius;
            ZZ2=6/radius;
            ZZ3=10/radius;
            ZZ4=20/radius;

            rmin = 8.0; %min r of control volumes
            rmax = 16.0; %max r of control volumes
            zmin = 2.0;   %min z of control volumes
            zmax = 4.0; %max z of control volumes

            %Weber's disc noundary(normalized) 
            bdry_x=[-1,-0.974927912181824,-0.930873748644204,-0.866025403784439,-0.781831482468030,-0.680172737770919,-0.563320058063622,-0.433883739117558,-0.294755174410904,-0.149042266176175,0,0.149042266176175,0.294755174410904,0.433883739117558,0.563320058063622,0.680172737770919,0.781831482468030,0.866025403784439,0.930873748644204,0.974927912181824,1];
            bdry_y=zeros(size(bdry_x));

            h=1;   %3;  %number of grid refinement level: h = 1 means starting out with grid spacing of hRR1


            %% The data is organized and ordered in this section
            RR_orgo=RRo;
            hRR1o=hRR1/radius;
            minRRo=min(RRo);
            maxRRo=max(RRo);

            hZZ1o=org_hZZ1/radius;
            hZZ2o=2*org_hZZ1/radius;
            hZZ3o=5*org_hZZ1/radius;

            ZZ1a=ZZ1:hZZ1o:ZZ2;
            ZZ2a=ZZ2+hZZ2o:hZZ2o:ZZ3;
            ZZo=[ZZ1a ZZ2a ZZ3+hZZ3o:hZZ3o:ZZ4]% z vector
            ZZ_org=ZZo;


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %This rearranges ZZ and RR according to gridfit format
            [Ro,Zo]=meshgrid(RRo,ZZo);
            R_org = Ro;
            Z_org = Zo;
            m=length(RRo);
            n=length(ZZo);

            r_org=reshape(Ro',n*m,1)';% original r vector
            z_org=reshape(Zo',n*m,1)';% original z vector

            for test_indx = 1: num_test
                disp('------------------Test Number---------------')
                test_indx

                RR = RRo;
                RR_org = RR_orgo;
                hRR1 = hRR1o;
                minRR = minRRo;
                maxRR = maxRRo;
                hZZ1 = hZZ1o;
                hZZ2=hZZ2o;
                hZZ3=hZZ3o;
                ZZ = ZZo;

                %% This creates an analytical vapor distribution with the CT output gri

                    %load([rootdir,'/resultpos_with_polyfit.mat']);


                rho_org =resultpos(1:71,:);
                %size(rho_org)
                %size(RR)
                %pause

                rho_org_matrix = rho_org';

                disp('Max and Min values of imposed concentration values BEFORE adding noise')
                max(max(rho_org_matrix))
                min(min(rho_org_matrix))
                disp('Max and Min values of imposed concentration values AFTER adding noise')
                rho_org_matrix = rho_org_matrix+((2.0*rand(size(rho_org_matrix))-1)*noise);
                max(max(rho_org_matrix))
                min(min(rho_org_matrix))

                %%%mirror data over z-axis and r-axis%%%%%%%%%%%%%%
                RR_org=unique([-RR_org RR_org]);
                minRR=min(RR_org);
                maxRR=max(RR_org);
                [R,Z]=meshgrid(RR_org,ZZ_org);
                R_org = R;
                Z_org = Z;
                m=length(RR_org);
                n=length(ZZ_org);
                r_org=reshape(R',n*m,1)';% mirrored r vector
                z_org=reshape(Z',n*m,1)';% original z vector
                size(rho_org_matrix)  %10   149
                rho_org_matrix = [fliplr(rho_org_matrix) rho_org_matrix(:,2:end)];
                rho_org=reshape(rho_org_matrix',n*m,1)';    

                %%%flip about r-axis
                ZZ_org_tmp=[fliplr(-ZZ_org) ZZ_org];
                [R_tmp,Z_tmp]=meshgrid(RR_org,ZZ_org_tmp);
                R_org_tmp = R_tmp;
                Z_org_tmp = Z_tmp;
                m_tmp=length(RR_org);
                n_tmp=length(ZZ_org_tmp);
                r_org_tmp=reshape(R_tmp',n_tmp*m_tmp,1)';% mirrored r vector
                z_org_tmp=reshape(Z_tmp',n_tmp*m_tmp,1)';% original z vector


                %%%%add saturated concentration along droplet boundary
                %r_org_tmp=[r_org_tmp,bdry_x];
                %z_org_tmp=[z_org_tmp,bdry_y];

                size(rho_org_matrix);  %10   297
                rho_org_matrix_tmp = [flipud(rho_org_matrix);rho_org_matrix];
                size(rho_org_matrix_tmp);  %20   297

                rho_org_tmp=reshape(rho_org_matrix_tmp',n_tmp*m_tmp,1)';
                size(rho_org_tmp) ;  %1        5940

                %%%add saturated concentration along the droplet boundary
                %rho_org_tmp = [rho_org_tmp,Cv*ones(size(bdry_x))];

                hRR1_org = hRR1;
                hZZ1_org = hZZ1;
                hZZ2_org = hZZ2;
                hZZ3_org = hZZ3;
                minRR_org = minRR;
                maxRR_org = maxRR;

                min_percent_error_cv_master = [];
                min_rms_error_cv_master = [];
                for smm=1:length(sm)

                    disp('*********Smoothness**********')
                    smooth_gf = sm(smm) 

                    for j=1:h

                        %disp('Normalized grid spacing: ')
                        if j ~= 1
                            hRR1=hRR1/2;% relative cell size in radial direction
                            hZZ1=hZZ1/2;
                            hZZ2=hZZ2/2;
                            hZZ3=hZZ3/2;
                        else
                            hRR1=hRR1_org;
                            hZZ1=hZZ1_org;
                            hZZ2=hZZ2_org;
                            hZZ3=hZZ3_org; 
                            minRR = minRR_org;
                            maxRR = maxRR_org;
                        end
                        disp('Grid size in mm ')
                        hRR1*radius

                        path0=[rootdir,'/Contour_Plots_Scaled_r_',num2str(scaler),'_Scaled_z_',num2str(scalez),'_Noise_',num2str(noise_flag),'_Smooth_',num2str(sm(smm)),'_grid_spacing_',num2str(hRR1*radius),'/'];        
                        path3=[rootdir,'/Evaporation_Rate_Scaled_r_',num2str(scaler),'_Scaled_z_',num2str(scalez),'_Noise_',num2str(noise_flag),'_Smooth_',num2str(sm(smm)),'_grid_spacing_',num2str(hRR1*radius),'/'];

                        if (exist(path0)~=0)
                            rmdir(path0,'s') 
                        end
                        mkdir(path0)

                        if (exist(path3)~=0)
                            rmdir(path3,'s') 
                        end
                        mkdir(path3)       

                        if (smm > 1) || (j > 1)
                            clear RR ZZ rho rho_gradr rho_gradz fit rgrid zgrid dif fitt rhot rgridt zgridt A max_num max_idx row col...
                                fit_gradr fit_gradz;     
                        end

                        %create grid for gridfit
                        RR=minRR-hRR1:hRR1:maxRR+hRR1; % radial vector
                        ZZ=(-0.0*hRR1:hRR1:ZZ4+hRR1); % z vector 

                        minRR = min(RR);
                        maxRR = max(RR);

                        minZZ = min(ZZ);
                        maxZZ = max(ZZ);        

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %analytic solutions for concentration and derivatives
                        [rho,rho_gradr,rho_gradz] = vap_dist_modification_generalization_sub_C50(RR,scaler,ZZ,scalez,C_v,Cf);
                        rho=rho';  %rows: along z, columns: along r
                        rho_gradz=rho_gradz';
                        rho_gradr=rho_gradr';       

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                        %using gridfit: http://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit
                        [rho_gridfit,rgrid,zgrid]=gridfit(r_org,z_org,rho_org,RR,ZZ,'smooth',smooth_gf); %,'tilesize',100);
                        [rho_gridfit_tmp,rgrid_tmp,zgrid_tmp]=gridfit(r_org_tmp,z_org_tmp,rho_org_tmp,RR,unique([fliplr(-ZZ) ZZ]),'smooth',smooth_gf); %,'tilesize',100);
                        rho_gridfit=rho_gridfit_tmp(round(size(rho_gridfit_tmp,1)/2):end,:);
                        zgrid=zgrid_tmp(round(size(zgrid_tmp,1)/2):end,:);

                        %%%set rho_gridfit to zero if negative
                        [row_sz,col_sz] = size(rho_gridfit);
                        %[row_sz,col_sz,neg_rho_gridfit] = find(rho_gridfit<0.0);
                        count = 0;
                        for ii = 1:row_sz
                            for jj = 1:col_sz
                                if rho_gridfit(ii,jj) < 0
                                   count = count + 1;
                                   rho_gridfit(ii,jj) = 0; 
                                end
                            end
                        end
                        %disp('Number of Negative Concentrations from gridfit')
                        %count


                        [fit_gradr,fit_gradz]=gradient(rho_gridfit,hRR1,hRR1);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        [R,Z]=meshgrid(RR,ZZ);

                       %%%%%%%compare gridfit results with theoretical results by visualization
                       %%%regarding concentration and gradient on the uniform grid
                       ControlVolumesLog = ['Producing images for control Volumes (',num2str(i),'/',num2str(length(OptimalSmoothness)),').'];
                        app.CVFluxRMSTextArea.Value = [RMSText';averageRMS;ControlVolumesLog];            
                       plot_conc_grad2(fsize,RR_org,ZZ_org,R_org,Z_org,rho_org_matrix,R,Z,rho,rho_gradr,rho_gradz,...
                                  rgrid,zgrid,rho_gridfit,fit_gradr,fit_gradz,...
                                  smooth_gf,hRR1,path0,rmin/radius,rmax/radius,zmin/radius,zmax/radius,...
                                  bdry_x,bdry_y);
                        ControlVolumesLog = ['Producing images for control Volumes (',num2str(i),'/',num2str(length(OptimalSmoothness)),')..'];
                        app.CVFluxRMSTextArea.Value = [RMSText';averageRMS;ControlVolumesLog];       
                        %pause
                        %%%compare concentration and gradient along r and z slices
%                         plot_conc_grad(fsize,RR_org,ZZ_org,R_org,Z_org,rho_org_matrix,R,Z,rho,rho_gradr,rho_gradz,...
%                                     rgrid,zgrid,rho_gridfit,fit_gradr,fit_gradz,...
%                                     smooth_gf,hRR1,path1,path2) 


                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %define control volume & calculate evaporation rate
                        Mu_array = (rmin:8.0:rmax)/radius;
                        %Mu_array = (7.0:0.5:7.5)/radius;
                        Mu_len = length(Mu_array);
                        %pause

                        height_array = (zmin:2.0:zmax)/radius; %assume that the height of the droplet is less than 1 mm.
                        %height_array = (1.0:0.5:1.5)/radius;
                        Height_len = length(height_array);

                        number_of_control_volume = Height_len*Mu_len;
                        %pause

                        test_array = zeros(number_of_control_volume,5);
                        P = 0;
                        for lc_ind = 1:length(height_array)

                            lc_array = height_array(lc_ind)*ones(size(Mu_array));  

                            trap_approx_int_cv = zeros(size(Mu_array));

                            exact_int_cv = zeros(size(Mu_array));

                            %%%%RMS error will be calculated for each control volume
                            error_rms = zeros(size(Mu_array));
                            avg_per_error = zeros(size(Mu_array));
                            for check=1:length(Mu_array)
                                P = P + 1;
                                Mu = Mu_array(check);
                                lc = lc_array(check);

                                evap_file1 = [path3,'Top_Evaporation_rate_height_',num2str(lc),'_radius_',num2str(Mu)];
                                evap_file2 = [path3,'Side_Evaporation_rate_height_',num2str(lc),'_radius_',num2str(Mu)];

                                %%%top integral
                                %disp('top')
    %                             [trap_approx_int_pointwise,exact_int_pointwise,trap_approx_int,exact_int,tcoords_top,fit_grad_t,exact_grad_t,...
    %                                 trap_approx_int_pointwise_at_0,exact_int_pointwise_at_0,fit_Ct,exact_Ct]
                                [trap_approx_int_pointwise, exact_int_pointwise, trap_approx_int, exact_int, tcoords_top, fit_grad_t, exact_grad_t,...
                                 trap_approx_int_pointwise_at_0, exact_int_pointwise_at_0, fit_C, exact_Ct] =top_integration(hRR1,Mu,lc,R,Z,rho_gridfit,fit_gradr,fit_gradz,...
                                                rho,rho_gradr,rho_gradz,radius,fsize,D);


                                %%%lateral side integral
                                %disp('side')
                                [trap_approx_int_pointwise_cv,exact_int_pointwise_cv,trap_approx_int_cv(check),exact_int_cv(check),tcoords_side,fit_grad_ts,exact_grad_ts,fit_Cs,exact_Cs] =...
                                    side_integration(hRR1,Mu,lc,R,Z,rho_gridfit,fit_gradr,fit_gradz,...
                                                rho,rho_gradr,rho_gradz,radius,fsize,D);


                                %%%%%%%%%%%%%%%%%%%%%
                                trap_approx_int_cv(check) = trap_approx_int_cv(check) + trap_approx_int;

                                exact_int_cv(check) = exact_int_cv(check) + exact_int;            


                               %%%plot control volume
                               ControlVolumesLog = ['Producing images for control Volumes (',num2str(i),'/',num2str(length(OptimalSmoothness)),')...'];
                                app.CVFluxRMSTextArea.Value = [RMSText';averageRMS;ControlVolumesLog]; 
                              plot_convol_Optimal(OptimalSmoothness,sm,P,MinRMSTot, fsize,path0,Mu,lc,rho,rho_gridfit,rgrid,zgrid,tcoords_top,tcoords_side,...
                                       fit_grad_t,fit_grad_ts,exact_grad_t,exact_grad_ts,...
                                       minRR,maxRR,minZZ,maxZZ,radius,bdry_x,bdry_y)
                                 ControlVolumesLog = ['Producing images for control Volumes (',num2str(i),'/',num2str(length(OptimalSmoothness)),')....'];
                                app.CVFluxRMSTextArea.Value = [RMSText';averageRMS;ControlVolumesLog];                                  
                                %pause

                                trap_approx_int_pointwise=trap_approx_int_pointwise(2:end);
                                exact_int_pointwise=exact_int_pointwise(2:end);
                                %abs discretization error
                                top_error_rms = [(trap_approx_int_pointwise_at_0-exact_int_pointwise_at_0)/exact_int_pointwise_at_0;...
                                                 (trap_approx_int_pointwise-exact_int_pointwise)./exact_int_pointwise];
                                indt = find(exact_Ct > -eps);
                                top_error_rms_reduced=top_error_rms(indt); %only count flux components at points with non-negative concentrations

                                side_error_rms = (trap_approx_int_pointwise_cv-exact_int_pointwise_cv)./exact_int_pointwise_cv;
                                inds = find(exact_Cs > -eps);
                                side_error_rms_reduced=side_error_rms(inds); %only count flux components at points with non-negative concentrations

                                cv_error_rms = [top_error_rms_reduced;side_error_rms_reduced];

                                error_rms(check) = 100*sqrt(mean(cv_error_rms.^2));
                                %sqrt(sum(cv_error_rms.^2)/sum([exact_int_pointwise;exact_int_pointwise_cv].^2));

                                avg_per_error(check) = mean(abs(cv_error_rms));
                            end   

                            test = [radius*Mu_array' radius*lc_array' -trap_approx_int_cv' -exact_int_cv'];

                            theory_no_trap = [Mu_array' lc_array' exact_int_cv'];

                            %%%%compute percent error from grid_fit and "gradient"
                            theory = compound_theory;
                            sol = abs(trap_approx_int_cv');
                            percent_error = 100*(sol - theory*ones(size(sol)))./(theory*ones(size(sol)));

                            sol_theo = abs(exact_int_cv');
                            percent_error_theo = 100*(sol - sol_theo)./(sol_theo);

                            test = [test percent_error_theo percent_error avg_per_error' error_rms'];

                            evap_file = [path3,'Evaporation_Rate_Smooth_',num2str(sm(smm)),'_grid_spacing_',num2str(hRR1*radius)];
                            dlmwrite(evap_file,test,'delimiter',' ','precision','%e','-append');
                            test_array((lc_ind-1)*Mu_len+1:lc_ind*Mu_len,1:8) = test;

                        end
                        %add control volume with the smallest percent error to the master list
                        %min_percent_error_cv_master = [min_percent_error_cv_master;sm(smm),hRR1*radius,test_sorted(1,:)];

                        %sort the rms errors in ascending orders
                        [~,idx]=sort(test_array(:,end));

                        test_sorted=test_array(idx,:);
                        
                        
                        for j = 1:length(test_sorted(:,1))
                        TextAreaText1{j} = [num2str(test_sorted(j,1)),'x',num2str(test_sorted(j,2)),', RMS = ',num2str(test_sorted(j,length(test_sorted(1,:)))),'%'];
                        end


                        
                       % averageRMS = ['Average RMS = ',num2str(sum(test_sorted(:,length(test_sorted)))/length(test_sorted(:,length(test_sorted)))),'%']
                        

                        evap_file = [path3,'RMS_Sorted_Evaporation_Rate_Smooth_',num2str(sm(smm)),'_grid_spacing_',num2str(hRR1*radius)];
                        dlmwrite(evap_file,test_sorted,'delimiter',' ','precision','%e');

                        %%add control volume with the smallest rms error to the master list
                        %min_rms_error_cv_master = [min_rms_error_cv_master;sm(smm),hRR1*radius,test_sorted(1,:)];
                    end
                end
            end
        
        
    end
    
    

end