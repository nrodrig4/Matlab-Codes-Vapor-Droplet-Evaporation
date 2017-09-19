function [resultsOutput,new_results_Size] = optimizedPolyFit(app,z_length,z_loc,mu,compound,M,C_v,R,x_spacing,noise_IA,rootdir,SpacingOption,powers,Fitting,error_rel_rms,error_rel2_rms,RMSresults)
    for i = 1:z_length
        for j = 1:length(powers)
            
            r_max = 50;
            z_o = z_loc(i);
            power = powers(j);
           
            SavesFolder = [rootdir,'/z_',num2str(z_o)];
            mkdir(SavesFolder)

            RMSdetails = [SavesFolder,'/RMS_details.txt'];
            RMS_mod_details = [SavesFolder,'/RMS_mod_details.txt'];
            neg_analytical_C = [SavesFolder,'/Neg_analytical_vals.txt'];
            
            data = IA_minus_Concentration_mod(noise_IA,z_o,mu,M,C_v,R,SavesFolder,x_spacing);
           

            if SpacingOption ~= 'Analytical..' 
                
                C_f = 0.0166%2*C_v/pi*atan(R/(sqrt(0.5*(rw^2+z_oo^2-R^2+ ...
                %5sqrt((rw^2+z_oo^2-R^2)^2+4*z_oo^2*R^2))))); %analytical C at (rw,z_oo)
                r = 50;

                x_loc = 0:0.5:50;

                for n = 1:length(x_loc)

                    x_o = x_loc(n);
                    a = sqrt(r^2-x_o^2);
                    fun = @(t) compute_C(t,C_v,x_o,z_o,R,C_f);
%                     disp('int1')
                    solution = integral(fun,-a,a); %line integral along the y-direction
%                     disp('int2')
                    %unit?
                    if(solution==0)
                         RampToZeroLocation = x_o;
                        break;
                    end
                   %solution = solution*1000/(1*10^6*M)*mu-(2.0*rand(size(solution))-1*ones(size(solution)))*noise_IA; %unit?


                end 
                
               data = addRamping(data,RampToZeroLocation);
            end 

  
            
            results=Variable_Fitting_CT_NR_OP(data,z_o,compound,power,SavesFolder,mu,Fitting);

            results_tmp = find(results(:,1)<35 + 1e-5 & results(:,1)>-eps);

            new_results_Size = [results_tmp(1) results_tmp(end)]; %still have all the points from CT analysis

            resultsOutput(:,i) = results(:,2)/(1/(1*10^3*M)); 

          
            
            
            t1=results(new_results_Size(1):new_results_Size(2),1);


            rw = r_max;  %sqrt(50^2 - z_o^2);  %%%assume that C(r,50) is also zero (= C(50,z)) 
            z_oo = 0;
            C_f = 2*C_v/pi*atan(R/(sqrt(0.5*(rw^2+z_oo^2-R^2+ ...
            sqrt((rw^2+z_oo^2-R^2)^2+4*z_oo^2*R^2))))); %analytical C at (rw,z_oo)

            temp = zeros(size(t1));
            rel_rms = zeros(size(temp));
            rel2_rms = [];
            count = 0;
            neg_ana_C = [];
            for k = 1:length(t1)
                x_o=t1(k);
                temp(k) = (2*C_v./pi*atan(R./(sqrt(0.5*(x_o.^2+z_o.^2-R.^2 + sqrt((x_o.^2+z_o.^2 ...
                                              -R.^2).^2 + 4*z_o.^2*R^2)))))-C_f)*1/(1*10^3*M);%units mole/cm^3
                rel_rms(k) = (results(new_results_Size(1)+k-1,2)-temp(k))/temp(k);

                if (temp(k) > 0.0) && (x_o < r_max + 1e-5)
                    rel2_rms = [rel2_rms;rel_rms(k)];
                end

            end
            
            error_tmp = [t1 results(new_results_Size(1):new_results_Size(2),2) temp rel_rms];

            %%%compute RMS error
            error_rel_rms_power(j) = 100*sqrt(mean(rel_rms.^2));
            error_rel2_rms_power(j) = 100*sqrt(mean(rel2_rms.^2))
      
            %error_tmp = [z_o count error_rel_rms(i) error_rel2_rms(i)];
            %save(RMSresults,'error_tmp','-ascii','-append') 

            
            TextAreaText{j} = ['z = ',num2str(z_o),', Power = ',num2str(power),', RMS = ',num2str(error_rel2_rms_power(j)),'%']
            
            [LowestRMS2,powerIndex2] = min(error_rel2_rms_power);  
            %[LowestRMS1,powerIndex1] = min(error_rel_rms_power);
            error_rel2_rms(i) = LowestRMS2
      
            %error_rel_rms(i) = LowestRMS1;
            %error_tmp = [z_o count error_rel_rms(i) error_rel2_rms(i)];
            optimalPower(i) = powers(powerIndex2);
            %error_rel_rms_power(i) = LowestRMS1;
            OP = optimalPower(i);
            optimalPowerForZ{i} = ['Power ',num2str(OP),' for z = ',num2str(z_o)]
            app.FittingCTRMSTextArea.Value = [TextAreaText';optimalPowerForZ{i}];     
            
            clc;
            close all;
            clearvars -except j i select_PFitcase new_results_Size p_val sm_val z_loc ...
            choice p_length z_length power M C_v R mu dirext rootdir z_o gs TextAreaText ...
            SavesFolder resultsOutput results C_f compound s_choice sm_length noise_flag powers....
            x_spacing resultPos data SpacingOption Spacing GridFit noise_IA Fitting app RMSresults...
            optimalPower optimalPowerForZ error_rel_rms error_tmp error_rel2_rms error_rel2_rms_power;       
           
        end
        clear error_rel2_rms_power;
        TextAreaText = {};
      optimalPowerForZ{i}
       averageRMS = ['Average RMS = ',num2str(sum(error_rel2_rms)/length(error_rel2_rms))]
      
       app.FittingCTRMSTextArea.Value = [optimalPowerForZ';averageRMS];
    end 
      
    
    
    
       for i = 1:z_length
            r_max = 50;
            z_o = z_loc(i);
            power = optimalPower(i);
            
            
            SavesFolder = [rootdir,'/z_',num2str(z_o)];
            mkdir(SavesFolder)

            RMSdetails = [SavesFolder,'/RMS_details.txt'];
            RMS_mod_details = [SavesFolder,'/RMS_mod_details.txt'];
            neg_analytical_C = [SavesFolder,'/Neg_analytical_vals.txt'];


            data = IA_minus_Concentration_mod(noise_IA,z_o,mu,M,C_v,R,SavesFolder,x_spacing);




            ax = app.UIAxes;
            app.UIAxes_2.NextPlot = 'replacechildren';
            fig=figure;           
            title(app.UIAxes, ['IA output: z = ',num2str(z_o)])
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'IA')

            scatter(ax,data(:,1),data(:,2));
            saveIAScatterFigure(data,SavesFolder,z_o)


         
            if SpacingOption ~= 'Analytical..' 
                
                C_f = 0.0166%2*C_v/pi*atan(R/(sqrt(0.5*(rw^2+z_oo^2-R^2+ ...
                %5sqrt((rw^2+z_oo^2-R^2)^2+4*z_oo^2*R^2))))); %analytical C at (rw,z_oo)
                r = 50;

                x_loc = 0:0.5:50;

                for n = 1:length(x_loc)

                    x_o = x_loc(n);
                    a = sqrt(r^2-x_o^2);
                    fun = @(t) compute_C(t,C_v,x_o,z_o,R,C_f);
                    %disp('int1')
                    solution = integral(fun,-a,a); %line integral along the y-direction
                    %disp('int2')
                    %unit?
                    if(solution==0)
                         RampToZeroLocation = x_o;
                        break;
                    end
                   %solution = solution*1000/(1*10^6*M)*mu-(2.0*rand(size(solution))-1*ones(size(solution)))*noise_IA; %unit?


                end 
                
               data = addRamping(data,RampToZeroLocation);
            end 
    
 

            ax = app.UIAxes;
            fig=figure;           
            title(app.UIAxes, ['IA output: z = ',num2str(z_o)])
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'IA')
            scatter(ax,data(:,1),data(:,2));

            saveIAScatterFigure(data,SavesFolder,z_o);

            % Starting PolyFit and CT
            results=Variable_Fitting_CT_NR(app,data,z_o,compound,power,SavesFolder,mu,Fitting);

            results_tmp = find(results(:,1)<35 + 1e-5 & results(:,1)>-eps);

            new_results_Size = [results_tmp(1) results_tmp(end)]; %still have all the points from CT analysis

            resultsOutput(:,i) = results(:,2)/(1/(1*10^3*M)); 

            fresults = [SavesFolder,'/CT_results_z_',num2str(z_o),'.mat'];
            save(fresults,'resultsOutput') 

            t1=results(new_results_Size(1):new_results_Size(2),1);


            rw = r_max;  %sqrt(50^2 - z_o^2);  %%%assume that C(r,50) is also zero (= C(50,z)) 
            z_oo = 0;
            C_f = 2*C_v/pi*atan(R/(sqrt(0.5*(rw^2+z_oo^2-R^2+ ...
            sqrt((rw^2+z_oo^2-R^2)^2+4*z_oo^2*R^2))))); %analytical C at (rw,z_oo)

            temp = zeros(size(t1));
            rel_rms = zeros(size(temp));
            rel2_rms = [];
            count = 0;
            neg_ana_C = [];
            for j = 1:length(t1)
                x_o=t1(j);
                temp(j) = (2*C_v./pi*atan(R./(sqrt(0.5*(x_o.^2+z_o.^2-R.^2 + sqrt((x_o.^2+z_o.^2 ...
                                              -R.^2).^2 + 4*z_o.^2*R^2)))))-C_f)*1/(1*10^3*M);%units mole/cm^3
                rel_rms(j) = (results(new_results_Size(1)+j-1,2)-temp(j))/temp(j);

                if (temp(j) > 0.0) && (x_o < r_max + 1e-5)
                    rel2_rms = [rel2_rms;rel_rms(j)];
                end

            end

    %         results(:,2)
    %         pause




            error_tmp = [t1 results(new_results_Size(1):new_results_Size(2),2) temp rel_rms];
            save(RMSdetails,'error_tmp','-ascii')  
            save(RMS_mod_details,'rel2_rms','-ascii') 
            if (count>0)
                save(neg_analytical_C,'neg_ana_C','-ascii')
            end

            %%%compute RMS error
            error_rel_rms(i) = 100*sqrt(mean(rel_rms.^2));
            error_rel2_rms(i) = 100*sqrt(mean(rel2_rms.^2));
            error_tmp = [z_o count error_rel_rms(i) error_rel2_rms(i)];
            save(RMSresults,'error_tmp','-ascii','-append') 

            errorCalc(i) = error_rel2_rms(i);
            TextAreaText{i} = ['z = ',num2str(z_o),', RMS = ',num2str(error_rel2_rms(i)),'%'];
            averageRMS = ['Average RMS = ',num2str(sum(errorCalc)/length(errorCalc)),'%'];
            TextAreaText
            Location = ['Location ',num2str(i),' of ',num2str(z_length)];
            averageRMS = ['Average RMS = ',num2str(sum(errorCalc)/length(errorCalc))];
            
            app.FittingCTRMSTextArea.Value = [optimalPowerForZ';Location;averageRMS];

            ax = app.UIAxes_2;
            fig=figure; 
            scatter(ax,results(:,1),results(:,2),'o');
            ax.NextPlot = 'add';
            scatter(ax,t1,temp,'x')
            if (count>0)
                scatter(ax,neg_ana_C(:,1),neg_ana_C(:,2),'+','MarkerEdgeColor',[0 0 0],'LineWidth',1)
            end

            xlabel(ax,'r axis [mm]');
            ylabel(ax,'C(r,z) [mole/cm^3]');
            if (count>0)
                legend(ax,'Fitted Solution', 'Exact Solution','Org. Neg. Analyt. Sol')
            else
                legend(ax,'Fitted Solution', 'Exact Solution')
            end
            title(ax,{['z = ',num2str(z_o),', Rel RMS (org.) = ',num2str(error_rel_rms(i)),'%, Rel RMS (mod.) = ',num2str(error_rel2_rms(i)),'%'];...
                        [' # points with neg analytical conc. = ',num2str(count)]});
            xlim(ax,[0 r_max]);
            saveCTComparison(results,temp,z_o,t1,SavesFolder,neg_ana_C,error_rel_rms(i),count,error_rel2_rms(i),r_max);
            clc;
            close all;
            clearvars -except j i select_PFitcase new_results_Size p_val sm_val z_loc optimalPower...
            choice p_length z_length power M C_v R mu dirext rootdir z_o gs TextAreaText...
            SavesFolder resultsOutput results C_f compound s_choice sm_length noise_flag powers optimalPowerForZ....
            x_spacing resultPos data SpacingOption Spacing GridFit noise_IA Fitting app RMSresults errorCalc;       

       end
       
    

end

function fun = compute_C(t,C_v,x_o,z_o,R,C_f)
               fun = 2.*C_v./pi.*atan(R./(sqrt(0.5.*(x_o.^2+t.^2+z_o.^2-R.^2+ ...
                    sqrt((x_o.^2+t.^2+z_o.^2-R.^2).^2+4.*z_o.^2.*R.^2)))))-C_f; %t represents y coordinates

   if (fun < 0.0)
       fun = 0.0*fun;
   end
end

