function [resultsOutput,new_results_Size] = IA_and_CT(app,z_length,z_loc,mu,compound,dirext,M,C_v,R,x_spacing,noise_IA,rootdir,SpacingOption,OptimalPower,power,Fitting,noise_flag)
        

        if exist( fullfile([rootdir,'/RMS.txt']), 'file' )
            delete([rootdir,'/RMS.txt'])
        end
    RMSresults = [rootdir,'/RMS.txt'];
    app.FittingCTRMSTextArea.Value = {''};

    error_rel_rms = zeros(z_length,1); %include all points from 0 to 50 (without modifying the analytical sol for comparison)
    error_rel2_rms = zeros(z_length,1); %include all points from 0 to the place where analytical sol becomes 0 or negative
    
    if OptimalPower == 1 && isequal(Fitting,'Polyfit')   
        powers = 6:1:20;
        [resultsOutput,new_results_Size] = optimizedPolyFit(app,z_length,z_loc,mu,compound,M,C_v,R,x_spacing,noise_IA,rootdir,SpacingOption,powers,Fitting,error_rel_rms,error_rel2_rms,RMSresults)
    else
       [resultsOutput,new_results_Size] = normalIA_and_CT(app,z_length,z_loc,mu,compound,M,C_v,R,x_spacing,noise_IA,rootdir,SpacingOption,power,Fitting,error_rel_rms,error_rel2_rms,RMSresults);
    end


end

