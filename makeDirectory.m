function rootdir = makeDirectory(SpacingOption,OptimalPower, power,Fitting,noise)
%Makes root directory to store images and data in

    if isequal(Fitting,'Polyfit')
        if OptimalPower == 1
            rootdir = [noise,'/',SpacingOption,'_',Fitting,'_Optimal'];
        else
            rootdir = [noise,'/',SpacingOption,'_',Fitting,'_Power_',num2str(power)];
        end 
    else
        rootdir = [noise,'/',SpacingOption,'_',Fitting];
    end
    if (exist(rootdir)~=1)
        mkdir(rootdir)
    end
        

end

