function [resultpos] = Format_Part1_to_Part2_mod(rootdir,resultsOutput,select_case,ind)
    
    resultpos = resultsOutput(ind:end,:);

    if select_case == 1
        save([rootdir,'/resultpos_without_polyfit.mat'],'resultpos')
    else
        save([rootdir,'/resultpos_with_polyfit.mat'],'resultpos')
    end
end

