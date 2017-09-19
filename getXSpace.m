function [x_spacing] = getXSpace(SpacingOption)
%Gets the user inputted spacing option and determines the correct spacing

    if SpacingOption == 'Analytical..'
        x_spacing = 0:0.5:50;
    elseif SpacingOption == 'Experimental'
        x_spacing = [0,2,4,6,8,12,16,20,25,30,35];
    end
end

