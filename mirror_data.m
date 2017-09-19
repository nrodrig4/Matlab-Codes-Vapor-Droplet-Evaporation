function D = mirror_data (data)
% function D = mirror_data (data) 
%   take 0--n data and flip about X axis (-x..0..x)
% KMN 6/16/2008

% 6/9/2011 - result is always odd, since input is N long.  Output is 2N-1, which is odd.
% 6/16/2009 - odd image size gives symmetric COP.
% 7/7/2008 - copy all cols, not just 2nd one
% 7/1/2008 - added zero to make center of projection an integer
% 6/16/2008 - bug - don't double-count 0 (FBP ignores X).

% The first version adds an extra point at the end, to ensure an odd image
% size - to give a symetric center of projection (COP).
extra_pt=0;
if (extra_pt)
    delta_x = data(end,1)-data(end-1,1);
    extra_x = data(end,1)+delta_x;
    %   regular data |  inverted data     |   extra pt
    %                    (minus 0)       
    x = [data(:,1);     -data(2:end,1);        extra_x        ];
    P = [data(:,2:end);  data(2:end,2:end);    data(end,2:end)];
else
% The second version omits the extra point
    x = [data(:,1);     -data(2:end,1);];      % radial position [mm]
    P = [data(:,2:end);  data(2:end,2:end);];  % rest of data copied
end

% Now, sort by x to put in order
D = sortrows([x P]);

