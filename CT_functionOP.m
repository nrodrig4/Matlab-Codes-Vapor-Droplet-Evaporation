function [results]=CT_functionOP(CT,filedir,detail,power,z,mu)
 figure('Visible','off');
% ct_data --- executive script for FilmLab Computerized Tomography
%{
UPDATE NOTICE: MAJOR UPDATE (June 1, 2015)
    Changed how mu was found by adding a GUI.
    Eliminated the process of using an excell spread sheet, instead pulls
    data from text document that is produced by
    Variable_Polynomial_Fitting. Utalizing new function, import_txt,
    instead of import_data.
    Changed how files are saved so user imputs what z level is being
    evaluated.
    Formatted it so that it is a function refrenced by the script
    Variable_Polynomial_Fitting_and_CT.
%}
%{
% Kevin Nickels 6/10/2008
% Cleanups to reflect current lab practices KMN 6/9/2011
%    6/12/12 Comments and review of above... kmn
%    * adjusted normalization of concentration (mu)
%    * now storing mu in input spreadsheet
%    * matlab now processes one PIC at a time.  You can save/rename/restore
%      results.mat to compare across datasets.
%    6/13/12 kmn
%    * fixed 1+ bug in series defn (match comments below)
%    * added explicit check to see if mu is defined in input spreadsheet,
%      issue warning (should be error?) if mu not defined.
%}

    %savepath = filedir;
    
    if length(power)~=1 %This if esle statement prompts the user to select which polynomial(s) they want analyzed. If there was only one polynomial run then it is assumed the user wants to run that one through the CT.
        expo=input('Which power polynomial(s) would you like to evaluate? (Enter as a vector. EX: [5 7 8]): ');
    else
        expo=power;
    end

for trials=expo
series=trials;    

    D=mirror_data(CT);
    x=D(:,1);
    
    %% Data is in, now proceed to format
    scaling = x(2)-x(1); % in mm/pixel
    % scaling = 1;    
          
% --------------------------------------------------------------------
    %% Pad/reformat for FBP
%{   
(from iradon() -  The projections are zero-padded to a power of 2
before filtering to prevent spatial domain aliasing and to speed up
the FFT.)
    
Here, I was playing with adding some more padding as a way of
filtering in the spatial domain.  kmn 7/2009
Around 100 points seems to work well for our data 6/2011
%}
    
    pad=100;
    PIC = D(:,series+1);
    PIC_pad = [ zeros(pad,1);          PIC;   zeros(pad,1)];
    %x_pad = [ min(x)-1:pad*scaling     x'    max(x)+1:pad*scaling]';
    
% --------------------------------------------------------------------
    %% make a circularly symmetric sinugraph
    Num = size(PIC_pad,1);
    Sinugraph = PIC_pad*ones(size(1:Num));% R is the full set of projections.  Here, we assume they are all the same.
  
    
    theta = 0:180/(Num-1):180; % theta is the view angles.  We say regular spacing from 0 to 180.
    
% --------------------------------------------------------------------
    %% Now, take inverse radon transform of Projections
    % This is where the MATLAB heavy-lifting is.  All the rest of the stuff is
    %  so that this one line can happen!
    Ir = iradon(Sinugraph,theta,'linear','Ram-Lak');
    %{
    % We played with different filters (Su'10), and concluded that the
    % RamLak is the proper filter to use in most cases.
    % 6/9/2011
    %    Ir = iradon(R,theta,'linear','Shepp-Logan');
    %    Ir = iradon(R,theta,'linear','Cosine');
    %    Ir = iradon(R,theta,'spline','Cosine');
    %    Ir = iradon(R,theta,'pchip','Cosine');
    %    Ir = iradon(R,theta,'linear','Hamming');
    %    Ir = iradon(R,theta,'linear','Hann');
    %    Ir = iradon(R,theta,'linear','None');
    %}
% --------------------------------------------------------------------
    %% Check - take radon transform of reconstructed Images
    % One way to verify a CT is to take the radon transform (projection) of
    % the reconstruction, then compare the original and reprojected PICs.
    % Ideally, P2==PIC
    %[PIC2,x2] = radon(Ir,theta);
    
    % This computes the spatial extent of the reconstruction.
    c = floor((size(Ir)+1)/2);
    Ir_x = (scaling:scaling:scaling*size(Ir))-scaling*c(1);

    
	abs_conc = Ir(:,c(2))/scaling; % remap from au/cell to au/mm
	conc = abs_conc / mu;


    % --------------------------------------------------------------------
    %% Display Various combinations and purmutations of result plots.
    
    %% Show reconstructed images
%     showIr=0;
%     if(showIr)
%         figure(1);
%         imshow(Ir,[]);
%     end
    
     
    %% show a slice through the reconstruction
    showIrslice=1;
    if(showIrslice)
      

        
        
    end
    
    %% show the original and reprojected PICs on same axes
    % PIC and P2 for data
    showdata=1;
    if(showdata)
        %figure(2);
        
        iptsetpref('ImshowAxesVisible','on');
        
       % for looking throught PICS in rapid code iteration mode 
       % index=315;
       % plot(x,R(:,index),'x',scaling*x2,P2(:,index),'o');
       
       
        % manually check the 3D reconstructed PICs (P2)
        %mesh(PIC2); title('Reconstructed PICs');
        
        % maybe plot errorbars plus mean? - kmn 6/11/12
        %saveas(gcf,[savepath,'/Reconstructed PICs'],'jpg'); %this saves the graphs
        %saveas(gcf,[savepath,'/Reconstructed PICs'],'fig');
        % for seeing a slice of the PICs
       
		  % add mu plus units to this label later - pkz 6/11/12
    end

%% Save results for later analysis   
results = [Ir_x' conc];

%save([filedir,'/resultsZ=',num2str(z),' ',detail,' ',num2str(series),'.txt'],'results','-ascii')
end
disp('the results of the CT are saved as a text file titled results Z= #')
end

