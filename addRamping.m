function data = addRamping(data,RampToZeroLocation)

        %interpolating IA to x = 50mm
       Y_2 = 0;
       Y_1 = data(length(data),2);
       X_2 = RampToZeroLocation;
       X_1 = data(length(data),1);
       
       m = (Y_2-Y_1)/(X_2-X_1)
       
       b = -m*X_2
       initial_length = length(data); 
       n = 0;
            for k = X_1+1:1:X_2
                n = n + 1;
                Interpolateddata(n,1) = m*k + b
                data(initial_length+n,1) = k;
                data(initial_length+n,2) = Interpolateddata(n,1);
            end
            
            if X_2<50
                
                for i = length(data)+1:2*length(data)
                    X_2 = X_2+1;
                    data(i,1) = X_2;
                    data(i,2) = 0;
                    if(X_2 == 50)
                        break;
                    end
                end
            end

end

