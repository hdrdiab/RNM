MPnew = MP(1:(pr+2-1)*mr+mr,1:(pr+2-1)*mr+mr);
for i=1:720
    for j = 1:720
        if(i+1 > 720 && j+1 <= 720)
            MPnew(i+1-720,j+1) = MP(i,j);
        elseif(j+1 > 720 && i+1 <= 720)
            MPnew(i+1,j+1-720) = MP(i,j);
        elseif(j+1 > 720 && i+1 > 720)
            MPnew(i+1-720,j+1-720) = MP(i,j);
        else
            MPnew(i+1,j+1) = MP(i,j);
        end
    end
    for j = 1+m:720+m
        if(i+1 > 720 && j+1 <= 720+m)
            MPnew(i+1-720,j+1) = MP(i,j);
        elseif(j+1 > 720+m && i+1 <= 720)
            MPnew(i+1,j+1-720) = MP(i,j);
        elseif(j+1 > 720+m && i+1 > 720)
            MPnew(i+1-720,j+1-720) = MP(i,j);
        else
            MPnew(i+1,j+1) = MP(i,j);
        end
    end
end
% start g with 1
for g = 1:pr
    for i = 1+g*m:720+g*m
        for j = 1+g*m-m:720+g*m-m
            if(i+1 > 720+g*m && j+1 <= 720+g*m-m)
                MPnew(i+1-720,j+1) = MP(i,j);
            elseif(i+1 <= 720+g*m && j+1 > 720+g*m-m)
                MPnew(i+1,j+1-720) = MP(i,j);
            elseif(i+1 > 720+g*m && j+1 > 720+g*m-m)
                MPnew(i+1-720,j+1-720) = MP(i,j);
            else
                MPnew(i+1,j+1) = MP(i,j);
            end
        end
        for j = 1+g*m:720+g*m
            if(i+1 > 720+g*m && j+1 <= 720+g*m)
                MPnew(i+1-720,j+1) = MP(i,j);
            elseif(i+1 <= 720+g*m && j+1 > 720+g*m)
                MPnew(i+1,j+1-720) = MP(i,j);
            elseif(i+1 > 720+g*m && j+1 > 720+g*m)
                MPnew(i+1-720,j+1-720) = MP(i,j);
            else
                MPnew(i+1,j+1) = MP(i,j);
            end
        end
        for j = 1+g*m+m:720+g*m+m
            if(i+1 > 720+g*m && j+1 <= 720+g*m+m)
                MPnew(i+1-720,j+1) = MP(i,j);
            elseif(i+1 <= 720+g*m && j+1 > 720+g*m+m)
                MPnew(i+1,j+1-720) = MP(i,j);
            elseif(i+1 > 720+g*m && j+1 > 720+g*m+m)
                MPnew(i+1-720,j+1-720) = MP(i,j);
            else
                MPnew(i+1,j+1) = MP(i,j);
            end
        end
    end
end
MP(1:(pr+2-1)*mr+mr,1:(pr+2-1)*mr+mr) = MPnew;