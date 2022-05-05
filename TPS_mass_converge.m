rigid = true;
if rigid
    %% convergence with 1st and 4th, 2nd and 3rd materials coupled RIGID
    %thickness = [0.015 0.015 0.015 0.015]; %thickness in m
    %thickness = [0.001 0.02 0.1 0.001]; % starting thicknesseses for rigid
    %fill with large arbitary #s .005 .01 .05 .005
    %thickness = [0.0005 0.0098 0.0459 0.0005]; % rigid optimized %51473,32805
    %thickness = [0.0010 0.0094 0.0437 0.0010];
    thickness = [0.000656 0.010 0.070 0.000656]; %a time  = 500
    [temperatures] = heat_transfer_new(thickness); % gets initial temps
    
    countmax = 4;
    thicknesses = thickness;
    prevthickness = thickness;
    savetemp = temperatures;
    k = [0.9 0.99];
    for j = 1:length(k)
        count = 0;
        while count < countmax && temperatures(4) < 343 %&& mass_sum_rigid < 3.5
            count = count + 1;
            if temperatures(1) < 1923 && temperatures(2) < 573
                thicknesses(1) = thickness(1)*k(j);
                thicknesses(4) = thickness(4)*k(j);
                prevthickness(1) = thickness(1);
                prevthickness(4) = thickness(4);
            end
            if thicknesses(1) == thickness(1) && thicknesses(4) == thickness(4)
                thickness(1) = prevthickness(1);
                thickness(4) = prevthickness(4);
                temperatures(2) = savetemp(2);
                temperatures(4) = savetemp(4);
                break;
            end
            if temperatures(4) > 343
                temperatures(4) = savetemp(4);
                thickness(1) = prevthickness(1);
                thickness(4) = prevthickness(4);
            end
            if temperatures(2) > 573
                temperatures(2) = savetemp(2);
                thickness(1) = prevthickness(1);
                thickness(4) = prevthickness(4);
            end
            savetemp(2) = temperatures(2);
            savetemp(4) = temperatures(4);
            [temperatures] = heat_transfer_new(thicknesses);
            thickness(1) = thicknesses(1);
            thickness(4) = thicknesses(4);
            if temperatures(4) > 343 || temperatures(2) > 573
                temperatures(2) = savetemp(2);
                temperatures(4) = savetemp(4);
                thickness(1) = prevthickness(1);
                thickness(4) = prevthickness(4);
                break;
            end
        end
        
        disp('converged')
        
        count = 0;
        %thicknesses = thickness;
        while count < countmax && temperatures(4) < 343 %&& mass_sum_rigid < 3.5
            count = count + 1;
            if temperatures(1) < 1650 && temperatures(2) < 573
                thicknesses(2) = thickness(2)*k(j);
                prevthickness(2) = thickness(2);
                thicknesses(3) = thickness(3)*k(j);
                prevthickness(3) = thickness(3);
            end
            if thicknesses(2) == thickness(2) && thicknesses(3) == thickness(3)
                thickness(2) = prevthickness(2);
                thickness(3) = prevthickness(3);
                temperatures(2) = savetemp(2);
                temperatures(4) = savetemp(4);
                break;
            end
            %assert(temperatures(2) > 573, 'greater than 573')
            if temperatures(4) > 343
                temperatures(4) = savetemp(4);
                thickness(2) = prevthickness(2);
                thickness(3) = prevthickness(3);
            end
            if temperatures(2) > 573
                temperatures(2) = savetemp(2);
                thickness(2) = prevthickness(2);
                thickness(3) = prevthickness(3);
            end
            savetemp(2) = temperatures(2);
            savetemp(4) = temperatures(4);
            [temperatures] = heat_transfer_new(thicknesses);
            thickness(2) = thicknesses(2);
            thickness(3) = thicknesses(3);
            if temperatures(4) > 343 || temperatures(2) > 573
                temperatures(2) = savetemp(2);
                temperatures(4) = savetemp(4);
                thickness(2) = prevthickness(2);
                thickness(3) = prevthickness(3);
                break;
            end
        end
    
        disp('converged')
    end
    
    tps_type = 'R';
    %prevthickeness if running script, otherwise just insert optimized
    %thickness here 
    [mass, total_mass] = get_TPS_mass(prevthickness, tps_type);
%% further convergenence, decoupling/increasing thicknesses to optimize

%while temperatures(2) < 573 && temperatures(1) < 1650 && temperatures(4) < 343
       
%end

else
%% convergence of FLEXIBLE
    %thickness = [0.00035 0.01 0.001 0.00025]; % starting thicknesses for %flexible mass 5.4550kg, but doesnt work 
    %thickness = [0.00025 0.0092 0.028 0.00032]; % optimized thickness %0.00031973
    %thickness = [0.00025 0.006 0.0219 0.0025]; 9.1316kg (1.38662194226018	1.12860707056498	5.15460097091327	1.46179893281638)
    %thickness = [0.00033 0.007 0.02 0.0005]; % start
    %thickness = [0.00027 0.015 0.015 0.00025]; % works for flexible
    %thickness = [0.0003 0.0054 0.015 0.00025];
    %thickness = [0.0025 0.0025 0.0025 0.0025];
    thickness = [0.00027 0.006 0.02 0.00025];

    [temperatures] = heat_transfer_new(thickness); % gets initial temps
    countmax = 10;
    thicknesses = thickness;
    prevthickness = thickness;
    savetemp = temperatures;
    k = [0.95 0.99];
    for j = 1:length(k)
    count = 0;
    while count < countmax && temperatures(4) < 343 %&& mass_sum_flexible < 5.5
        count = count + 1;
        if temperatures(1) < 923 && thickness(1) > 0.00025
            thicknesses(1) = thickness(1)*k(j);
            prevthickness(1) = thickness(1);
        end
        if thicknesses(1) == thickness(1)
            thickness(1) = prevthickness(1);
            temperatures(1) = savetemp(1);
            temperatures(2) = savetemp(2);
            temperatures(4) = savetemp(4);
            break;
        end
        if temperatures(4) > 343
            temperatures(4) = savetemp(4);
            thickness(1) = prevthickness(1);
        end
        if temperatures(1) > 923
            temperatures(1) = savetemp(1);
            thickness(1) = prevthickness(1);
        end
        savetemp(1) = temperatures(1);
        savetemp(2) = temperatures(2);
        savetemp(4) = temperatures(4);
        if thicknesses(1) < 0.00025
           thicknesses(1) = prevthickness(1);
            break;
        end
        [temperatures] = heat_transfer_new(thicknesses);
        thickness(1) = thicknesses(1);
        if temperatures(4) > 343 || temperatures(1) > 923
            temperatures(1) = savetemp(1);
            temperatures(2) = savetemp(2);
            temperatures(4) = savetemp(4);
            thickness(1) = prevthickness(1);
            break;
        end
    end
    
    disp('converged')
    
 count = 0;
    while count < countmax && temperatures(4) < 343 %&& mass_sum_flexible < 5.5
        count = count + 1;
        if temperatures(1) < 923 && thickness(4) > 0.00025
            thicknesses(4) = thickness(4)*k(j);
            prevthickness(4) = thickness(4);
        end
        if thicknesses(4) == thickness(4)
            thickness(4) = prevthickness(4);
            temperatures(1) = savetemp(1);
            temperatures(2) = savetemp(2);
            temperatures(4) = savetemp(4);
            break;
        end
        if temperatures(4) > 343
            temperatures(4) = savetemp(4);
            thickness(4) = prevthickness(4);
        end
        if temperatures(1) > 923
            temperatures(1) = savetemp(1);
            thickness(4) = prevthickness(4);
        end
        savetemp(1) = temperatures(1);
        savetemp(2) = temperatures(2);
        savetemp(4) = temperatures(4);
        if thicknesses(4) < 0.00025
            thicknesses(4) = prevthickness(4);
            break;
        end
        [temperatures] = heat_transfer_new(thicknesses);
        thickness(4) = thicknesses(4);
        if temperatures(4) > 343 || temperatures(1) > 923
            temperatures(1) = savetemp(1);
            temperatures(2) = savetemp(2);
            temperatures(4) = savetemp(4);
            thickness(4) = prevthickness(4);
            break;
        end
    end
    
    disp('converged')

    count = 0;
    while count < countmax && temperatures(4) < 343 %&& mass_sum_flexible < 5.5
        count = count + 1;
        if temperatures(1) < 923 && thickness(2) > 0.0050
            thicknesses(2) = thickness(2)*k(j);
            prevthickness(2) = thickness(2);
        end
        if thicknesses(2) == thickness(2)
            thickness(2) = prevthickness(2);
            temperatures(1) = savetemp(1);
            temperatures(2) = savetemp(2);
            temperatures(4) = savetemp(4);
            break;
        end
        if temperatures(4) > 343
            temperatures(4) = savetemp(4);
            thickness(2) = prevthickness(2);
        end
        if temperatures(1) > 923
            temperatures(1) = savetemp(1);
            thickness(2) = prevthickness(2);
        end
        if thicknesses(2) < 0.0050
            break;
        end
        savetemp(1) = temperatures(1);
        savetemp(2) = temperatures(2);
        savetemp(4) = temperatures(4);
        [temperatures] = heat_transfer_new(thicknesses);
        thickness(2) = thicknesses(2);
        if temperatures(4) > 343 || temperatures(1) > 923
            temperatures(2) = savetemp(2);
            temperatures(4) = savetemp(4);
            thickness(2) = prevthickness(2);
            break;
        end
    end
    
    disp('converged')
    
    count = 0;
    while count < countmax && temperatures(4) < 343 %&& mass_sum_flexible < 5.5
        count = count + 1;
        if temperatures(1) < 923 && thickness(3) > 0.00015
            thicknesses(3) = thickness(3)*k(j);
            prevthickness(3) = thickness(3);
        end
        if thicknesses(3) == thickness(3)
            thickness(3) = prevthickness(3);
            temperatures(1) = savetemp(1);
            temperatures(2) = savetemp(2);
            temperatures(4) = savetemp(4);
            break;
        end
        if temperatures(4) > 343
            temperatures(4) = savetemp(4);
            thickness(3) = prevthickness(3);
        end
        if temperatures(1) > 923
            temperatures(1) = savetemp(1);
            thickness(3) = prevthickness(3);
        end
        if thicknesses(3) < 0.00015
            break;
        end
        savetemp(1) = temperatures(1);
        savetemp(2) = temperatures(2);
        savetemp(4) = temperatures(4);
        [temperatures] = heat_transfer_new(thicknesses);
        thickness(3) = thicknesses(3);
        if temperatures(4) > 343 || temperatures(1) > 923
            temperatures(2) = savetemp(2);
            temperatures(4) = savetemp(4);
            thickness(3) = prevthickness(3);
            break;
        end
    end
    
    disp('converged')
    
    end
    
    tps_type = 'F';
    [mass, masses] = get_TPS_mass(prevthickness, tps_type);
end
