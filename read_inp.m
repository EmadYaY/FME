%--------------------------------------------------------------------------
% Developed by Dr. Hojjat Badnava  
% contact: badnava@bkatu.ac.ir or ho.badnava@gmail.com
%--------------------------------------------------------------------------


function [GEO] = read_inp(file_name)



CON.OutName                 = file_name(1:length(file_name)-4);

[n_particles, n_element]    = read_n_particles(file_name);

GEO.NP                      = n_particles;
GEO.NE                      = n_element;

% Initial coordinates
GEO.XP          = zeros(2,GEO.NP);

% Element connectivity
GEO.CONN        = zeros(4,GEO.NE);

in_step         = 0;

n_BC            = 0;
n_IC            = 0;
n_MAT           = 0;
n_STEP          = 0;
n_Section       = 0;
n_ALEC          = 0;

% -------------------------------------------------------------------------
%  Initializing defaault values

cfl                     = 0.3;

output_frequency        = 10;

CON.Time_Integrator     = 1;

CON.Kernel_Type         = 1;

DAT.Upwind_p            = 1.0;

DAT.Upwind_J            = 0.0;

DAT.Alfa                = 1.0;

CON.Algorithm           = 1;

CON.Kernel_Type_Plot    = 1;

CON.Restart_Run         = 1;

% -------------------------------------------------------------------------
fid =   fopen(file_name);
line = upper (fgets(fid));

while ~feof(fid)
     %# read line by line
%     x=sscanf(fgetl(fid),'%f');

%     disp(x.')
  
    
    % skip comment lines
    if strcmp( line(1:2),'**' )
%          disp(line);
         line =  upper (fgets(fid));
%     elseif length (line) < 5
%         line =  upper (fgets(fid));
    elseif length(line) > 4 && strcmp( line(1:5),'*NODE' )

        
        while(true)
            
            line = upper (fgets(fid));
            
            if strcmp( line(1:1),'*' )
                break;
            end
            
            data = strsplit(line,',');
            
            node_lable = str2double (data(1));
            
            node_coords(1) = str2double (data(2));
            node_coords(2) = str2double (data(3));
%             node_coords(3) = str2double (data(4));
            

            GEO.XP(:,node_lable) = node_coords(:);
            
                
        end
    elseif length(line) > 7 && strcmp( line(1:8),'*ELEMENT' )
        
        line = strrep(line,' ','');
        line = strrep(line,'\n','');
        line = regexprep(line,'[\n\r]+','');
        data = strsplit(line,',');
        
        temp2 = data{2};
        
        type = strrep(temp2(6:length(temp2)),'-','_');        

        
        if contains(type , 'C3D8R')
            GEO.ElementType = 'C3D8R';         
        elseif contains(type , 'C3D8T')
            GEO.ElementType = 'C3D8T'; 
        elseif contains(type , 'C3D8')
            GEO.ElementType = 'C3D8'; 
        elseif contains(type , 'CPS4')
            GEO.ElementType = 'CPS4';             
        else
            error('Invalid Element type');
        end 
        
        while(true)
            
            line = upper (fgets(fid));
            
            if strcmp( line(1:1),'*' )
                break;
            end
            
            data = strsplit(line,',');
            
            el_lable = str2double (data(1));
            
            connectivity(1) = str2double (data(2));
            connectivity(2) = str2double (data(3));
            connectivity(3) = str2double (data(4));
            connectivity(4) = str2double (data(5));
%             connectivity(5) = str2double (data(6));
%             connectivity(6) = str2double (data(7));
%             connectivity(7) = str2double (data(8));
%             connectivity(8) = str2double (data(9));
            

            GEO.CONN(:,el_lable) = connectivity(:);
            
                
        end  
    elseif length(line) > 4 && strcmp( line(1:5),'*NSET' )

        line = strrep(line,' ','');
        line = strrep(line,'\n','');
        line = regexprep(line,'[\n\r]+','');
        data = strsplit(line,',');
        
        temp = data{2};
        
        set_name = temp(6:length(temp));
        set_name = strrep(set_name,'-','_');
        
        if length(data) > 2
            
            if contains(line, 'GENERATE')
             
                line = upper (fgets(fid));
                
                data = strsplit(line,',');
                
                start = str2double (data(1));
                end_ = str2double (data(2));
                step_ = str2double (data(3));
                
                n_value = (end_ - start)/step_ + 1;
                
                set_value = zeros(n_value,1);

                for i=1:n_value
                    set_value(i) = start + step_ * (i-1);
                end
                                
                GEO.NSET.(set_name) = set_value;
                
                line = upper (fgets(fid));
                
            else
                
                set_value = zeros(1,1);
                
                ns_value = 0;
                
                while(true)
                    
                    line = upper (fgets(fid));
                    
                    if strcmp( line(1:1),'*' )
                        GEO.NSET.(set_name) = set_value;
                        break;
                    end
                    
                    data = strsplit(line,',');
                    
                    for i=1:length(data)
                        ns_value = ns_value + 1;
                        set_value(ns_value) = str2double (data(i));
                    end
                    
                end
            end
        else
            
            set_value = zeros(1,1);

            ns_value = 0;
            
            while(true)
                
                line = upper (fgets(fid));
                
                if strcmp( line(1:1),'*' )
                    GEO.NSET.(set_name) = set_value;
                    break;
                end
                
                data = strsplit(line,',');
                
                for i=1:length(data)
                    ns_value = ns_value + 1;
                    set_value(ns_value) = str2double (data(i));
                end                
                
            end
        end
        
        
    elseif length(line) > 8 && strcmp( line(1:9),'*BOUNDARY' )
        
        
        while(true)
            
            line = upper (fgets(fid));
            line = strrep(line,' ','');
            line = strrep(line,'\n','');
            line = regexprep(line,'[\n\r]+','');
            
            if strcmp( line(1:1),'*' )
                break;
            end
            
            n_BC = n_BC + 1 ;
            
            data = strsplit(line,',');
            
            temp = data{1};
            set_name = strrep(temp,'-','_');
            
            GEO.BC(n_BC).set = set_name;
            GEO.BC(n_BC).dof = str2double (data(2));
            
            if(in_step == 0)
                GEO.BC(n_BC).value = 0.0;
            else
                GEO.BC(n_BC).value = str2double (data(4));
            end
            
        end
    elseif length(line) > 18 && strcmp( line(1:19),'*INITIAL CONDITIONS' )
                        
        line = strrep(line,' ','');
        line = strrep(line,'\n','');
        line = regexprep(line,'[\n\r]+','');
        data = strsplit(line,',');
        
        temp2 = data{2};
        
        type = strrep(temp2(6:length(temp2)),'-','_');        

        
        if contains(type , 'VELOCITY')
            type = 'VELOCITY';
        elseif contains(type , 'TEMPERATURE')    
            type ='TEMPERATURE';
        else
            error('Invalid INITIAL CONDITION type');
        end    
        
        
        if strcmp(type , 'VELOCITY')
            while(true)
                
                line = upper (fgets(fid));
                line = strrep(line,' ','');
                line = strrep(line,'\n','');
                line            = regexprep(line,'[\n\r]+','');
                
                if strcmp( line(1:1),'*' )
                    break;
                end
                
                n_IC = n_IC + 1 ;
                
                data = strsplit(line,',');
                
                temp = data{1};
                set_name = strrep(temp,'-','_');
                
                GEO.IC(n_IC).set = set_name;
                GEO.IC(n_IC).dof = str2double (data(2));
                
                GEO.IC(n_IC).value = str2double (data(3));
                
                GEO.IC(n_IC).type = type;
                
            end
        elseif  strcmp(type , 'TEMPERATURE')
            while(true)
                
                line = upper (fgets(fid));
                line = strrep(line,' ','');
                line = strrep(line,'\n','');
                line = regexprep(line,'[\n\r]+','');
                
                if strcmp( line(1:1),'*' )
                    break;
                end
                
                n_IC = n_IC + 1 ;
                
                data = strsplit(line,',');
                
                temp = data{1};
                set_name = strrep(temp,'-','_');
                
                GEO.IC(n_IC).set = set_name;
                GEO.IC(n_IC).dof = 11;
                
                GEO.IC(n_IC).value = str2double (data(2));
                
                GEO.IC(n_IC).type = type;
                
            end            
        end

    elseif length(line) > 24 && strcmp( line(1:25),'*ADAPTIVE MESH CONSTRAINT' )
                        
        line            = strrep(line,' ','');
        line            = strrep(line,'\n','');
        data            = strsplit(line,',');
        
        temp2           = data{2};
        temp3           = data{3};
        
        constraint_type = strrep(temp2(16:length(temp2)),'-','_');
        type            = strrep(temp3(6:length(temp3)),'-','_');

        if contains(constraint_type , 'LAGRANGIAN')
            constraint_type = 'LAGRANGIAN';
        else
            error('Invalid ADAPTIVE MESH CONSTRAINT type. Just LAGRANGIAN is valid');
        end         
        
        
        if contains(type , 'VELOCITY')
            type = 'VELOCITY';
        else
            error('Invalid ADAPTIVE MESH CONSTRAINT type. Just velocity type is implemented');
        end
        
        n_ALEC          = n_ALEC + 1;
        
        line            = upper (fgets(fid));
        line            = strrep(line,' ','');
        line            = strrep(line,'\n','');
        line            = regexprep(line,'[\n\r]+','');
        data            = strsplit(line,',');
        temp            = data{1};
        
        set_name        = strrep(temp,'-','_');
        
        
        
        GEO.ALEC(n_ALEC).set               = set_name;
        GEO.ALEC(n_ALEC).type              = type;
        GEO.ALEC(n_ALEC).constraint_type   = constraint_type;
        
        
    elseif length(line) > 13 && strcmp( line(1:14),'*SOLID SECTION' )
        n_Section = n_Section + 1;
        
        line = strrep(line,' ','');
        line = strrep(line,'\n','');
        line = regexprep(line,'[\n\r]+','');
        data = strsplit(line,',');
        
        temp = data{2};
        set_name = strrep(temp(6:length(temp)),'-','_');
        temp2 = data{3};
        material_name = strrep(temp2(10:length(temp2)),'-','_');
        
        DAT.SECTION(n_Section).set = set_name;
        
        DAT.SECTION(n_Section).mat = material_name;
        
        line = upper (fgets(fid));
        
    elseif length(line) > 8 && strcmp( line(1:9),'*MATERIAL' )
        
        n_MAT = n_MAT + 1;
        
        line = strrep(line,' ','');
        line = strrep(line,'\n','');
        line = regexprep(line,'[\n\r]+','');
        data = strsplit(line,',');
        
        temp2 = data{2};
        material_name = strrep(temp2(6:length(temp2)),'-','_');
        
        DAT.MAT(n_MAT).name = material_name;
        
        line = upper (fgets(fid));
        
    elseif length(line) > 7 && strcmp( line(1:8),'*DENSITY' )
        
        line = upper (fgets(fid));
        line = strrep(line,' ','');
        line = strrep(line,'\n','');
        line = regexprep(line,'[\n\r]+','');
        data = strsplit(line,',');
        
        densiy = str2double (data(1));
        
        DAT.MAT(n_MAT).density = densiy;
    elseif length(line) > 10 && strcmp( line(1:11),'*NEOHOOKEAN' )
        
        DAT.MAT(n_MAT).type = 1;
        
        props_value = zeros(1,1);
        
        ns_value = 0;
        
        while(true)
            
            line = upper (fgets(fid));
            
            if strcmp( line(1:1),'*' )
%                 DAT.MAT(n_MAT).props = props_value;
                
                DAT.MAT(n_MAT).props.E = props_value(1);
                DAT.MAT(n_MAT).props.Poisson = props_value(2);                
                
                
                break;
            end
            
            data = strsplit(line,',');
            
            for i=1:length(data)
                ns_value = ns_value + 1;
                props_value(ns_value) = str2double (data(i));
            end
            
        end
        
    elseif length(line) > 7 && strcmp( line(1:8),'*ELASTIC' )
        
        DAT.MAT(n_MAT).type = 3;
        
        props_value = zeros(1,1);
        
        ns_value = 0;
        
        while(true)
            
            line = upper (fgets(fid));
            
            if strcmp( line(1:1),'*' )
%                 DAT.MAT(n_MAT).props = props_value;
                
                DAT.MAT(n_MAT).props.E = props_value(1);
                DAT.MAT(n_MAT).props.Poisson = props_value(2);                
                
                
                break;
            end
            
            data = strsplit(line,',');
            
            for i=1:length(data)
                ns_value = ns_value + 1;
                props_value(ns_value) = str2double (data(i));
            end
            
        end
        
    elseif length(line) > 7 && strcmp( line(1:8),'*PLASTIC' )
        
        DAT.MAT(n_MAT).type = 2;
        
        props_value = zeros(1,1);
        
        ns_value = 0;
        
        while(true)
            
            line = upper (fgets(fid));
            
            if strcmp( line(1:1),'*' )
%                 DAT.MAT(n_MAT).props = props_value;
                
                DAT.MAT(n_MAT).props.E              = props_value(1);
                DAT.MAT(n_MAT).props.Poisson        = props_value(2);                
                DAT.MAT(n_MAT).props.YieldStress    = props_value(3);
                DAT.MAT(n_MAT).props.Hardening      = props_value(4);
                
                
                break;
            end
            
            data = strsplit(line,',');
            
            for i=1:length(data)
                ns_value = ns_value + 1;
                props_value(ns_value) = str2double (data(i));
            end
            
        end
        
        
        
    elseif length(line) > 4 && strcmp( line(1:5),'*STEP' )
        
        n_STEP = n_STEP + 1;
        in_step = 1;
        line = strrep(line,' ','');
        line = strrep(line,'\n','');
        line = regexprep(line,'[\n\r]+','');
        data = strsplit(line,',');
        
        temp2 = data{2};
        step_name = strrep(temp2(6:length(temp2)),'-','_');
        
        DAT.STEP(n_STEP).name = step_name;
        
        line = upper (fgets(fid));
        line = strrep(line,' ','');
        line = strrep(line,'\n','');
        line = regexprep(line,'[\n\r]+','');
        data = strsplit(line,',');
        
        temp2 = data{1};        
        
        step_type = strrep(temp2(2:length(temp2)),'-','_');
        
%         if contains(step_type , 'DYNAMIC')
%             DAT.STEP(n_STEP).type = 1;           
%         else
%             error('Invalid Step type');
%         end
        
        % Computation of time increment
        % [1] Fixed Delta t
        % [2] Fixed CFL
        
        flag_fix = 0;
        CON.Type_dt = 2;
        
        if length (data) > 2
            
            temp2 = data{3}; 
            
            if (temp2(1:6) ==  "DIRECT")
                flag_fix = 1;
                CON.Type_dt = 1;
            end
            
        end
        

        line = upper (fgets(fid));
        line = strrep(line,' ','');
        line = strrep(line,'\n','');
        line = regexprep(line,'[\n\r]+','');
        
        step_value = zeros(1,1);  
        
        data = strsplit(line,',');
        
        if flag_fix == 1
            for i=1:length(data)
                step_value(i) = str2double (data(i));
            end    
            DAT.STEP(n_STEP).Increment = fix(step_value(2)/step_value(1));
        else            
            step_value(1) = str2double (data(2));
            DAT.STEP(n_STEP).Increment = 10000;
        end
        
        DAT.STEP(n_STEP).VALUES = step_value;
        
        DAT.STEP(n_STEP).cfl = cfl;
        DAT.STEP(n_STEP).output_frequency = output_frequency;
        
        line = upper (fgets(fid));        
        
    elseif length(line) > 3 && strcmp( line(1:4),'*CFL' )
        line =  upper (fgets(fid));
                        
        data = strsplit(line,',');
        
        cfl = str2double (data(1));
        
        DAT.STEP(n_STEP).cfl = cfl;
        
        line =  upper (fgets(fid));
        
    elseif length(line) > 9 && strcmp( line(1:10),'*INCREMENT' )
        line =  upper (fgets(fid));
                        
        data = strsplit(line,',');
        
        inc = str2double (data(1));
        
        DAT.STEP(n_STEP).Increment = fix(inc);
        
        line =  upper (fgets(fid));
        
    elseif length(line) > 16 && strcmp( line(1:17),'*OUTPUT FREQUENCY' )
        
        line =  upper (fgets(fid));
        line = strrep(line,' ','');
        line = strrep(line,'\n','');
        line = regexprep(line,'[\n\r]+','');
        data = strsplit(line,',');                        
        
        output_frequency = str2double (data(1));
        
        DAT.STEP(n_STEP).output_frequency = output_frequency;
        
        line =  upper (fgets(fid));        
        
    elseif length(line) > 15 && strcmp( line(1:16),'*TIME INTEGRATOR' )
        
        line = strrep(line,' ','');
        line = strrep(line,'\n','');
        line = regexprep(line,'[\n\r]+','');
        data = strsplit(line,',');
        
        temp2 = data{2};
        
        type = strrep(temp2(6:length(temp2)),'-','_');        

        if (type == 'TRK')
            CON.Time_Integrator = 1;
        else
            error('Invalid Time Integrator type');
        end       
        
        line = upper (fgets(fid)); 
        
        
    elseif length(line) > 15 && strcmp( line(1:16),'*KERNEL FUNCTION' )
        
        line = strrep(line,' ','');
        line = strrep(line,'\n','');
        line = regexprep(line,'[\n\r]+','');
        data = strsplit(line,',');
        
        temp2 = data{2};
        
        type = temp2(6:length(temp2));
        
        if length(type) > 8 && strcmp(type(1:9), 'QUADRATIC')
            CON.Kernel_Type = 2;
        else
            error('Invalid Kernel Function type');
        end        
        
        line = upper (fgets(fid));         
        
    elseif length(line) > 8 && strcmp( line(1:9),'*UPWIND_P' )

        line =  upper (fgets(fid));
        line = strrep(line,' ','');
        line = strrep(line,'\n','');
        line = regexprep(line,'[\n\r]+','');
        data = strsplit(line,',');                        
        
        DAT.Upwind_p = str2double (data(1));
        
        line =  upper (fgets(fid));
        
    elseif length(line) > 8 && strcmp( line(1:9),'*UPWIND_J' )

        line =  upper (fgets(fid));
        line = strrep(line,' ','');
        line = strrep(line,'\n','');
        line = regexprep(line,'[\n\r]+','');
        data = strsplit(line,',');                        
        
        DAT.Upwind_J = str2double (data(1));
        
        line =  upper (fgets(fid)); 
        
    elseif length(line) > 16 && strcmp( line(1:17),'*SMOOTHING LENGTH' )
        % Coefficient of smoothing length
        line            =  upper (fgets(fid));
        line            = strrep(line,' ','');
        line            = strrep(line,'\n','');
        line            = regexprep(line,'[\n\r]+','');
        data            = strsplit(line,',');                        
        
        DAT.h_Coeff     = str2double (data(1));
        
        line            =  upper (fgets(fid));   
        
    elseif length(line) > 4 && strcmp( line(1:5),'*ALFA' )
        % Coefficient of alfa for mesh motion
        line        =  upper (fgets(fid));
        line        = strrep(line,' ','');
        line        = strrep(line,'\n','');
        line        = regexprep(line,'[\n\r]+','');
        data        = strsplit(line,',');                        
        
        DAT.Alfa    = str2double (data(1));
        
        line        =  upper (fgets(fid));        

    elseif length(line) > 11 && strcmp( line(1:12),'*EXECUTATION' )
        % Execution
        % [1]: Background run
        % [2]: Command Window
        
        line    = strrep(line,' ','');
        line    = strrep(line,'\n','');
        line    = regexprep(line,'[\n\r]+','');
        data    = strsplit(line,',');
        
        temp2   = data{1};
        
        type    = temp2(6:length(temp2));
        
        if length(type) > 9 && strcmp(type(1:10), 'BACKGROUND')
            CON.Execution = 1;
        else
            CON.Execution = 2;
        end        
        
        line    = upper (fgets(fid));         
      
        
    else
        line    =  upper (fgets(fid));  
            
    end
       

end

fclose(fid);


GEO.n_BC        = n_BC;       % Number of boundary condition
GEO.n_ALEC      = n_ALEC;     % Number of ALE constraint
GEO.n_IC        = n_IC;

DAT.n_MAT       = n_MAT;      % Number of material
DAT.n_STEP      = n_STEP;     % Number of step
DAT.n_Section   = n_Section;  % Number of materia section


% ALE coordinates
GEO.Xi          = GEO.XP;

end










