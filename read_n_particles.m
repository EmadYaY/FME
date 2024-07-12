

function [n_particles, n_element] = read_n_particles(file_name)

fid=fopen(file_name);

n_particles = 0;
n_element   = 0;

line = upper (fgets(fid));

while ~feof(fid)
    % read line by line
    % skip comment lines
%     disp(line);
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
            
            n_particles = n_particles + 1;
                
        end
    elseif length(line) > 7 && strcmp( line(1:8),'*ELEMENT' )
        
        while(true)
            
            line = upper (fgets(fid));
            
            if strcmp( line(1:1),'*' )
                break;
            end
            
            n_element = n_element + 1;
                
        end        
    else
        line =  upper (fgets(fid));
            
    end
    

end

fclose(fid);

end




