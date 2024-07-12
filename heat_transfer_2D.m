% A plane stress code for FEM
clc

%% Input reading

% Reading mesh information from the input file
%--------------------------------------------------------------------------
% Input file reading
%--------------------------------------------------------------------------
file_name                       = 'heat_transfer_final.inp';

[GEO]                 = read_inp(file_name);


% Defining node set for fixed node and imposing boundary conditions at nodes

set_name = "LEFT";
nset_left = GEO.NSET.(set_name);
%[ 15,  17, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116];
nset_right = GEO.NSET.("RIGHT");
%[14, 16, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68];

%% initialisation
n_node = GEO.NP;    % number of nodes
n_element = GEO.NE;    % number of nodes
n_dof = n_node; % number of degree of freedom

K = zeros(n_dof, n_dof);
RHS = zeros(n_dof, 1);

thikness = 1.0;
k_conduction = 1.0; 
h_convection = 0;
%%
%--------------------------------------------------------------------------
%   Define positions of integration points
%--------------------------------------------------------------------------
xi = zeros(2,4);
xi(1,1) = -0.5773502692;
xi(2,1) = xi(1,1);
xi(1,2) = -xi(1,1);
xi(2,2) = xi(1,1);
xi(1,3) = xi(1,1);
xi(2,3) = -xi(1,1);
xi(1,4) = -xi(1,1);
xi(2,4) = -xi(1,1);

%%

% loop over element
for i_el = 1: n_element
    
    
    conn_el = GEO.CONN(:,i_el);    % element connectivity
    
        % nodal coordinates
        COORDS = zeros(2,4);
        for ii=1:4
            COORDS(:,ii) = GEO.XP(:,conn_el(ii));
        end    
    
    k__T = zeros(4,4);
        
    for igp =1:4 % loop over iontegration points
        
        xi_gp(1) = xi(1,igp);
        xi_gp(2) = xi(2,igp);
        
        dN_rs = shapefunctionderivs(xi_gp);
        N = shapefunction(xi_gp);
        J = COORDS * dN_rs;
        
        J_inv = inv(J);
        
        det_J = det(J);
        
        dN_xy = dN_rs * J_inv;
        
        dN_x = dN_xy(:,1);
        dN_y = dN_xy(:,2);
        
        k__T = k__T + k_conduction * (dN_x * dN_x' + dN_y * dN_y') * det_J * thikness;
        k__T = k__T + 2.0 * h_convection * (N * N') * det_J * thikness;
    end
    
    % Assembly 
    L = [conn_el(1) conn_el(2) conn_el(3) conn_el(4)];
    
    K(L,L) =  K(L,L) + k__T;
   
end

%
% boundary condition
%

% loop over BCs
for i_BC = 1:GEO.n_BC

    set_name = GEO.BC(i_BC).set ;

    set_value = GEO.NSET.(set_name);

    dof = GEO.BC(i_BC).dof;

    value = GEO.BC(i_BC).value;


    for ip=1:length(set_value)
        dof_number = set_value(ip);

        temp_value = K(dof_number,dof_number);

        K(dof_number,:) = 0.0 ;

        for kk=1:n_dof
            RHS(kk) = RHS(kk) - K(kk,dof_number) * value;
            K(kk,dof_number) = 0.0 ;
        end

        K(dof_number,dof_number) = temp_value ;
        RHS(dof_number) = temp_value * value;

    end



end


% Solve
temperature = K\RHS;

x=0:0.5:10;
y=0:0.5:10;

count = 0;
for i=1:21
    for j=1:21
        count = count +1;
        PP(i,j)=temperature(count);
    end
end


[c,h] = contourf(x,y,PP);
clabel(c,h);
colorbar
colormap(jet)
xlabel('X Axis');
ylabel('Y Axis');


%%
%
%================= SHAPE FUNCTION DERIVATIVES ======================
%
function dNdxi = shapefunctionderivs(xi)

  dNdxi = zeros(4,2);


    dNdxi(1,1) = -0.25*(1.-xi(2));
    dNdxi(1,2) = -0.25*(1.-xi(1));
    dNdxi(2,1) = 0.25*(1.-xi(2));
    dNdxi(2,2) = -0.25*(1.+xi(1));
    dNdxi(3,1) = 0.25*(1.+xi(2));
    dNdxi(3,2) = 0.25*(1.+xi(1));
    dNdxi(4,1) = -0.25*(1.+xi(2));
    dNdxi(4,2) = 0.25*(1.-xi(1));

end

function N = shapefunction(xi)

  N = zeros(4,1);


    N(1,1) = 0.25*(1.0-xi(2))*(1.0-xi(1));
    N(2,1) = 0.25*(1.0+xi(1))*(1.0-xi(2));
    N(3,1) = 0.25*(1.0+xi(1))*(1.0+xi(2));
    N(4,1) = 0.25*(1.0-xi(1))*(1.0+xi(2));

end
