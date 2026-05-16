classdef laminateCLT_700
    properties(GetAccess = 'public', SetAccess = 'public')
        % public read access, public write access.
        myMaterialList; % [E1,E2,G12,v12,thickness] nx5, where n is the number of fiber types,        
        myOrientations; % [...] ply orientations, 1xn, n is number of plies, top to bottom.
        myFibers;       % [...] fiber types, 1xn, n is number of plies.

        myPlyNumber;    % n, where n is the number of plies
        myQ;            % non-rotated Q matrix.
        myQbar;         % rotated Q matrix, Qbar, for each ply;.
        myT;            % associated rotation matrix, T, for each Qbar.
        myZ;            % [1:n+1], ply layer locations, where n is the number of plies.
        myABD;          % computed ABD matrix.
        myE1;           % laminate E1.
        myE2;           % laminate E2.
        myG12;          % laminate G12.
        mynu12;         % laminate nu12.
        myThickness;    % laminate thickness.
    end
    methods

%%%%%%%%%%%%%%%%%%%%%%%%%% USER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% creating the laminate
        function obj = laminateCLT_700(materialList,orientations,fibers)
            % defining the ply,
            obj.myMaterialList = materialList;
            obj.myOrientations = orientations;
            obj.myPlyNumber = length(orientations);

            % ERROR HANDLING:
            if (nargin == 3)
                obj.myFibers = fibers;
                if length(fibers) ~= length(orientations)
                    % fiber length must match orientations
                    fprintf('Warning: fiber list does not match orientations')
                end
            else
                % if layers is not given, then assign 1 to all plies
                obj.myFibers = ones(1,obj.myPlyNumber);
            end
            if length(materialList) ~= 5
                % if materialList is not fully populated
                fprintf('Warning: materialList has missing properties')
            end
        end

%%%%%%% calculating the laminate properties 
      % obj = laminateCLT_700(materialList,orientations,fibers);
        function obj = properties(obj)
            obj = obj.Qmatrices();
            obj = obj.zPositions();
            obj = obj.ABDmatrix();
            obj = obj.stiffness();
        end

%%%%%%% generate a report with laminate properties
      % obj = laminateCLT(materialList,orientations,fibers);
      % obj = obj.properties;
        function obj = ResultsFile(obj)
            % Open file for writing (overwrite if exists)
            filename = 'Results.txt';
            fid = fopen(filename, 'w');
            % laminate properties report
            fprintf(fid,'LIST OF FIBERS AND PROPERTIES\n');
            fprintf(fid,'-----------------------------\n');
            fprintf(fid,'%-12s %-12s %-12s %-12s %-12s %-12s\n',...
                'FIBER','E1','E2','G12','nu12','THICKNESS');
            for i = 1:height(obj.myMaterialList)
                data = obj.myMaterialList(i,:);
                fprintf(fid,'%-12d %-12.0f %-12.0f %-12.0f %-12.2f %-12.4f\n', ...
                    i, data(i,1), data(i,2), data(i,3), data(i,4), data(i,5));
            end

            fprintf(fid,'\n\nLAMINATE ORIENTATIONS AND FIBERS\n');
            fprintf(fid, '--------------------------------\n');
            fprintf(fid,'%-6s %-6s %-8s\n','PLY','FIBER','THETA');
            for i = 1:obj.myPlyNumber
                fprintf(fid,'%-6d %-6d %-8.0f\n', i, obj.myFibers(i), obj.myOrientations(i));
            end

            fprintf(fid,'\n\nLAMINATE PROPERTIES\n');
            fprintf(fid, '-------------------\n');
            fprintf(fid,'%-12s %-12s %-12s %-12s %-12s\n',...
            'E1','E2','G12','nu12','THICKNESS');
            fprintf(fid,'%-12.0f %-12.0f %-12.0f %-12.2f %-12.4f\n', ...
            obj.myE1,obj.myE2,obj.myG12,obj.mynu12,obj.myThickness);
            fprintf(fid,'\n\nZ-POSITIONS\n');
            fprintf(fid,'-----------\n');
            for i = 1:obj.myPlyNumber+1
                fprintf(fid,'%-12.5f\n',obj.myZ(i));
            end
            fprintf(fid,'\n\nABD MATRIX\n');
            fprintf(fid, '----------\n');
            fprintf(fid, '%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n', obj.myABD');

            fprintf(fid,'\n\nCOMPUTED PLY MATRICES: Q, T, QBAR\n');
            fprintf(fid,'---------------------------------\n');
            for i = 1:obj.myPlyNumber
                fprintf(fid, 'PLY %d:\n', i);
                fprintf(fid, 'Q:\n');
                fprintf(fid, '%8.0f %8.0f %8.0f\n', obj.myQ{i}');
                fprintf(fid, '\nT:\n');
                fprintf(fid, '%8.2f %8.2f %8.2f\n', obj.myT{i}');
                fprintf(fid, '\nQBAR:\n');
                fprintf(fid, '%8.0f %8.0f %8.0f\n', obj.myQbar{i}');
                fprintf(fid, '\n--------------------------\n\n');
            end

            % Close file
            fclose(fid);
            disp('Data successfully written to Results.txt');
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUB-FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Q, T, and Qbar matrices
        function obj = Qmatrices(obj)
            % allocating empty space
            obj.myQ = cell(1,obj.myPlyNumber);
            obj.myT = cell(1,obj.myPlyNumber);
            obj.myQbar = cell(1,obj.myPlyNumber);
           
            % calling ply properties from material list
            for i = 1:obj.myPlyNumber            
                properties = obj.myMaterialList(obj.myFibers(i),1:4);
                E1 = properties(1);
                E2 = properties(2);
                G12 = properties(3);
                v12 = properties(4);

            % Q matrix for each ply
                v21 = v12*E2/E1;
                Q11 = E1/(1-v12*v21);
                Q12 = v12*E2/(1-v12*v21);
                Q22 = E2/(1-v12*v21);
                Q66 = G12;
                obj.myQ{i} = [Q11 Q12  0 ;
                              Q12 Q22  0 ;
                              0   0  Q66];

            % T and Qbar matrix for each ply
            theta = obj.myOrientations(i)*pi/180; % deg -> rad
            obj.myT{i} = [cos(theta)^2,sin(theta)^2, 2*sin(theta)*cos(theta);
                          sin(theta)^2,cos(theta)^2,-2*sin(theta)*cos(theta);
                          -sin(theta)*cos(theta),sin(theta)*cos(theta), ...
                          (cos(theta)^2 - sin(theta)^2)];
            obj.myQbar{i} = obj.myT{i}\obj.myQ{i}*transpose(inv(obj.myT{i}));
            end
        end

%%%%%%% Z vector
        function obj = zPositions(obj)
            % allocating empty space
            obj.myZ = zeros(1,obj.myPlyNumber+1);
           
            n = obj.myPlyNumber + 1;
            bigZ = zeros(1,n);
            for i = 1:n-1
                % calling thickness from myMaterialList
                thickness = obj.myMaterialList(obj.myFibers(i),5);
                bigZ(i+1) = thickness + bigZ(i);
            end
            obj.myZ = bigZ - bigZ(n)/2;
        end

%%%%%%% ABD matrix
        function obj = ABDmatrix(obj)
            % allocating empty space
            obj.myABD = zeros(6);
            A = zeros(3,3);
            B = zeros(3,3);
            D = zeros(3,3);

            % calculating A, B, and D 3x3 matrices
            for i = 1:obj.myPlyNumber
                A = A + obj.myQbar{i}*(obj.myZ(i+1) - obj.myZ(i));
                B = B + 1/2*obj.myQbar{i}*(obj.myZ(i+1)^2 - obj.myZ(i)^2);
                D = D + 1/3*obj.myQbar{i}*(obj.myZ(i+1)^3 - obj.myZ(i)^3);
            end
            obj.myABD = [A B; B D];
        end

%%%%%%% laminate properties and total thickness
        function obj = stiffness(obj)
            % properties derived from inverse ABD matrix
            abd = inv(obj.myABD);
            obj.myE1 = 1/abd(1,1)/obj.myZ(end)/2;
            obj.myE2 = 1/abd(2,2)/obj.myZ(end)/2;
            obj.myG12 = 1/abd(3,3)/obj.myZ(end)/2;
            obj.mynu12 = -abd(1,2)/abd(1,1);
            obj.myThickness = obj.myZ(end) - obj.myZ(1);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%% NOT FUNCTIONAL %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end