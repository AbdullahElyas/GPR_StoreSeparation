classdef SimData < handle
    properties
        % Simulation data
        Flow_Mach
        Flow_AOA
        Flow_Sideslip
        Flow_Altitude
        Num_TimeStep
        Size_TimeStep
        Max_Time
        BodyAxis_Loads
        InertialAxs_Loads
        Data_Folder
        Root_Folder = '.\Sim_Data\';
        Pressure_Data_Folder = 'pressure_data';
        SixDOF_Raw_Data
        SixDOF_Processed_Data
        Store_STL

    end

    % build a constructor
    methods
        function obj = SimData(Flow_Mach, Flow_AOA, Flow_Sideslip, Flow_Altitude)

            obj.Flow_Mach = Flow_Mach;
            obj.Flow_AOA = Flow_AOA;
            obj.Flow_Sideslip = Flow_Sideslip;
            obj.Flow_Altitude = Flow_Altitude;
            Mach = regexprep(num2str(obj.Flow_Mach),'\.','');
            AOA  = obj.Flow_AOA;
            FormatSpec = 'M%s_AOA%d';
            Data_Folder = sprintf(FormatSpec,Mach,AOA);
            obj.Data_Folder = Data_Folder;

            obj.LoadStoreSTL();
            obj.ImportSixDOF();
            obj.EstimateStates();
            obj.EvaluateLoads();
        end

        function  ImportSixDOF(obj)
            % remove decimal frpm mach number using regex
            Mach = regexprep(num2str(obj.Flow_Mach),'\.','');

            AOA  = obj.Flow_AOA;
            FormatSpec = 'M%s_AOA%d_store.6dof';
            filename = sprintf(FormatSpec,Mach,AOA);
            full_filename = fullfile(obj.Root_Folder,obj.Data_Folder,filename);

            %% Initialize variables.
            if nargin<=2
                startRow = 6;
                endRow = inf;
            end

            %% Format for each line of text:
            %   column1: double (%f)
            %	column2: double (%f)
            %   column3: double (%f)
            %	column4: double (%f)
            %   column5: double (%f)
            %	column6: double (%f)
            %   column7: double (%f)
            % For more information, see the TEXTSCAN documentation.
            formatSpec = '%12f%13f%13f%13f%13f%13f%f%[^\n\r]';

            %% Open the text file.
            fileID = fopen(full_filename,'r');

            %% Read columns of data according to the format.
            % This call is based on the structure of the file used to generate this code. If an error occurs for a different file, try regenerating the code from the Import Tool.
            dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            for block=2:length(startRow)
                frewind(fileID);
                dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
                for col=1:length(dataArray)
                    dataArray{col} = [dataArray{col};dataArrayBlock{col}];
                end
            end

            %% Close the text file.
            fclose(fileID);

            %% Post processing for unimportable data.
            % No unimportable data rules were applied during the import, so no post processing code is included. To generate code which works for unimportable data, select unimportable cells in a file and regenerate the script.

            %% Create output variable
            obj.SixDOF_Raw_Data = table(dataArray{1:end-1}, 'VariableNames', {'time', 'CG_X', 'CG_Y', 'CG_Z', 'THETA_X', 'THETA_Y', 'THETA_Z'});


        end

        function EstimateStates(obj)
            % Flight Conditions
            flowMach = obj.Flow_Mach;
            flowAltitude = obj.Flow_Altitude*0.3048;
            [~,sound_speed,~,~] = atmosisa(flowAltitude);
            flowAOA = obj.Flow_AOA;
            flowSideslip = obj.Flow_Sideslip;
            flowVel = flowMach*sound_speed;

            %wind flow vector
            flowVector = [flowVel;0;0];

            %Conversion of flow from wind to body axis
            windBetaMatrix = [cosd(flowSideslip) sind(flowSideslip) 0;-sind(flowSideslip) cosd(flowSideslip) 0;0 0 1];
            windAOAMatrix = [cosd(flowAOA) 0 sind(flowAOA); 0 1 0; -sind(flowAOA) 0 cosd(flowAOA)];
            windVel = [windBetaMatrix*windAOAMatrix]'*flowVector;

            %import of data from 6dof file
            data=obj.SixDOF_Raw_Data;
            data.Properties.VariableNames = {'time','X','Y','Z','Roll','Yaw','Pitch'};
            % delete the data with Nan values
            data = data(~any(isnan(data{:,:}),2),:);


            %Estimation of Velocities in inertial axis
            data.Northvel = -[0;diff(data.X)./diff(data.time)];
            data.Downvel = -[0;diff(data.Y)./diff(data.time)];
            data.Eastvel = -[0;diff(data.Z)./diff(data.time)];

            %Conversion of axis and angle to euler angles
            rotation=[ -data.Roll -data.Pitch -data.Yaw  ]./57.3;
            rotmag=((rotation(:,1).^2+rotation(:,2).^2+rotation(:,3).^2).^(1/2));
            % deal with the case of zero rotation
            rotdirectionvec=rotation./rotmag;
            % replace NaN values with 0
            rotdirectionvec(isnan(rotdirectionvec))=0;
            rotm=axang2rotm([rotdirectionvec rotmag]);
            eul = rotm2eul(rotm)*57.3;
            % replace NaN values with 0
            eul(isnan(eul))=0;

            pitch = eul(:,2);
            roll = eul(:,3);
            yaw = eul(:,1);
            E_PitchRate = [0;diff(pitch)./diff(data.time)];
            E_RollRate = [0;diff(roll)./diff(data.time)];
            E_YawRate = [0;diff(yaw)./diff(data.time)];
            EulerRates = [E_RollRate E_PitchRate E_YawRate];

            data.Euler_Pitch = pitch;
            data.Euler_Roll = roll;
            data.Euler_Yaw = yaw;

            for i = 1:height(data)

                %conversion of euler rates to body body angular rates
                euler_to_pqr = [1 0                   -sin(pitch(i)/57.3);
                    0 cos(roll(i)/57.3)  sin(roll(i)/57.3)*cos(pitch(i)/57.3);
                    0 -sin(roll(i)/57.3) cos(roll(i)/57.3)*cos(pitch(i)/57.3)];

                BodyRates(:,i) = euler_to_pqr* EulerRates(i,:)';


                %conversion of inertial velocities to body velocities
                rollMat = [1 0 0;
                    0 cosd(eul(i,3)) sind(eul(i,3));
                    0 -sind(eul(i,3)) cosd(eul(i,3))];

                pitchMat = [cosd(eul(i,2)) 0 -sind(eul(i,2));
                    0 1 0;
                    sind(eul(i,2)) 0 cosd(eul(i,2))];

                yawMat = [cosd(eul(i,1)) sind(eul(i,1)) 0;
                    -sind(eul(i,1)) cosd(eul(i,1)) 0;
                    0 0 1];
                finalrotmatrix = rollMat*pitchMat*yawMat;


                bodyVel(:,i) = finalrotmatrix * [data.Northvel(i);data.Eastvel(i);data.Downvel(i)];

                WindVelCompBody(:,i) = finalrotmatrix*windVel;



            end

            data.PitchRate = BodyRates(2,:)';
            data.RollRate = BodyRates(1,:)';
            data.YawRate  = BodyRates(3,:)';

            airspeedComp = bodyVel'+WindVelCompBody';
            airspeedmag=((airspeedComp(:,1).^2+airspeedComp(:,2).^2+airspeedComp(:,3).^2).^(1/2));

            Alpha = atan(airspeedComp(:,3)./airspeedComp(: ,1))*57.3;
            Sideslip = asin(airspeedComp(:,2)./airspeedmag)*57.3;
            data.Alpha = Alpha;
            data.Sideslip = Sideslip;
            X1 = data.X-data.X(1);
            Y1 = -(data.Z-data.Z(1));
            Z1 = data.Y-data.Y(1);

            data.BodyX = X1;
            data.BodyY = Y1;
            data.BodyZ = Z1;
            obj.SixDOF_Processed_Data = data;
        end

        function LoadStoreSTL(obj)
            % LoadStoreSTL Load and store the STL files
            Full_FileName = fullfile(obj.Root_Folder,'store.stl');
            TR1 = stlread(Full_FileName);
            % TR1 = stlread('box.stl');
            TR.ConnectivityList = TR1.ConnectivityList;
            TR.Points = TR1.Points./1000;
            obj.Store_STL = TR;
        end


        function EvaluateLoads(obj)
            data  = obj.SixDOF_Processed_Data;
            TR = obj.Store_STL;
            % Estimate Body Forces and Moments at each time step and new store CG location with new transformation matrix
            for iter=1:length(data.time)

                % for i=1:200:201
                iter
                TIMESTEP = iter-1;
                XCG = data.X(iter);
                YCG = data.Y(iter);
                ZCG = data.Z(iter);

                %Conversion of axis and angle to euler angles
                rotation=[data.Roll(iter)  data.Yaw(iter) data.Pitch(iter)]./57.3;
                rotmag=((rotation(:,1).^2+rotation(:,2).^2+rotation(:,3).^2).^(1/2));
                % deal with the case of zero rotation
                rotdirectionvec=rotation./rotmag;
                % replace NaN values with 0
                rotdirectionvec(isnan(rotdirectionvec))=0;
                TRANSFORM_MATRIX=axang2rotm([rotdirectionvec rotmag]);
                % if TRANSFORM_MATRIX is Nan then replace it with identity matrix
                if isnan(TRANSFORM_MATRIX)
                    TRANSFORM_MATRIX = eye(3);
                end
                if iter==1
                    index_cell = 0;
                end
                [Force1(:,iter),moment1(:,iter),index_cell] = Estimate_BodyForcemoments(TIMESTEP,XCG,YCG,ZCG,TRANSFORM_MATRIX,TR,index_cell);
                FORCE_BODY(:,iter) = TRANSFORM_MATRIX*Force1(:,iter);
                MOMENT_BODY(:,iter) = TRANSFORM_MATRIX*moment1(:,iter);
            end

            obj.BodyAxis_Loads.Force = FORCE_BODY;
            obj.BodyAxis_Loads.Moment = MOMENT_BODY;
            obj.InertialAxs_Loads.Force = Force1;
            obj.InertialAxs_Loads.Moment = moment1;

            function [Force1,moment1,index_cell] = Estimate_BodyForcemoments(TIMESTEP,XCG,YCG,ZCG,TRANSFORM_MATRIX,TR,index_cell)

                % remove decimal frpm mach number using regex
                Mach = regexprep(num2str(obj.Flow_Mach),'\.','');

                AOA  = obj.Flow_AOA;
                FormatSpec = 'M%s_AOA%d_%04i';
                filename = sprintf(FormatSpec,Mach,AOA,TIMESTEP);

                full_filename = fullfile(obj.Root_Folder,obj.Data_Folder,obj.Pressure_Data_Folder,filename);
                % formatspec = './pressure_data/0,95M_AOA0-%04i';
                % formatspec = './pressure_data_2/M095_AOA0_%04i';
                % filename = sprintf(formatspec,TIMESTEP);
                pressure_cell = readmatrix(full_filename);
                xyz_cell = pressure_cell(:,2:4);
                %  find index of the closest point with points in TR.Points

                for ii = 1:size(TR.ConnectivityList, 1)
                    i = TR.ConnectivityList(ii,1);
                    j = TR.ConnectivityList(ii,2);
                    k = TR.ConnectivityList(ii,3);
                    pointA = TR.Points(i,:);
                    pointB = TR.Points(j,:);
                    pointC = TR.Points(k,:);

                    %calculate the centroid of the triangle
                    centroid(ii,:) = Centroid(pointA,pointB,pointC);
                    centroid(ii,:) = TRANSFORM_MATRIX*centroid(ii,:)';
                    centroid(ii,:) = centroid(ii,:)+ [XCG,YCG,ZCG];
                end

                if mod(TIMESTEP+1,20) == 0 || TIMESTEP == 0
                    index_cell = dsearchn(xyz_cell,centroid);
                end
                pressure_sorted = pressure_cell(index_cell,5);
                Force1 = [0,0,0];
                moment1 = [0,0,0];
                for ii = 1:size(TR.ConnectivityList, 1)
                    i = TR.ConnectivityList(ii,1);
                    j = TR.ConnectivityList(ii,2);
                    k = TR.ConnectivityList(ii,3);
                    A = TR.Points(i,:);
                    B = TR.Points(j,:);
                    C = TR.Points(k,:);


                    [Force_tri,moment_tri] = TriArea(pressure_cell(index_cell(ii),:),A,B,C,XCG,YCG,ZCG,TRANSFORM_MATRIX);
                    Force1 = Force1 + Force_tri;
                    moment1 = moment1 + moment_tri;

                end

                function [Force_tri,moment_tri] = TriArea(pressure_cell_value,A,B,C,XCG,YCG,ZCG,TRANSFORM_MATRIX)

                    % moment center and transformation matrix is updated at each time step
                    moment_center = [XCG,YCG,ZCG];
                    % negative sign so the area normal is pointing in the correct direction (outward direction)
                    Area = -0.5*cross(B-A,C-A);
                    Transform_Area = TRANSFORM_MATRIX*Area';
                    Force_tri = Transform_Area'*pressure_cell_value(5);
                    moment_arm = [pressure_cell_value(2:4)]-moment_center;
                    moment_tri = cross(moment_arm,Force_tri);

                end

                function centroid = Centroid(A,B,C)
                    centroid = (A+B+C)/3;
                end

            end
        end

    end
end