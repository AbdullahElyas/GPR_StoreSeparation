classdef RegressionModel < handle
    properties
        SimData
        Train_Ratio
        Test_Ratio
        XTRAIN_Matrix
        YTRAIN_Matrix
        XTEST_Matrix
        YTEST_Matrix
        GPR_Model
    end

    methods
        function obj = RegressionModel(Train_Ratio, Test_Ratio)

            % Define the Mach and AOA of cases
            Mach = [0.85; 0.85; 0.95; 0.95];
            AOA = [0; 4; 0; 4];

            % load SimData.mat if it exists otherwise create it
            if exist('simData.mat', 'file')
                load('simData.mat');
            else
                % Preallocate cell array to store SimData instances
                simData = cell(length(Mach), 1);
                parfor i = 1:length(Mach)
                    % Create instance of class SimData using Mach and AOA
                    simData{i} = SimData(Mach(i), AOA(i), 0, 26000);
                end

                % Convert cell array to array of SimData objects
                simData = [simData{:}];

                % save the workspace
                save('simData.mat', 'simData');
            end
            obj.SimData = simData;
            obj.Train_Ratio = Train_Ratio;
            obj.Test_Ratio = Test_Ratio;

            % Process the data
            Process_Data(obj);
            % Build Regression Model
            Build_GPR(obj);

        end

        function Process_Data(obj)
            % determine size of obj.SImData
            [~,len_SimData] = size(obj.SimData);
            for i = 1:len_SimData
                [X_train, Y_train, XTest, YTest] = Data_Preprocess_Instance(obj.SimData(i),obj.Train_Ratio,obj.Test_Ratio);
                obj.XTRAIN_Matrix = [obj.XTRAIN_Matrix; X_train];
                obj.YTRAIN_Matrix = [obj.YTRAIN_Matrix; Y_train];
                obj.XTEST_Matrix = [obj.XTEST_Matrix; XTest];
                obj.YTEST_Matrix = [obj.YTEST_Matrix; YTest];
            end

            function [XTrain, YTrain, XTest, YTest] = Data_Preprocess_Instance(simData_instance,Train_Ratio,Test_Ratio)
                data = simData_instance.SixDOF_Processed_Data;
                FORCE_BODY = simData_instance.BodyAxis_Loads.Force;
                Mach = simData_instance.Flow_Mach;
                AOA_FreeStream = simData_instance.Flow_AOA;
                Len_Data = height(data);
                % repmat Mach and AOA
                Mach = repmat(Mach,Len_Data,1);
                AOA_FreeStream = repmat(AOA_FreeStream,Len_Data,1);
                MOMENT_BODY = simData_instance.BodyAxis_Loads.Moment;
                Input_Matrix = [data.X data.Y data.Z data.Alpha data.Sideslip, Mach, AOA_FreeStream];
                Output_Matrix = [FORCE_BODY' MOMENT_BODY'];

                % Prepare the data for training and testing
                rng(1); % For reproducibility
                n = size(Input_Matrix,1);
                idx = randperm(n);
                idxTrain = idx(1:round(n*Train_Ratio));
                idxTest = idx(round(n*Test_Ratio)+1:end);

                Time_Train = data.time(idxTrain);
                Time_Test = data.time(idxTest);
                XTrain = Input_Matrix(idxTrain,:);
                YTrain = Output_Matrix(idxTrain,:);
                XTest = Input_Matrix(idxTest,:);
                YTest = Output_Matrix(idxTest,:);
            end
        end
        function Build_GPR(obj)
            XTrain = obj.XTRAIN_Matrix;
            YTrain = obj.YTRAIN_Matrix;

            % % use gaussian process regression to estimate the forces and moments
            regression_model.Mdl_force_body_x = fitrgp(XTrain,YTrain(:,1));
            regression_model.Mdl_force_body_y = fitrgp(XTrain,YTrain(:,2));
            regression_model.Mdl_force_body_z = fitrgp(XTrain,YTrain(:,3));
            regression_model.Mdl_moment_body_x = fitrgp(XTrain,YTrain(:,4));
            regression_model.Mdl_moment_body_y = fitrgp(XTrain,YTrain(:,5));
            regression_model.Mdl_moment_body_z = fitrgp(XTrain,YTrain(:,6));
            obj.GPR_Model = regression_model;
        end
    end
end