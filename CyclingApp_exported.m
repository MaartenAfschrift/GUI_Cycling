classdef CyclingApp_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                  matlab.ui.Figure
        clearButton               matlab.ui.control.Button
        ReadyLamp                 matlab.ui.control.Lamp
        ReadyLampLabel            matlab.ui.control.Label
        NameEditField             matlab.ui.control.EditField
        NameEditFieldLabel        matlab.ui.control.Label
        freqSlider                matlab.ui.control.Slider
        freqrotminSliderLabel     matlab.ui.control.Label
        CranklengthSliderLabel_3  matlab.ui.control.Label
        SaddleRySlider            matlab.ui.control.Slider
        LogTextArea               matlab.ui.control.TextArea
        LogTextAreaLabel          matlab.ui.control.Label
        VisualiseButton           matlab.ui.control.Button
        SimulateButton            matlab.ui.control.Button
        CranklengthSliderLabel_2  matlab.ui.control.Label
        SaddleRxSlider            matlab.ui.control.Slider
        CranklengthmSlider        matlab.ui.control.Slider
        CranklengthmSliderLabel   matlab.ui.control.Label
        Visualisation             matlab.ui.control.UIAxes
        ScoreAxis                 matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        SimResults % structure with simulation results
        ScoreVect % vector that keeps track of results
        resultsdir % directory with restuls
        savename % output name 
        res_store % Description
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function start_cyclingApp(app)
            % init of cycling app

            % add path to simulation code
            % MainFolderPath = which('CyclingApp.mlapp');
            MainFolderPath = which('CyclingApp_exported.m'); % or point to readme
            [MainPath, MainFolder] = fileparts(fileparts(MainFolderPath));
            MainPath = fullfile(MainPath, MainFolder);
            addpath(genpath(fullfile(MainPath)));
            app.LogTextArea.Value = {'Start GUI'};
            app.resultsdir = fullfile(MainPath,'res');
            if ~isfolder(app.resultsdir)
                mkdir(app.resultsdir);
            end

            % add the casadipath
            install_casadi();
            % CasFilePath = which('casadiMEX.mexw64');
            % [CasPath,~,~] = fileparts(CasFilePath);
            % addpath(genpath(CasPath));
            
            if exist(fullfile(app.resultsdir, 'allres.mat'),'file')
                res_prev_sessions = load(fullfile(app.resultsdir, 'allres.mat'),'res_store');
                res_store1 = res_prev_sessions.res_store;
                
            else
                res_store1.res_mat = [];
                res_store1.res_header = {'Power','crank_length','saddle_rx','saddle_ry','freq'};
                res_store1.res_names = {};
            end
            app.res_store = res_store1;


            %  % model properties
            % app.model.crank_length = app.CranklengthmSlider.Value; % crank length in m
            % app.model.saddle_coord = [app.SaddleRxSlider.Value, ...
            %     app.SaddleRySlider.Value]; % x and y position crankc w.r.t bb
            % app.model.freq = [120]; % rotations per minute 

            % init visualisation based on default settings
            % find feasible solution with cranck vertical 
            % and muscle fiber lengths (assuming rigid tendons)
            % closest to optimal length

            % current fix is to run simulation a startup and visualise this
            % result
            % app.SimulateButtonPushed()
            % app.VisualiseSim()



        end

        % Value changed function: CranklengthmSlider
        function CranklengthmSliderValueChanged(app, event)
            value = app.CranklengthmSlider.Value;

        end

        % Button pushed function: SimulateButton
        function SimulateButtonPushed(app, event)
            % Pushed button run simulation

            % select muscle model (contraction dynamics)
            S.MuscleModel = 'Leuven'; % options are (1) VU, (2) Leuven
            S.specific_tension = 25; % 25 N/cm2

            % select MSK model
            S.MSK_Model = 'Kistemaker'; % (1) model Kistemaker et al ()

            % Task settings
            S.Cycling.cf = app.freqSlider.Value/60; % cycling frequency [Hz]
            S.Cycling.isokin = true;
            S.Cycling.Saddle = [app.SaddleRxSlider.Value, ...
                app.SaddleRySlider.Value]; % x and y position crankc w.r.t bb
            S.Cycling.FixedPower = NaN; % impose net power generation [note power of 1 leg]
            S.Cycling.Opt_cf = false; % True if you want to opimize pedalling frequency (instead of imposing it)
            S.Cycling.CrankLength = app.CranklengthmSlider.Value;
            % Time discretization
            S.Coll.scheme = 'Trapezoidal'; % options are: (1) Trapezoidal
            S.Coll.N = 50; % mesh intervals

            % objective function
            % options are (1) Multi_a_E, (2) stim, (3) maxpower, (4) MinNegFiberWork
            S.Objective.type = {'Multi_a_E'};
            S.Objective.w_metab = 0.1;
            S.Objective.w_stim = 100;
            S.Objective.w_qdd = 0.001;
            S.Objective.w_vMtilde = 0.0001;
            S.Objective.w_M1 = 0.0001;
            S.Objective.scale = 10^-2;
            S.Objective.C_forces = 10^-5;
            S.Objective.w_cranckP = 0.1;
            S.Objective.w_minNegWork = 10;

            % scaling of optimization problem
            scaling.fi = 1;
            scaling.fid = 10;
            scaling.fidd = 1000;
            scaling.Fr5x = 100;
            scaling.Fr5y = 400;
            scaling.M1 = 100;
            S.scaling = scaling;

            % settings metabolic energy equations
            S.Metab.scaleRate = 1; % scaling heat rate coeffs
            S.Metab.Model = 'Bhargava2004'; % options are (1) Bhargava2004, (2) Umberger2003 (3) Umberger2010 (4) Umberger2016

            % Solver settings
            SolverSetup.nlp.solver = 'ipopt';
            SolverSetup.derivatives.derivativelevel = 'second';
            %     SolverSetup.optionssol.ipopt.hessian_approximation = 'limited-memory';
            SolverSetup.optionssol.ipopt.nlp_scaling_method = 'none';
            SolverSetup.optionssol.ipopt.linear_solver = 'mumps'; %(ma57 is a bit faster for default problem)
            SolverSetup.optionssol.ipopt.tol = 1e-4;
            SolverSetup.optionssol.ipopt.max_iter = 1000;
            % SolverSetup.optionssol.detect_simple_bounds = true;
            SolverSetup.optionssol.expand = true; % expands MX variabels to SX (faster)
            S.SolverSetup = SolverSetup;

            % plot bool
            S.BoolPlot = false;

            % print to log info
            app.LogTextArea.Value = [app.LogTextArea.Value; {'   '}];
            app.LogTextArea.Value = [app.LogTextArea.Value; {'Start simulation with settings:'}];
            app.LogTextArea.Value = [app.LogTextArea.Value; {['    crank length ' num2str(S.Cycling.CrankLength) ' m'] }];
            app.LogTextArea.Value = [app.LogTextArea.Value; {['    saddle rx ' num2str(S.Cycling.Saddle(1)) ' m']}];
            app.LogTextArea.Value = [app.LogTextArea.Value; {['    saddle ry ' num2str(S.Cycling.Saddle(2)) ' m']}];
            app.LogTextArea.Value = [app.LogTextArea.Value; {['    freq' num2str(S.Cycling.cf*60) ' rot/min'] }];
            app.LogTextArea.Value = [app.LogTextArea.Value; {'    ....simulation running....'}];


            % Run simulations
            S.Objective.type = {'maxpower'};% options are (1) Multi_a_E, (2) stim, (3) maxpower, (4) MinNegFiberWork
            try
                app.ReadyLamp.Color = 'red';
                [app.SimResults] = PredSim_Cycling(S);

                % mean power
                mean_power = mean(app.SimResults.CranckPower,'omitnan')*2;

                % print start simulation
                if app.SimResults.stats.success

                    app.LogTextArea.Value = [app.LogTextArea.Value; {'Simulation finished'}];
                    app.LogTextArea.Value = [app.LogTextArea.Value;...
                        { ['average mechanical power is : ' num2str(mean_power) 'W']}];
                    if ~isempty(app.res_store.res_mat)
                        [max_power, imax] = max(app.res_store.res_mat(:,1));
                        record_holder = app.res_store.res_names{imax};
                    else
                        max_power = 0;
                    end
                    
                    if mean_power>max_power
                        app.LogTextArea.Value = [app.LogTextArea.Value; {'New Record!'}];
                    else
                        app.LogTextArea.Value = [app.LogTextArea.Value; {[' record is ' num2str(max_power) ' W from ', record_holder]}];
                    end


                    % update score vect
                    if isempty(app.ScoreVect)
                        app.ScoreVect = mean_power;
                    else
                        app.ScoreVect = [app.ScoreVect, mean_power];
                    end
                    resdatvis(1) = plot(app.ScoreAxis,app.ScoreVect,'o','MarkerSize',6,'Color',[0 0 1],'MarkerFaceColor',[0 0 1]);
                    hold(app.ScoreAxis,'on')
                    % add max result
                    xrange = [1 length(app.ScoreVect)];
                    if ~isempty(app.res_store.res_mat)
                        [max_power, imax] = max(app.res_store.res_mat(:,1));
                        min_power = min(app.res_store.res_mat(:,1));
                        resdatvis(2) = plot(app.ScoreAxis,xrange,[max_power max_power],'--k');
                        % add text with max power
                        resdatvis(3) = text(app.ScoreAxis,mean(xrange), ...
                            max_power + (max_power-min_power)*0.02, ...
                            app.res_store.res_names{imax}, ...
                            'Color', [0 0 0], ...
                            'FontSize', 10, 'FontWeight', 'bold');

                    end
                    hold(app.ScoreAxis,'off');

                else
                    app.LogTextArea.Value = [app.LogTextArea.Value; {'Optimizer failed'}];
                    app.LogTextArea.Value = [app.LogTextArea.Value; {'My worst nightmare just occured :), the optimization did not converge ...'}];
                    app.LogTextArea.Value = [app.LogTextArea.Value; {'Please adjust the model parameters, execute the simulation again'}];
                end

                % save results
                % adjust the savename if this file already exists
                app.savename = app.NameEditField.Value;
                if exist(fullfile(app.resultsdir,[app.savename '.mat']),'file')
                    NameNew = app.savename;
                    ct = 1;
                    while exist(fullfile(app.resultsdir,[NameNew '.mat']),'file')
                        NameNew = [app.savename '_' num2str(ct)];
                        ct = ct+1;
                    end
                    app.LogTextArea.Value = [app.LogTextArea.Value; {[app.savename ' was already used, file will be saved as ' NameNew]}];
                    app.savename = NameNew;                    
                end

                % save results to a matlfile
                simres = app.SimResults;
                save(fullfile(app.resultsdir,[app.savename '.mat']),'simres')

                % update structure with results
                if isempty(app.res_store.res_mat)
                    ct = 1;
                else
                    ct = length(app.res_store.res_mat(:,1))+1;
                end
                app.res_store.res_mat(ct,1) = mean_power;
                app.res_store.res_mat(ct,2) = app.CranklengthmSlider.Value;
                app.res_store.res_mat(ct,3) = app.SaddleRxSlider.Value;
                app.res_store.res_mat(ct,4) = app.SaddleRySlider.Value;
                app.res_store.res_mat(ct,5) = app.freqSlider.Value;
                app.res_store.res_names{ct} = app.savename;

                % maybe save res_store
                res_store = app.res_store;
                save(fullfile(app.resultsdir, 'allres.mat'),'res_store');
            catch
                app.SimResults = [];
                app.LogTextArea.Value = [app.LogTextArea.Value; {'Optimizer failed'}];
                app.LogTextArea.Value = [app.LogTextArea.Value; {'My worst nightmare just occured :), the optimization did not converge ...'}];
                app.LogTextArea.Value = [app.LogTextArea.Value; {'Please adjust the model parameters, execute the simulation again'}];
            end
            app.ReadyLamp.Color = 'green';


        end

        % Button pushed function: VisualiseButton
        function VisualiseSim(app, event)
            % update app.Visualisation here
            if ~isempty(app.SimResults)

                % rJoints is 2x5xN: [x;y] x joints (crank→foot→lower leg→upper leg→saddle) x frames
                rJoints = app.SimResults.rJoint;
                time = app.SimResults.t;
                nFrames = size(rJoints, 3);

                % Axis limits
                % allX = reshape(rJoints(1,:,:), [], 1);
                % allY = reshape(rJoints(2,:,:), [], 1);
                % margin = 0.1 * max([range(allX), range(allY)]);
                xLimits = [-0.5 0.5];
                yLimits = [-0.3  1.1];

                % Colors: crank (black), foot (dark green), lower leg (blue), upper leg (red)
                segmentColors = {'k', [0 0.5 0], 'b', 'r'};
                segmentNames = {'Crank', 'Foot', 'Lower Leg', 'Upper Leg'};

                % Setup figure
                ax = app.Visualisation;
                app.Visualisation.DataAspectRatio = [1 1 1];
                %grid on;
                % xlabel('X (m)');
                % ylabel('Y (m)');
                % title('VU Cycling Model');
                app.Visualisation.XLim = xLimits;
                app.Visualisation.YLim = yLimits;
                set(app.Visualisation, 'XLimMode', 'manual', 'YLimMode', 'manual');

                % Draw segments (4 lines between 5 joints)
                hLinks = gobjects(1, 4);
                hold(app.Visualisation,'off')
                for i = 1:4
                    hLinks(i) = plot(ax, [rJoints(1,i,1), rJoints(1,i+1,1)], ...
                        [rJoints(2,i,1), rJoints(2,i+1,1)], ...
                        'o-', 'Color', segmentColors{i}, ...
                        'LineWidth', 3, 'MarkerSize', 6, ...
                        'MarkerFaceColor', segmentColors{i}); 
                    if i == 1
                        hold(app.Visualisation,'on')
                    end
                end
                

                % Highlight saddle (last joint)
                hSaddle = plot(ax, rJoints(1,5,1), rJoints(2,5,1), 's', ...
                    'MarkerSize', 10, 'MarkerEdgeColor', 'k', ...
                    'MarkerFaceColor', 'y');

                % Optional: draw dashed crank circle (based on crank length)
                crankRadius = norm(rJoints(:,1,1) - rJoints(:,2,1));
                theta = linspace(0, 2*pi, 100);
                crankCircle = plot(ax,...
                    rJoints(1,1,1) + crankRadius*cos(theta), ...
                    rJoints(2,1,1) + crankRadius*sin(theta), ...
                    '--', 'Color', [0.6 0.6 0.6]);

                % Optional: text labels (fixed positions)
                hLabels = gobjects(1, 4);
                for i = 1:4
                    hLabels(i) = text(ax,mean([rJoints(1,i,1), rJoints(1,i+1,1)]), ...
                        mean([rJoints(2,i,1), rJoints(2,i+1,1)]), ...
                        ['  ', segmentNames{i}], ...
                        'Color', segmentColors{i}, ...
                        'FontSize', 10, 'FontWeight', 'bold');
                end

                % Animation
                dt = mean(diff(time),'omitnan');
                
                t_main = tic;
                for frame = 1:nFrames
                    % Update links
                    t_local = tic;
                    for i = 1:4
                        set(hLinks(i), ...
                            'XData', [rJoints(1,i,frame), rJoints(1,i+1,frame)], ...
                            'YData', [rJoints(2,i,frame), rJoints(2,i+1,frame)]);
                        set(hLabels(i), ...
                            'Position', mean([rJoints(1:2,i,frame), rJoints(1:2,i+1,frame)], 2));
                    end

                    % Update saddle
                    %set(hSaddle, 'XData', rJoints(1,5,frame), 'YData', rJoints(2,5,frame));

                    drawnow limitrate;
                    dt_clmp = toc(t_local);
                    pause((dt-dt_clmp)*0.92); % 0.92 needed because calling pause function als takes some time apparantly
                end
                dt_vis = toc(t_main);
                % disp(['visualisation time is', num2str(dt_vis)]);
                % disp(['actual time cycle is ', num2str(time(end))]);
            end



        end

        % Callback function
        function SaveResultsPushed(app, event)
            
            % ID for URL to spread sheet
            % spreadsheetID = '1ZUZ-RShKnZ43-_WUhOOVf0SlH6AwzJSM7-rFp5P3bOc';
            spreadsheetID =   '1_NEsSKkFNPmdqTsbI4dVG6Ddh6lEwRDihkeiFRqP0PY';
            % read content of google sheet to not overwrite existing info
            sheetcontent = GetGoogleSpreadsheet(spreadsheetID);
            [nrow,~] = size(sheetcontent);
            
            % check if the last name is equal to
            namesheet = app.savename;
            if strcmp(namesheet,sheetcontent{nrow,1})
                % Create a unique name
                ct = 1;
                while strcmp(namesheet,sheetcontent{nrow,1})
                    NameNew = [namesheet '_' num2str(ct)];
                    ct = ct+1;
                end
                namesheet = NameNew;
            end     
            datainput = cell(1,6);
            sheetID = '0';
            pos = [nrow+1,1];
            
            datainput{1} = namesheet;
            mean_power = mean(app.SimResults.CranckPower,'omitnan')*2;
            datainput{2} = round(mean_power,4);
            datainput{3} = app.CranklengthmSlider.Value;
            datainput{4} = app.SaddleRxSlider.Value;
            datainput{5} = app.SaddleRySlider.Value;
            datainput{6} = app.freqSlider.Value;

            
            mat2sheets(spreadsheetID, sheetID, pos, datainput);
            app.LogTextArea.Value = [app.LogTextArea.Value; {['Your results were printed on the google sheets as ' S.namesheet]}];

        end

        % Button pushed function: clearButton
        function clear_text(app, event)
            app.LogTextArea.Value = {''};
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [1 1 1];
            app.UIFigure.Position = [100 100 643 543];
            app.UIFigure.Name = 'MATLAB App';

            % Create ScoreAxis
            app.ScoreAxis = uiaxes(app.UIFigure);
            title(app.ScoreAxis, 'Average mechanical power')
            xlabel(app.ScoreAxis, 'attempts')
            ylabel(app.ScoreAxis, 'power [W]')
            zlabel(app.ScoreAxis, 'Z')
            app.ScoreAxis.GridLineWidth = 0.5;
            app.ScoreAxis.MinorGridLineWidth = 0.5;
            app.ScoreAxis.Position = [328 392 298 139];

            % Create Visualisation
            app.Visualisation = uiaxes(app.UIFigure);
            title(app.Visualisation, 'VU Cycling Model')
            xlabel(app.Visualisation, 'X [m]')
            ylabel(app.Visualisation, 'Y [m]')
            zlabel(app.Visualisation, 'Z')
            app.Visualisation.GridLineWidth = 0.5;
            app.Visualisation.MinorGridLineWidth = 0.5;
            app.Visualisation.Position = [16 95 286 395];

            % Create CranklengthmSliderLabel
            app.CranklengthmSliderLabel = uilabel(app.UIFigure);
            app.CranklengthmSliderLabel.HorizontalAlignment = 'right';
            app.CranklengthmSliderLabel.Position = [329 353 93 22];
            app.CranklengthmSliderLabel.Text = 'Crank length [m]';

            % Create CranklengthmSlider
            app.CranklengthmSlider = uislider(app.UIFigure);
            app.CranklengthmSlider.Limits = [0.01 0.25];
            app.CranklengthmSlider.ValueChangedFcn = createCallbackFcn(app, @CranklengthmSliderValueChanged, true);
            app.CranklengthmSlider.Position = [443 362 150 3];
            app.CranklengthmSlider.Value = 0.175;

            % Create SaddleRxSlider
            app.SaddleRxSlider = uislider(app.UIFigure);
            app.SaddleRxSlider.Limits = [-0.2 0.2];
            app.SaddleRxSlider.Position = [443 312 150 3];
            app.SaddleRxSlider.Value = -0.05;

            % Create CranklengthSliderLabel_2
            app.CranklengthSliderLabel_2 = uilabel(app.UIFigure);
            app.CranklengthSliderLabel_2.HorizontalAlignment = 'right';
            app.CranklengthSliderLabel_2.Position = [346 302 76 22];
            app.CranklengthSliderLabel_2.Text = 'Saddle rx [m]';

            % Create SimulateButton
            app.SimulateButton = uibutton(app.UIFigure, 'push');
            app.SimulateButton.ButtonPushedFcn = createCallbackFcn(app, @SimulateButtonPushed, true);
            app.SimulateButton.Position = [52 57 88 23];
            app.SimulateButton.Text = 'Simulate';

            % Create VisualiseButton
            app.VisualiseButton = uibutton(app.UIFigure, 'push');
            app.VisualiseButton.ButtonPushedFcn = createCallbackFcn(app, @VisualiseSim, true);
            app.VisualiseButton.Position = [179 57 87 23];
            app.VisualiseButton.Text = 'Visualise';

            % Create LogTextAreaLabel
            app.LogTextAreaLabel = uilabel(app.UIFigure);
            app.LogTextAreaLabel.HorizontalAlignment = 'right';
            app.LogTextAreaLabel.Position = [602 70 25 22];
            app.LogTextAreaLabel.Text = 'Log';

            % Create LogTextArea
            app.LogTextArea = uitextarea(app.UIFigure);
            app.LogTextArea.Position = [316 7 311 148];

            % Create SaddleRySlider
            app.SaddleRySlider = uislider(app.UIFigure);
            app.SaddleRySlider.Limits = [0.6 1.1];
            app.SaddleRySlider.Position = [443 267 150 3];
            app.SaddleRySlider.Value = 0.87;

            % Create CranklengthSliderLabel_3
            app.CranklengthSliderLabel_3 = uilabel(app.UIFigure);
            app.CranklengthSliderLabel_3.HorizontalAlignment = 'right';
            app.CranklengthSliderLabel_3.Position = [346 257 76 22];
            app.CranklengthSliderLabel_3.Text = 'Saddle ry [m]';

            % Create freqrotminSliderLabel
            app.freqrotminSliderLabel = uilabel(app.UIFigure);
            app.freqrotminSliderLabel.HorizontalAlignment = 'right';
            app.freqrotminSliderLabel.Position = [351 205 72 22];
            app.freqrotminSliderLabel.Text = 'freq [rot/min]';

            % Create freqSlider
            app.freqSlider = uislider(app.UIFigure);
            app.freqSlider.Limits = [50 170];
            app.freqSlider.Position = [444 214 150 3];
            app.freqSlider.Value = 90;

            % Create NameEditFieldLabel
            app.NameEditFieldLabel = uilabel(app.UIFigure);
            app.NameEditFieldLabel.HorizontalAlignment = 'right';
            app.NameEditFieldLabel.Position = [51 509 37 22];
            app.NameEditFieldLabel.Text = 'Name';

            % Create NameEditField
            app.NameEditField = uieditfield(app.UIFigure, 'text');
            app.NameEditField.Position = [103 509 183 22];
            app.NameEditField.Value = 'MaartenAfschrift';

            % Create ReadyLampLabel
            app.ReadyLampLabel = uilabel(app.UIFigure);
            app.ReadyLampLabel.HorizontalAlignment = 'right';
            app.ReadyLampLabel.Position = [317 163 50 22];
            app.ReadyLampLabel.Text = 'Ready ?';

            % Create ReadyLamp
            app.ReadyLamp = uilamp(app.UIFigure);
            app.ReadyLamp.Position = [370 164 20 20];

            % Create clearButton
            app.clearButton = uibutton(app.UIFigure, 'push');
            app.clearButton.ButtonPushedFcn = createCallbackFcn(app, @clear_text, true);
            app.clearButton.Position = [590 156 51 20];
            app.clearButton.Text = 'clear';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = CyclingApp_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @start_cyclingApp)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end