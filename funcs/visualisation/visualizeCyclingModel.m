function visualizeCyclingModel(rJoints, time, varargin)

    % rJoints is 2x5xN: [x;y] x joints (crank→foot→lower leg→upper leg→saddle) x frames
    nFrames = size(rJoints, 3);

    % Axis limits
    allX = reshape(rJoints(1,:,:), [], 1);
    allY = reshape(rJoints(2,:,:), [], 1);
    margin = 0.1 * max([range(allX), range(allY)]);
    xLimits = [min(allX)-margin, max(allX)+margin];
    yLimits = [min(allY)-margin, max(allY)+margin];

    % Colors: crank (black), foot (dark green), lower leg (blue), upper leg (red)
    segmentColors = {'k', [0 0.5 0], 'b', 'r'};
    segmentNames = {'Crank', 'Foot', 'Lower Leg', 'Upper Leg'};

    % Setup figure
    figure;
    ax = gca;
    hold on;
    axis equal;
    grid on;
    xlabel('X (m)');
    ylabel('Y (m)');
    title('VU Cycling Model');
    xlim(ax, xLimits);
    ylim(ax, yLimits);
    set(ax, 'XLimMode', 'manual', 'YLimMode', 'manual');

    % Draw segments (4 lines between 5 joints)
    hLinks = gobjects(1, 4);
    for i = 1:4
        hLinks(i) = plot([rJoints(1,i,1), rJoints(1,i+1,1)], ...
                         [rJoints(2,i,1), rJoints(2,i+1,1)], ...
                         'o-', 'Color', segmentColors{i}, ...
                         'LineWidth', 3, 'MarkerSize', 6, ...
                         'MarkerFaceColor', segmentColors{i});
    end

    % Highlight saddle (last joint)
    hSaddle = plot(rJoints(1,5,1), rJoints(2,5,1), 's', ...
                   'MarkerSize', 10, 'MarkerEdgeColor', 'k', ...
                   'MarkerFaceColor', 'y');    

    % Optional: draw dashed crank circle (based on crank length)
    crankRadius = norm(rJoints(:,1,1) - rJoints(:,2,1));
    theta = linspace(0, 2*pi, 100);
    crankCircle = plot(...
        rJoints(1,1,1) + crankRadius*cos(theta), ...
        rJoints(2,1,1) + crankRadius*sin(theta), ...
        '--', 'Color', [0.6 0.6 0.6]);

    % Optional: text labels (fixed positions)
    hLabels = gobjects(1, 4);
    for i = 1:4
        hLabels(i) = text(mean([rJoints(1,i,1), rJoints(1,i+1,1)]), ...
                          mean([rJoints(2,i,1), rJoints(2,i+1,1)]), ...
                          ['  ', segmentNames{i}], ...
                          'Color', segmentColors{i}, ...
                          'FontSize', 10, 'FontWeight', 'bold');
    end

    % Animation
    dt = nanmean(diff(time));
    for frame = 1:nFrames
        % Update links
        for i = 1:4
            set(hLinks(i), ...
                'XData', [rJoints(1,i,frame), rJoints(1,i+1,frame)], ...
                'YData', [rJoints(2,i,frame), rJoints(2,i+1,frame)]);
            set(hLabels(i), ...
                'Position', mean([rJoints(1:2,i,frame), rJoints(1:2,i+1,frame)], 2));
        end

        % Update saddle
        set(hSaddle, 'XData', rJoints(1,5,frame), 'YData', rJoints(2,5,frame));

        drawnow;
        pause(dt);
    end
end