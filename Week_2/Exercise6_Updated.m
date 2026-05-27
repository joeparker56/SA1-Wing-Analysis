clear;
close all;
global re_l dudx
dudx = 0 ;

re_l = 1e6;
[theta1, hes1, int1, ils1, itr1, its1,x] = calc_bl(dudx, re_l);

re_l = 1e7;
[theta2, hes2, int2, ils2, itr2, its2, x] = calc_bl(dudx, re_l);

figure(1);
hold on;
p1 = plot(x, theta1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Re = 1e6');
p2 = plot(x, theta2, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Re = 1e7');
add_bl_labels(x, theta1, int1, ils1, itr1, its1, 'b');
add_bl_labels(x, theta2, int2, ils2, itr2, its2, 'r');

hold off;

xlabel('x/L');
ylabel('\theta/L');
legend([p1, p2], 'Location', 'best');
grid on;
print -deps2c Week_2/Exercise6a.eps;

figure(2);
hold on;
p1 = plot(x, hes1, 'b-',  'LineWidth', 1.5, 'DisplayName', 'Re = 1e6');
p2 = plot(x, hes2, 'r-',  'LineWidth', 1.5, 'DisplayName', 'Re = 1e7');
add_bl_labels(x, hes1, int1, ils1, itr1, its1, 'b');
add_bl_labels(x, hes2, int2, ils2, itr2, its2, 'r');
hold off;

xlabel('x/L');
ylabel('He');
legend([p1, p2], 'Location', 'best');
grid on;
print -deps2c Week_2/Exercise6b.eps;


%now add velocity gradient 
dudx = -0.25;

re_l=1e4;
[theta1, hes1, int1, ils1, itr1, its1,x] = calc_bl(dudx, re_l);

re_l=1e5;
[theta2, hes2, int2, ils2, itr2, its2, x] = calc_bl(dudx, re_l);

re_l=1e6;
[theta3, hes3, int3, ils3, itr3, its3, x] = calc_bl(dudx, re_l);

figure(3);
hold on;
p1 = plot(x, theta1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Re = 1e4');
p2 = plot(x, theta2, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Re = 1e5');
p3 = plot(x, theta3, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Re = 1e6');
add_bl_labels(x, theta1, int1, ils1, itr1, its1, 'b');
add_bl_labels(x, theta2, int2, ils2, itr2, its2, 'r');
add_bl_labels(x, theta3, int3, ils3, itr3, its3, 'g');


hold off;

xlabel('x/L');
ylabel('\theta/L');
legend([p1, p2,p3], 'Location', 'best');
grid on;
print -deps2c Week_2/Exercise6c.eps;

figure(4);
hold on;
p1 = plot(x, hes1, 'b-',  'LineWidth', 1.5, 'DisplayName', 'Re = 1e4');
p2 = plot(x, hes2, 'r-',  'LineWidth', 1.5, 'DisplayName', 'Re = 1e5');
p3 = plot(x, hes3, ['g' ...
    '-'],  'LineWidth', 1.5, 'DisplayName', 'Re = 1e6');

add_bl_labels(x, hes1, int1, ils1, itr1, its1, 'b');
add_bl_labels(x, hes2, int2, ils2, itr2, its2, 'r');
add_bl_labels(x, hes3, int3, ils3, itr3, its3, 'g');

hold off;

xlabel('x/L');
ylabel('He');
legend([p1, p2,p3], 'Location', 'best');
grid on;
print -deps2c Week_2/Exercise6d.eps;



function add_bl_labels(x, y, i_nt, i_ls, i_tr, i_ts, color)
    % Plots markers and text labels for boundary layer events.
    % Uses 'HandleVisibility', 'off' so markers don't clutter the legend.
    
    % Natural Transition (Circle)
    if i_nt > 0 && ~isnan(i_nt)
        plot(x(i_nt), y(i_nt), [color 'o'], 'MarkerFaceColor', color, 'HandleVisibility', 'off');
        text(x(i_nt), y(i_nt), '    Natural transition          ', 'Color', color, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment','right');
    end
    
    % Laminar Separation (Square)
    if i_ls > 0 && ~isnan(i_ls)
        plot(x(i_ls), y(i_ls), [color 'o'], 'MarkerFaceColor', color, 'HandleVisibility', 'off');
        text(x(i_ls), y(i_ls), '  Laminar Separation      ', 'Color', color, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment','right');
    end
    
    % Turbulent Reattachment (Triangle)
    if i_tr > 0 && ~isnan(i_tr)
        plot(x(i_tr), y(i_tr), [color 'o'], 'MarkerFaceColor', color, 'HandleVisibility', 'off');
        text(x(i_tr), y(i_tr), '  Turbulent Reattachment       ', 'Color', color, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment','left');
    end
    
    % Turbulent Separation (Cross)
    if i_ts > 0 && ~isnan(i_ts)
        plot(x(i_ts), y(i_ts), [color 'o'], 'MarkerFaceColor',color, 'HandleVisibility', 'off');
        text(x(i_ts), y(i_ts), '  Turbulent Seperation     ', 'Color', color, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment','right');
    end
end