data = load('Double_Slit_100x_0.04mm_slow_Treated.txt'); % Load data into a variable
x_data = data(:, 1);
y_data = data(:, 2);
x_error = data(:, 3);
y_error = data(:, 4);

initial_guess = [1.24, 1200, 190, 0, 0, 0]; % Initial parameter guesses

% Perform curve fitting with increased MaxIter
options = optimset('MaxIter', 350);  % Adjust the value as needed
fitted_params = lsqcurvefit(@double_slit_fit, initial_guess, x_data, y_data, [], [], options);

% Calculate uncertainties
[~, R, J] = nlinfit(x_data, y_data, @double_slit_fit, fitted_params);
ci = nlparci(fitted_params, R, 'Jacobian', J);
uncertainties = diff(ci) / 2;

% Calculate residuals
residuals = y_data - double_slit_fit(fitted_params, x_data);

% Define x_fit and y_fit with the same length
x_fit = linspace(min(x_data), max(x_data), 10000);
y_fit = double_slit_fit(fitted_params, x_fit);

% Plot the data with error bars and the fitted curve
figure;
errorbar(x_data, y_data, y_error, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
hold on;
plot(x_fit, y_fit, 'r', 'LineWidth', 2);
legend('Data','Best Fit')
xlabel('Position along pattern (m)');
ylabel('Intensity (W/m^2)');
title('Double Slit Interference and Diffraction Pattern Intensity vs Position Along Pattern');

% Plot residuals with error bars
figure
errorbar(x_data, residuals, y_error, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
hold on;
plot([min(x_data), max(x_data)], [0, 0], 'r', 'LineWidth', 2); % Add a red horizontal line at y = 0
xlabel('Position along pattern (m)');
ylabel('Residuals');
title('Residuals vs Position');
for i = 1:numel(fitted_params)
    fprintf('Parameter %d: %.4f ± %.4f\n', i, fitted_params(i), uncertainties(i));
end

saveas(gcf, 'your_figure.png');

% Calculate chi-square statistic
chi_square = sum((residuals ./ y_error).^2);

% Calculate R-squared
y_mean = mean(y_data);
ss_total = sum((y_data - y_mean).^2);
ss_residual = sum(residuals.^2);
r_squared = 1 - (ss_residual / ss_total);

% Root Mean Square Error (RMSE)
rmse = sqrt(mean(residuals.^2));


fprintf('Chi-Square: %.4f\n', chi_square);
fprintf('R-squared: %.4f\n', r_squared);
fprintf('RMSE: %.4f\n', rmse);

% QQ Plot
figure
qqplot(residuals);

% Example fitting function with parameters 'a', 'b', 'c', and 'd'
function y_fit = double_slit_fit(params, x)
    a = params(1);
    b = params(2);
    c = params(3);
    d = params(4);
    e = params(5);
    f = params(6);
    y_fit = a * (cos(b * x + e).^2) .* (sinc(c * x + f).^2) + d;
end


