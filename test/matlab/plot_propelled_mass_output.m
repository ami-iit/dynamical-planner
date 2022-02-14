function plot_propelled_mass_output(inputStruct, position, velocity, force, propeller, forceDerivative, t, costValue, elapsedTime, freeFalling, expectedForce)

figure

subplot(2,2,1)
plot(t, position)
hold on
plot(t, freeFalling, '--')
title("x")
ylim([-0.01, 1.1 *inputStruct.x0])

subplot(2,2,2)
plot(t, force)
hold on
plot(t, expectedForce, '--')
title("f")

subplot(2,2,3)
plot(t, position .* force)
title("complementarity")

subplot(2,2,4)
plot(t, propeller)
title(['p (', num2str(costValue), ')'])

sgtitle([inputStruct.complementairity, ' Complementarity, ( ', num2str(elapsedTime), 's)'], 'Interpreter', 'none')

end