xd = [0.005; 0];
run pid_control

xref = ones(length(out.x),1)*xd(1);

plot(out.t,out.x(:,1))
hold on
plot(out.t,out.x(:,2))
plot(out.t,xref,'--','color',[.5 .5 .5])
legend(["c_x","c_y","desired"])
xlabel("time (s)")
ylabel("c_x (milirad)")
title("PID Closed-Loop Response")
hold off
saveas(gca,"figures/pid_sim.png")