import matplotlib.pyplot as plt

simulation_steps = []
cpu_time_used = []

# read the result from results.txt
with open('results1.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        step, time = line.strip().split(", ")
        simulation_steps.append(int(step))
        cpu_time_used.append(float(time))

# draw plot
plt.plot(simulation_steps, cpu_time_used)
plt.xlabel('Simulation Steps')
plt.ylabel('CPU Time Used (seconds)')
plt.title('CPU Time Used vs. Simulation Steps')
# plt.show()
plt.savefig('CPU_Time_vs_Steps.png')