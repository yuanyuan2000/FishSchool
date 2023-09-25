import matplotlib.pyplot as plt

data = [83.9135995, 45.992330, 19.791143 ]

base = data[0]
speedup = [base / x for x in data]

labels = ['Base', 'Stack', 'Stack And O3']
plt.barh(labels, speedup, color='skyblue')

plt.xlim(0.98, max(speedup)+0.3)

plt.xlabel('Speedup')
plt.ylabel('Test Case (run with 16 Threads)')
plt.title('Speedup Compared to Baseline')

for i, v in enumerate(speedup):
    plt.text(v, i, " {:.2f}".format(v), va='center',
             color='black', fontweight='bold')

plt.savefig('speedup_chart_stack_opt.png')

plt.show()
