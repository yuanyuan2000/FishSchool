import matplotlib.pyplot as plt

data = [99.221689, 98.385476, 94.646056, 83.9133495, 94.541680, 83.557980]

base = data[0]  
speedup = [base / x for x in data]

labels = ['Single Thread', '2 Threads', '4 Threads', '16 Threads', '32 Threads', '64 Threads']
plt.barh(labels, speedup, color='skyblue')

plt.xlim(0.99, max(speedup) + 0.1)

plt.xlabel('Speedup')
plt.ylabel('Test Case')
plt.title('Speedup Compared to Baseline (Single Thread)')

for i, v in enumerate(speedup):
    plt.text(v, i, " {:.2f}".format(v), va='center', color='black', fontweight='bold')

plt.savefig('speedup_chart_without_opt.png')

plt.show()
