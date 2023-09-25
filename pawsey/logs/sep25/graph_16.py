import matplotlib.pyplot as plt

data = [0.00144034375, 0.000091625, 0.000002562,
        0.000000020, 0.000000199]

base = data[0]
speedup = [base / x for x in data]

labels = ['2^6', '2^10', '2^14', '2^18', '2^22']
plt.barh(labels, speedup, color='skyblue')

plt.xlim(0, max(speedup)+15000)

plt.xlabel('Speedup')
plt.ylabel('Test Case (run with 16 Threads)')
plt.title('Speedup (Second per fish) Compared to Baseline (2^6)')

for i, v in enumerate(speedup):
    plt.text(v, i, " {:.2f}".format(v), va='center',
             color='black', fontweight='bold')

plt.savefig('speedup_chart_time_per_fish.png')

plt.show()
