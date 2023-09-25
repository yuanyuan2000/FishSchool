import matplotlib.pyplot as plt

data = [83.9135995, 35.77747975, 31.65407975,
        31.7201316667]

base = data[0]
speedup = [base / x for x in data]

labels = ['Base', 'O_1', 'O_2', 'O_3']
plt.barh(labels, speedup, color='skyblue')

plt.xlim(0.95, max(speedup)+0.3)

plt.xlabel('Speedup')
plt.ylabel('Test Case (run with 16 Threads, 2^22 Fishes)')
plt.title('Speedup Compared to Baseline (Without Compiler Optimization)')

for i, v in enumerate(speedup):
    plt.text(v, i, " {:.5f}".format(v), va='center',
             color='black', fontweight='bold')

plt.savefig('speedup_chart_diff_cmpl_opts.png')

plt.show()
