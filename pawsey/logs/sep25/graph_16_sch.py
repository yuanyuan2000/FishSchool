import matplotlib.pyplot as plt

data = [83.9135995, 84.1932735, 83.894083,
        88.3546285, 88.02525625]

base = data[0]
speedup = [base / x for x in data]

labels = ['Base', 'Static, chunk:1', 'Static, chunk:2^18', 'Dynamic, chunk:2^18', 'Guided']
plt.barh(labels, speedup, color='skyblue')

plt.xlim(0.94, max(speedup)+0.01)

plt.xlabel('Speedup')
plt.ylabel('Test Case (run with 16 Threads, 2^22 Fishes)')
plt.title('Speedup Compared to Baseline (Without Schedule Method)')

for i, v in enumerate(speedup):
    plt.text(v, i, " {:.5f}".format(v), va='center',
             color='black', fontweight='bold')

plt.savefig('speedup_chart_diff_sch.png')

plt.show()
