# 图1
- 实验配置与描述：当没有任何优化时，每次实验仅增加参与完成计算的线程数目，线程的数量为单线程、2线程、4线程、16线程、32线程、64线程。这些线程都是在同一节点上的。因此不存在不同节点中内存访问的问题。在实验期间，可能存在访问冲突的变量已使用reduction避免冲突。

- 实验结果为：单线程为基准，2线程加速比为1.01,4线程加速比为1.05,16线程加速比为1.18,32线程加速比为1.05,64线程加速比为1.19

- 分析：实验结果可得，当完成相同迭代数、相同鱼总数（2^22尾）的任务时，线程数数增加，加速比增加。当线程数为16、64时加速比分别达到1.18和1.19。通过观察可得，32线程的情况下获得的加速比低于16线程、64线程的加速比，这有可能由以下几点原因造成：
1. reduction 的引用导致一个线程更新共享变量时，其他线程需要等待更新的完成。
2. 非均匀的内存访问，当使用32线程时可能存在有些线程所在的核心与存放数据的内存的距离较远，产生了一定的访问延迟。
3. 更多的线程可能导致更频繁的缓存未命中，继而需要从内存中读取数据，导致了性能的降低。

# 图2
- 实验配置与描述：在这个实验中使用16线程完成任务，迭代500轮，每次实验中的鱼的总数不同，分别为2^6, 2^10, 2^14, 2^18, 2^22。由于鱼群总数不同，为了实验结果可以在同一尺度下对比，将完成任务所消耗的时间除以当前实验中鱼群的总数，得到完成一条鱼迭代500轮所需要的时间。

- 实验结果为：以池塘中存在2^6条鱼，每条鱼完成500轮迭代所需的时间为基准值。2^10条鱼时加速比为15.72， 2^14条鱼时加速比为562.20, 2^18条鱼时加速比为72017.19, 2^22条鱼时加速比为7237.91。

- 分析：通过观察实验结果，可以得到如下的现象：
1. 当问题规模较大时，计算任务更容易均匀分布在所有处理器或线程上，从而提高资源利用率。
2. 随着问题规模增加，加速比先显著增加，但在2^22条鱼的实验条件下，加速比有所降低。
3. 对于大规模问题，计算通常会占用更大比例的时间，而与其他线程或处理器进行数据交换和通信的时间相对会减少。
4. 由于没有使用栈对部分堆上的数据进行复制，导致了加速比不能进一步上升。并且由于数据量的上升，内存访问的次数多了，产生了瓶颈，继而导致了加速比下降。

# 图3
- 实验配置：实验中使用了五种不同的调度方法，包括"Guided"、“Dynamic, chunk=2^18”、“Static, chunk=2^18”、“Static, chunk=1"和"Basic”。
- 实验结果：实验结果以条形图的形式展示，x轴表示相对于基线（没有调度方法）的加速比，y轴表示每个调度方法的加速比。
- 分析：从结果来看，"Static, chunk=2^18"调度方法与基线相当，没有提升。其余的调度方法都导致了不同程度的性能损失。这是由于Fish数组声明在堆上导致了每次进行计算的时候都需要从堆中访问数组。

# 图4
- 实验配置：在实验中，我们使用了四种不同的编译器优化级别，分别是"Base"（无优化）、"O1"、"O2"和"O3"。每个测试用例均以16个线程和2^22条鱼进行运行。

- 实验结果：实验的结果以条形图形式展示。其中，x轴表示各测试用例相对于没有编译器优化的基准（Base）的加速比，y轴表示不同的测试用例。在所有测试用例中，使用"O2"和"O3"编译器优化级别的加速比最高，达到了约2.65倍，而未进行任何编译器优化的"Base"用例的加速比为1.00。

- 分析：实验结果表明，编译器优化能够显著提升程序性能。特别是使用了"O2"和"O3"优化级别后，程序的运行速度比没有任何优化的基准快了约165%。这一提升可能源于这些优化级别能更有效地利用硬件资源和减少运行时的计算开销。

# 图5
- 实验配置：在实验中，每个线程完成任务的时候需要将分配到的任务从堆复制到栈中，然后再进行计算。这个实验使用了16个线程，2^22条鱼。

- 实验结果：其中，x轴表示各测试用例相对于没有编译器优化的基准（Base）的加速比，y轴表示不同的测试用例。在所有测试用例中，其中，优化后的两个测试用例“仅使用栈优化”，“栈与O3优化”加速比分别达到了1.82和4.24。

- 分析：实验结果表明，使用栈保存当前需要运算的数据和编译器优化选项能够显著提升程序性能。

# Figure 1
- **Experimental Configuration and Description**: Without any optimization, each experiment increases the number of threads involved in the calculation. The number of threads are: single-thread, 2 threads, 4 threads, 16 threads, 32 threads, and 64 threads. All these threads are on the same node, so there are no issues related to memory access across different nodes. During the experiments, variables that may have access conflicts have been avoided using `reduction`.

- **Results**: Using a single thread as the baseline, the speedup ratio for 2 threads is 1.01, for 4 threads is 1.05, for 16 threads is 1.18, for 32 threads is 1.05, and for 64 threads is 1.19.

- **Analysis**: 
  1. Using `reduction` can cause one thread to wait for another when updating shared variables.
  2. Non-uniform memory access: when using 32 threads, some threads may have longer access times to the memory where data is stored.
  3. More threads may lead to more frequent cache misses, resulting in decreased performance.

# Figure 2
- **Experimental Configuration and Description**: In this experiment, the task is performed with 16 threads for 500 iterations. The total number of fish in each experiment varies: \(2^6\), \(2^{10}\), \(2^{14}\), \(2^{18}\), \(2^{22}\).

- **Results**: The baseline is the time required for each fish to complete 500 iterations with \(2^6\) fish. The speedup ratios are: \(2^{10}\) - 15.72, \(2^{14}\) - 562.20, \(2^{18}\) - 72017.19, \(2^{22}\) - 7237.91.

- **Analysis**:
  1. Larger problem sizes result in more evenly distributed computation tasks across all processors or threads, thus improving resource utilization.
  2. The speedup ratio first increases significantly with the size of the problem, but then decreases for \(2^{22}\) fish.
  3. For larger problem sizes, computation generally takes up a larger proportion of the time, while the time for data exchange and communication with other threads or processors is relatively reduced.
  4. Speedup ratio declines due to increased memory accesses and not copying part of the heap data to the stack.

# Figure 3
- **Experimental Configuration**: Five different scheduling methods were used, including "Guided", "Dynamic, chunk=\(2^{18}\)", "Static, chunk=\(2^{18}\)", "Static, chunk=1", and "Basic".

- **Results**: The results are shown in a bar graph. The x-axis represents the speedup ratio compared to the baseline (no scheduling method), and the y-axis represents the speedup ratio for each scheduling method.

- **Analysis**: The "Static, chunk=\(2^{18}\)" scheduling method performs comparably to the baseline. All other methods resulted in some performance loss due to heap access for the Fish array.

# Figure 4
- **Experimental Configuration**: Four different compiler optimization levels were used: "Base" (no optimization), "O1", "O2", and "O3". Each test case is run with 16 threads and \(2^{22}\) fish.

- **Results**: The results are displayed in the form of a bar graph. The x-axis represents the speedup ratio for each test case relative to the baseline, which has no compiler optimization. The highest speedup ratios are achieved with "O2" and "O3", approximately 2.65 times faster, while the "Base" case has a speedup ratio of 1.00.

- **Analysis**: Compiler optimizations can significantly improve program performance. Specifically, using "O2" and "O3" optimization levels, the program runs approximately 165% faster than the unoptimized baseline. This increase is likely due to these optimization levels being more effective at utilizing hardware resources and reducing runtime computational overhead.

# Figure 5
- **Experimental Configuration**: In this experiment, each thread is required to copy its allocated task from the heap to the stack before performing the computation. The experiment uses 16 threads and \(2^{22}\) fish. According to the cache architecture of the Setonix supercomputer, the optimal number of fish per thread is determined to be 1.5M.

- **Results**: The x-axis represents the speedup ratio for each test case relative to the baseline, which has no compiler optimization. The y-axis represents different test cases. Among all test cases, the optimized test cases "Stack-only Optimization" and "Stack with O3 Optimization" achieved speedup ratios of 1.82 and 4.24, respectively.

- **Analysis**: The results indicate that using the stack to store data currently required for computation, along with compiler optimization options, can significantly improve program performance.
