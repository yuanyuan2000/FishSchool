当没有任何优化时，每次实验仅增加参与完成计算的线程数目，线程的数量为单线程、2线程、4线程、16线程、32线程、64线程。这些线程都是在同一节点上的。因此不存在不同节点中内存访问的问题。在实验期间，可能存在访问冲突的变量已使用reduction

实验结果可得，当完成相同迭代数、相同鱼总数（2^22尾）的任务时，线程数数增加，加速比增加。当线程数为16、64时加速比分别达到1.18和1.19。通过观察可得，32线程的情况下获得的加速比低于16线程、64线程的加速比，这有可能由以下几点原因造成：
1. 