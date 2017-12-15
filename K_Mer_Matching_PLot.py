import matplotlib.pyplot as plt
import numpy as np

# M. tuberculosis
test1 = [82, 69, 3, 2, 2]
leg1 = 'k = 17, min. 2 kmer matches, ref. genomes restricted to top 0.5 SD # matches'
test2 = [48, 2, 1, 1, 1]
leg2 = 'k = 31, min. 1 kmer matches, ref. genomes restricted to top 0.5 SD # matches'
test4 = [100, 100, 0, 0, 0]
leg4 = 'k = 11, min. 2 kmer matches, ref. genomes restricted to top 0.5 SD # matches [BLAT]'
test5 = [75, 48, 1, 1, 1]
leg5 = 'k = 16, min. 4 kmer matches, ref. genomes restricted to top 0.5 SD # matches'
test6 = [56, 24, 9, 6, 6]
leg6 = 'k = 21, min. 1 kmer matches, all ref. genomes considered'
labels = ['Precision', 'Domain \n Specificity', 'Order \n Specificity', 'Family \n Specificity', 'Genus \n Specificity']

plt.scatter(np.arange(0, 5), test1, label=leg1)
plt.scatter(np.arange(0, 5), test2, label=leg2)
plt.scatter(np.arange(0, 5), test5, label=leg5)
plt.scatter(np.arange(0, 5), test6, label=leg6)
plt.xticks(np.arange(0, 5), labels)
plt.ylabel('Percent of Total Reads')
plt.subplots_adjust(left=None, bottom=0.35, right=None, top=None, wspace=None, hspace=None)
plt.legend(bbox_to_anchor=(-.15, -0.65), loc="lower left")
plt.show()
