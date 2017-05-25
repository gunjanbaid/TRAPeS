import subprocess
import sys

# lowQ = [True, False]
# num_iterations = [10, 20, 30, 40, 50, 60, 70]
#score = [30, 45, 60, 75]
score = [20, 25, 30, 35, 40, 45]
overlap = [20]
bases = [-10]
trim = [30]


def main(path, cell):
	i = 0
	for b in score:
	        for c in overlap:
	                for d in bases:
	                        for e in trim:
	                                print("---------------RUNNING ITERATION NUMBER----------------" + str(i))
					print("score ", b, "overlap ", c, "bases ", d, "trim ", e)
                                        subprocess.call(["python", "trapes.py",
	                                "-path", path, "-bam", "sorted.bam",
	                                "-unmapped", "unmapped.bam", "-output", "output" + str(i),
	                                "-sumF", path + "/" + str(i), "-genome", "mm10_ncbi",
	                                "-overlap", str(c), "-score", str(b), "-lowQ", "-bases", str(d), "-trim", str(e)])
	                                i += 1              
        results.close() 

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])
 

