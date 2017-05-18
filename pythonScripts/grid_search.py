import subprocess
import sys

# lowQ = [True, False]
# num_iterations = [10, 20, 30, 40, 50, 60, 70]
#score = [30, 45, 60, 75]
score = [65, 75, 50]
overlap = [16, 20]
bases = [200, 300]
trim = [30]


def main(path, cell):
	results = open(path + "the_results.txt", "w")
	i = 0
	for b in score:
	        for c in overlap:
	                for d in bases:
	                        for e in trim:
	                                print("---------------RUNNING ITERATION NUMBER----------------" + str(i))
	                                results.write("score " + str(b) +  " overlap " + str(c) + " bases " + str(d) + " trim " + str(e) + '\n')
					results.flush()
					print("score ", b, "overlap ", c, "bases ", d, "trim ", e)
                                        subprocess.call(["python", "trapes.py",
	                                "-path", path, "-bam", "sorted.bam",
	                                "-unmapped", "unmapped.bam", "-output", "output",
	                                "-sumF", path, "-genome", "mm10_ncbi",
	                                "-overlap", str(c), "-score", str(b), "-lowQ", "-bases", str(d), "-trim", str(e)], stdout=None)
          	                        results.write(subprocess.check_output(["cat", path + cell  + "/output.summary.txt"]))
	                                i += 1              
        results.close() 

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])
 

