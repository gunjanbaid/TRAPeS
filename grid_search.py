import subprocess

lowQ = [True, False]
# num_iterations = [10, 20, 30, 40, 50, 60, 70]
score = [15, 25, 30, 35, 45, 50, 55, 60, 70]
overlap = [5, 15, 20, 30, 40]
# bases = [-10, 20, 30, 40, 50, 60]

productive_runs = []

# with lowQ

i = 0
for b in score:
	for c in overlap:
			print(True, b, c)
			print("RUNNING ITERATION NUMBER " + str(i))
			subprocess.call(["python", "trapes.py",
			"-path", "single_end_cells/", "-bam", "sorted.bam",
			"-unmapped", "unmapped.bam", "-output", "output",
			"-sumF", "single_end_cells/" + str(i), "-genome", "hg38",
			"-overlap", str(c), "-score", str(b), "-lowQ"])
			i += 1

# without lowQ

for b in score:
	for c in overlap:
			print(False, b, c)
			print("RUNNING ITERATION NUMBER " + str(i))
			subprocess.call(["python", "trapes.py",
			"-path", "single_end_cells/", "-bam", "sorted.bam",
			"-unmapped", "unmapped.bam", "-output", "output",
			"-sumF", "single_end_cells/" + str(i), "-genome", "hg38",
			"-overlap", str(c), "-score", str(b), "-lowQ"])
			i += 1
