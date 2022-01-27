
import swiftstat7 as swift
import numpy as np

print (swift.version())

if swift.mpislave():
	print "I was an MPI slave and will exit now."
	swift.mpifinalize()
	exit(0)

path = "/Users/max/OneDrive/Uni/RScripts/B03-PROENS/DATA"
corpus = "b03_proens"
seqid = "all"
seed = None
#"outputPath", "parmPath", "corpusFile", "fixseqName", "seed"
swift.mpiload(outputPath = path, parmPath = "%s/swpar_%s_%s.par" % (path, corpus, seqid), corpusFile = "%s/corpus_%s.dat" % (path, corpus), fixseqName = "%s/fixseqin_%s.dat" % (path, seqid), neval = 1, seed = np.random.randint(0,2**48) if seed is None else seed)

print "Begin likelihood evaluations..."
print swift.mpieval()
swift.mpiupdate("delta0", 10.0)
print swift.mpieval()
swift.mpiupdate("delta0", 15.0)
print swift.mpieval()
#print swift.mpievallb()

print "I was the MPI master and will exit now."

swift.mpifinalize()

