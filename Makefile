all:
	g++ -g -O2 -W -Wall twoChoices_openQueue.cc -o simulate

clean:
	rm simulate *~

get_results_from_cluster_epfl:
	rsync -avb icsil1-access1.epfl.ch:velib/svn-velib/twoChoicesSimu/results .
