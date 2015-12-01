# MPI-Floyd-Warshall-C
 Parallel implementation (in C)  of the Floyd-Warshall algorithm using Fox algorithm in MPI to solve the "All-Pairs Shortest Paths" problem.

####Commands:
 mpirun -np 1 -hostfile mycluster program  < input12 > out12_4p_4m_1np 

####Cluster File Example:
localhost slots=2 <br/>
blabla@ssh.dcc.bla.bla@t0107 cpu=2 <br/>
^(Specify your machine)

####Input Example:
   6 <br />
   0 2 0 5 0 0 <br />
   0 0 0 0 0 0 <br />
   0 2 0 0 0 5 <br />
   0 0 0 0 1 0 <br />
   3 9 3 0 0 0 <br />
   0 0 0 0 1 0 <br />

####Bibliography:
 1.[Floyd-Warshall algorithm](http://math.mit.edu/~rothvoss/18.304.1PM/Presentations/1-Chandler-18.304lecture1.pdf) <br />
 2.[Fox algorithm](http://www.lac.inpe.br/~stephan/CAP-372/Fox_example.pdf)<br />
