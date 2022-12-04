# mc
This Git- repository include sourse to run experiments  for the improvement of the algorithm in the paper 
X.Wang,Improving the rejection sampling method in quasi-Monte Carlo methods,J.Comput.Appl.Math.114(2000),
231–246.

lattice rules
It contains three ways to generated lattice rules，the basic component by component（cbc）way to generate the generator，
improved fast cbc algrithom（proposed by Dirk Nuyens), and the pod weight lattice rule

importance sampling
It contains an example that importance sampling and lattice rule are combined to calculate a Randleman-Bartter model.

AR
It contains algorithms which use different qmc point set in reject sampling to calculate an integral.We use three smooth
methods: linear smooth(proposed by Wang,2000),cubic curve smooth and conditional smooth methods.
