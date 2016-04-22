msgm
====

An implementation of the framework described in

"A Multiscale Variable-grouping Framework for MRF Energy Minimization",
by Omer Meir, Meirav Galun, Stav Yagev, Ronen Basri and Irad Yavneh

You are welcome to use to use, modify and redistribute the code,
but please cite our paper if you do.

Omer Meir, 2015 (contact: omer.meir@weizmann.ac.il)
========================================================================


USAGE
=====
Set parameters in msgmParams() (can leave default),
select an optimization method (can replace with your own!),
and pass a graphical model G (see description in msgm()) to msgm().

See 'msgmDemo()' for a demo script comparing the performance
of multiscale and single-scale inference on a square grid with
random energy potentials.


OPTIMIZATION MODULES
====================

The framework can be run without an optimization module by setting
(params.optimization = 'NONE'), or with an external optimization module.
Available is a wrapper for two optimization algorithms, which must be
downloaded and added to path. QPBO must be wrapped in order to call
it from Matlab.

	- QPBO, by Vladimir Kolmogorov
	  http://pub.ist.ac.at/~vnk/software/QPBO-v1.4.src.zip
	
	- LSA, by Lena Gorelick
	  http://www.csd.uwo.ca/~ygorelic/LSA_TR_v2.02.zip	
		
	  LSA needs to be patched to handle a degenerate case
	  that may occur when using it within this framework.
	  Replace Gorelick's 'reparamEnergy.m' with the version supplied herein.


ADDITIONAL NOTES
================

The code is structured such that the main components can easily be replaced,
e.g. the variable-grouping algorithm and the selection of an interpolation rule.
 
A single-scale optimization scheme can be wrapped and placed in msgmOptimizeScale().
 


	
	
	
	
	