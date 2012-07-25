## npar
8
## ######################################################################### ##


## ######################################################################### ##

## _________________________________________________________________________ ##
## ival		lb		ub		phz		prior		p1		p2		#Parameter   ##
   7.5		0		15		1		0			-10		20		#log_ddot_r
   7.5		0		15		1		0			-10		20		#log_bar_r
   0.20		0.01	0.50	4		2			-2.5257	0.50	#m_infty
   36.0		20.0	60.0	3		1			36.46	38.58	#l_infty
   0.18		0.01	0.90	3		1			0.243	0.23	#vbk
   0.65		0.00	2.00	-4		0			0.00	2.00	#beta
   10.0		3.00	20.0	2		2			2.30258	0.20	#mu_r
   0.20		0.00	1.00	2		3			20.0	80.0    #cv_r
## _________________________________________________________________________ ##


## nflags
5
## flags
0			# 1) Verbose (1==on) (2==more detail)
15.0		# 2) Minimum size (cm) of individual fish that can be tagged.
0.100		# 3) Std deviation in total catch in first phase.
0.050		# 4) Std deviation in total catch in last phase.
0			# 5) Case value of Size Transition (0=Estimate single P, 1=Estimate time varying linf, 2=use Size transitions)