# FLASH_NTCP
Code and data for "Effect of the dose distribution and organ response on the toxicity in FLASH radiotherapy: a modeling study"

VARIABLES IN THE DATA.MAT WORKSPACE

D_het: vector of simulated doses, heterogeneous case.

D_hom: vector of simulated doses, homogeneous case.

d_rate:  vector of simulated dose rates.

dose_R: heterogeneous dose distribution.

F_*: interpolants for the calculation of SF_FLASH as a function of dose and dose rate:

	hom/het refers to homogeneous or heterogeneous oxygenations.
	het_33 and het_90 refers to the two heterogeneous oxygenations studied in this work.
	“n”mmHg refers th the oxygenation “n” in the homogeneous oxygenation scenarios.
	ab10 refers to and alpha/beta of 10 Gy.	

FMF: flash modifying factor for the values dose and dose rate in D_hom and d_rate.

n: vector of volume response parameter values.

NTCP_ROD*: Nx3 matrices containing the results of the NTCP calculations using the ROD model, where N is the number of n values (volume parameter). The first column contains the NTCP values for conventional RT, the second column contains the NTCP values for FLASH, the third column contains the values of D_50 for each value of n:

	hom/het refers to homogeneous or heterogeneous oxygenations.
	het_33 and het_90 refers to the two heterogeneous oxygenations studied in this work.
	“n”mmHg refers th the oxygenation “n” in the homogeneous oxygenation scenarios.
	ab10 refers to and alpha/beta of 10 Gy.
	ref02 refers to the reference NTCP value for conventional delivery 0.2.
	D50_20 refers to the reference D_50=20Gy value for the calculation of the NTCP with the 	LKB model for n=1.

NTC_mod2: Nx3 matrix containing the results of the NTCP calculations using the phenomenological model, where N is the number of n values (volume parameter). The first column contains the NTCP values for conventional RT, the second column contains the NTCP values for FLASH, the third column contains the values of D_50 for each value of n. 

pp_conv*: interpolants for the calculation of SF_CONV as a function of dose:

	hom/het refers to homogeneous or heterogeneous oxygenations.
	het_33 and het_90 refers to the two heterogeneous oxygenations studied in this work.
	“n”mmHg refers th the oxygenation “n” in the homogeneous oxygenation scenarios.
	ab10 refers to and alpha/beta of 10 Gy.

SF_conv*: surviving fractions for conventional delivery for different dose values and (in the case of homogeneous oxygenation) oxygenation values:

	hom/het refers to homogeneous or heterogeneous oxygenations.
	het_33 and het_90 refers to the two heterogeneous oxygenations studied in this work.
	ab10 refers to and alpha/beta of 10 Gy.

SF_flash*: surviving fractions for FLASH delivery for different dose values, dose rate values, and (in the case of homogeneous oxygenation) oxygenation values:

	hom/het refers to homogeneous or heterogeneous oxygenations.
	het_33 and het_90 refers to the two heterogeneous oxygenations studied in this work.
	ab10 refers to and alpha/beta of 10 Gy.

u_hom:  vector of simulated oxygenations, homogeneous case.
