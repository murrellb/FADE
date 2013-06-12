ExecuteAFile (PROMPT_FOR_FILE);
LikelihoodFunction lf 		= 	(filteredData, givenTree);

alpha := 1;
VERBOSITY_LEVEL				= 1;
AUTO_PARALLELIZE_OPTIMIZE	= 1;
Optimize 						(res0,lf);
AUTO_PARALLELIZE_OPTIMIZE	= 0;
VERBOSITY_LEVEL				= -1;
ClearConstraints(alpha);

/* export baseline model LF */
LIKELIHOOD_FUNCTION_OUTPUT = 7;
outPath = basePath + ".base";
fprintf (outPath, CLEAR_FILE, lf);
