pThresh = 0.001;

function AddABiasFADE (ModelMatrixName&, ModelMatrixName2&, biasedBase)
{
	ModelMatrixName2 = {20,20};
	
	t = 1;	/* branch length, local parameter */
	_numericRateMatrix = ModelMatrixName;
	
	/* the probability that site is undergoing biased substitution rates */
	/*
	global	  P_bias = 0.1;  P_bias :< 0.5;
	category catVar = (2,{{1-P_bias,P_bias}},MEAN,,{{0,1}},0,1);
	*/
	global alpha = 1;
	global beta = 1;
	
	for (ri = 0; ri < 20; ri = ri+1)
	{
		for (ci = ri+1; ci < 20; ci = ci+1)
		{
			ModelMatrixName2[ri][ci] := _numericRateMatrix__[ri__][ci__] * t * alpha;
			ModelMatrixName2[ci][ri] := _numericRateMatrix__[ri__][ci__] * t * alpha;
		
		}
	}
	
	if (biasedBase < 20)
	{
		global rateBiasTo 	 := beta;
		global rateBiasFrom	 := 1/rateBiasTo;
			
		rateBiasTo    :>1;
		relBias       :>1;	/* UNUSED ?!? */
		for (ri = 0; ri < 20; ri = ri+1)
		{
			if (ri != biasedBase)
			{
				ModelMatrixName2[ri][biasedBase] := _numericRateMatrix__[ri__][biasedBase__] * t * alpha * rateBiasTo;
				ModelMatrixName2[biasedBase][ri] := _numericRateMatrix__[ri__][biasedBase__] * t * alpha * rateBiasFrom;
			} 
		}
	}

	return 1;
}

function AddABiasREL (ModelMatrixName&, ModelMatrixName2&, biasedBase)
{
	ModelMatrixName2 = {20,20};
	
	t = 1;	/* branch length, local parameter */
	c = 1;	/* rate variation */
	_numericRateMatrix = ModelMatrixName;
	
	/* the probability that site is undergoing biased substitution rates */
	global	  P_bias = 0.1;  P_bias :< 0.5;
	
	
	category catVar = (2,{{1-P_bias,P_bias}},MEAN,,{{0,1}},0,1);
	
	for (ri = 0; ri < 20; ri = ri+1)
	{
		for (ci = ri+1; ci < 20; ci = ci+1)
		{
			ModelMatrixName2[ri][ci] := _numericRateMatrix__[ri__][ci__] * t * c;
			ModelMatrixName2[ci][ri] := _numericRateMatrix__[ri__][ci__] * t * c;
		
		}
	}

	if (biasedBase < 20)
	{
		global rateBiasTo 	  = 5.0;
		global rateBiasFrom	 := 1/rateBiasTo;
			
		rateBiasTo    :>1;
		relBias       :>1;	/* UNUSED ?!? */
		for (ri = 0; ri < 20; ri = ri+1)
		{
			if (ri != biasedBase)
			{
				ModelMatrixName2[ri][biasedBase] := _numericRateMatrix__[ri__][biasedBase__] * t * c * ((catVar==1)*rateBiasTo+(catVar==0));
				ModelMatrixName2[biasedBase][ri] := _numericRateMatrix__[ri__][biasedBase__] * t * c * ((catVar==1)*rateBiasFrom+(catVar==0));
			}
		}
	}

	return 1;
}

function defineFadeGrid (gammaPoints, biasPoints) { // gamma = site-to-site rate variation
	
	fprintf(stdout,"defineFadeGrid","\n");
	
	fraction_gamma_points_less_than_one = 0.5;
	fraction_bias_points_less_than_one = 0;
	
	minGammaValue = 0.1;
	maxGammaValue = 50;
	maxBiasValue = 50;
	
	gammaAndBiasGrid = {gammaPoints * biasPoints,2}; 
	
	gammaPoints = Max(gammaPoints, 10);
	biasPoints = Max(biasPoints, 10);
	
	gamma_points_less_than_one = gammaPoints * fraction_gamma_points_less_than_one $ 1;
	bias_points_less_than_one = biasPoints * fraction_bias_points_less_than_one $ 1;
	
	maxGammaStep = (maxGammaValue-1)^(1/3)/(gammaPoints-gamma_points_less_than_one-1);
	maxBiasStep = (maxBiasValue-1)^(1/3)/(biasPoints-bias_points_less_than_one-1);
	
	index = 0;
	for(i = 0 ; i < gammaPoints ; i += 1)
	{
		for(j = 0 ; j < biasPoints ; j += 1)
		{
			if(i <= gamma_points_less_than_one)
			{
				gammaAndBiasGrid[index][0] = minGammaValue + (1-minGammaValue)*(i / gamma_points_less_than_one);
				if(gamma_points_less_than_one == 0)
				{
					gammaAndBiasGrid[index][0] = 1;
				}
			}
			else
			{
				gammaAndBiasGrid[index][0] = 1+(maxGammaStep*(i - gamma_points_less_than_one))^3;
			}
			
			if(j <= bias_points_less_than_one)
			{
				gammaAndBiasGrid[index][1]  = j / bias_points_less_than_one;
				if(bias_points_less_than_one == 0)
				{
					gammaAndBiasGrid[index][1] = 1;
				}
			}
			else
			{
				gammaAndBiasGrid[index][1] = 1+(maxBiasStep*(j - bias_points_less_than_one))^3;
			}
			
			index +=1;
		}
	}


	index = 0;
	for(i = 0 ; i < gammaPoints ; i += 1)
	{
		for(j = 0 ; j < biasPoints ; j += 1)
		{
			//gammaAndBiasGrid[index][0] = 1;
			//gammaAndBiasGrid[index][1] = 1;
			index +=1;
		}
	}

	/*
    alphaBetaGrid = {one_d_points^2,2}; // (alpha, beta) pair
    oneDGrid      = {one_d_points,1};
   
    one_d_points    = Max (one_d_points, 10);
    neg_sel         = 0.7;
    neg_sel_points  = ((one_d_points)*neg_sel+0.5)$1;
    pos_sel_points  = (one_d_points-1)*(1-neg_sel)$1;
    if (neg_sel_points + pos_sel_points != one_d_points) {
        pos_sel_points = one_d_points - neg_sel_points; 
    }
    _neg_step = 1/neg_sel_points;
    for (_k = 0; _k < neg_sel_points; _k += 1) {
        oneDGrid [_k][0] =  _neg_step * _k;
    }
    oneDGrid [neg_sel_points-1][0] = 1;
    _pos_step = 49^(1/3)/pos_sel_points;
    for (_k = 1; _k <= pos_sel_points; _k += 1) {
        oneDGrid [neg_sel_points+_k-1][0] = 1+(_pos_step*_k)^3;
    }
    
    _p = 0;
    for (_r = 0; _r < one_d_points; _r += 1) {
        for (_c = 0; _c < one_d_points; _c += 1) {
           alphaBetaGrid[_p][0] = oneDGrid[_r];
           alphaBetaGrid[_p][1] = oneDGrid[_c];
           _p += 1;
        }
    }*/
	
	fprintf(stdout,"END","\n");
    
    return gammaAndBiasGrid;   
}


/*--------------------------------------------------------------------------------------------*/

function AddABiasREL (ModelMatrixName&, ModelMatrixName2&, biasedBase)
{
	ModelMatrixName2 = {20,20};
	
	t = 1;	/* branch length, local parameter */
	c = 1;	/* rate variation */
	_numericRateMatrix = ModelMatrixName;
	
	/* the probability that site is undergoing biased substitution rates */
	global	  P_bias = 0.1;  P_bias :< 0.5;
	
	
	category catVar = (2,{{1-P_bias,P_bias}},MEAN,,{{0,1}},0,1);
	
	for (ri = 0; ri < 20; ri = ri+1)
	{
		for (ci = ri+1; ci < 20; ci = ci+1)
		{
			ModelMatrixName2[ri][ci] := _numericRateMatrix__[ri__][ci__] * t * c;
			ModelMatrixName2[ci][ri] := _numericRateMatrix__[ri__][ci__] * t * c;
		
		}
	}

	if (biasedBase < 20)
	{
		global rateBiasTo 	  = 5.0;
		global rateBiasFrom	 := 1/rateBiasTo;
			
		rateBiasTo    :>1;
		relBias       :>1;	/* UNUSED ?!? */
		for (ri = 0; ri < 20; ri = ri+1)
		{
			if (ri != biasedBase)
			{
				ModelMatrixName2[ri][biasedBase] := _numericRateMatrix__[ri__][biasedBase__] * t * c * ((catVar==1)*rateBiasTo+(catVar==0));
				ModelMatrixName2[biasedBase][ri] := _numericRateMatrix__[ri__][biasedBase__] * t * c * ((catVar==1)*rateBiasFrom+(catVar==0));
			}
		}
	}

	return 1;
}

/*--------------------------------------------------------------------------------------------*/

function promptModel (dummy)
{
	ChoiceList	     (pickAModel,"Subsitution Model",1,SKIP_NONE, "HIV Within","HIV Within",
							 								      "HIV Between","HIV Between",
							 								      "JTT","JTT",
							 								      "Flu H5N1", "Empirical model for H5N1 Influenza",
																  "LG", "Le-Gasquel 2008",
							 								      "REV", "Use general time-reversible model (WARNING: 189 rate parameters will be estimated from your alignment)");
							 								     
							 								    
	if (pickAModel < 0)
	{
		return 0;
	}

	modelSTDINoverload = {};
	modelSTDINoverload["2"] = "Rate variation";
	modelSTDINoverload["3"] = "Beta-Gamma";
	modelSTDINoverload["4"] = "4";
	modelNameString = "_customAAModelMatrix";

	if (pickAModel == 0)
	{
		modelSTDINoverload["0"] = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TemplateModels" + DIRECTORY_SEPARATOR + "EmpiricalAA" + 
								  DIRECTORY_SEPARATOR + "HIVWithin";
	}

	if (pickAModel == 1)
	{
		modelSTDINoverload["0"] = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TemplateModels" + DIRECTORY_SEPARATOR + "EmpiricalAA" + 
								  DIRECTORY_SEPARATOR + "HIVBetween";
	}

	if (pickAModel == 2)
	{
		modelSTDINoverload["0"] = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TemplateModels" + DIRECTORY_SEPARATOR + "EmpiricalAA" + 
								  DIRECTORY_SEPARATOR + "JTT";
	}

	if (pickAModel == 3)
	{
		modelSTDINoverload["0"] = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TemplateModels" + DIRECTORY_SEPARATOR + "EmpiricalAA" + 
								  DIRECTORY_SEPARATOR + "H5N1";
	}
	
	if (pickAModel == 4)
	{
		modelSTDINoverload["0"] = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TemplateModels" + DIRECTORY_SEPARATOR + "EmpiricalAA" + 
								  DIRECTORY_SEPARATOR + "LG";
	}

	if (pickAModel < 5)
	{
		/* estimate base frequencies as model parameters */
		modelSTDINoverload["1"]	= "ML Estimates";
		modelPath = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TemplateModels" + DIRECTORY_SEPARATOR + "Custom_AA_empirical.mdl";
	}
	else
	{
		modelNameString = "mtREVMatrix";
		modelPath = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TemplateModels" + DIRECTORY_SEPARATOR + "mtREV.mdl";	
	}
	
	ExecuteAFile (modelPath,modelSTDINoverload);

	if (dummy)
	{
		ChoiceList	     (pickATarget,"Target residue",1,SKIP_NONE, "Fixed","Fixed Sequence",
							 								    "Inferred","Inferred Sequence");
	}						 								    
	return pickAModel;
}