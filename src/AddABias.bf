pThresh = 0.001;

/*
 * @param ModelMatrixName a string specifying the matrix you want to copy and and add a bias to.
 */
function AddABiasFADE (ModelMatrixName&, ModelMatrixName2&, biasedBase)
{
	ModelMatrixName2 = {20,20};
	
	t = 1;	/* branch length, local parameter */
	_numericRateMatrix = ModelMatrixName;

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

/*
 * @param ModelMatrixName points to the matrix variable you want to copy and and add a bias to.
 */
function AddABiasFADE2 (ModelMatrixName, ModelMatrixName2&, biasedBase)
{
	ModelMatrixName2 = {20,20};
	
	t = 1;	/* branch length, local parameter */
	_numericRateMatrix = ModelMatrixName;

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

function defineFadeGrid (alphaPoints, biasPoints) { // alpha = site-to-site rate variation

	fraction_alpha_points_less_than_one = 0.5;
	fraction_bias_points_less_than_one = 0;
	
	minAlphaValue = 0.1;
	maxAlphaValue = 50;
	maxBiasValue = 50;
	
	alphaAndBiasGrid = {alphaPoints * biasPoints,2}; 
	
	alphaPoints = Max(alphaPoints, 10);
	biasPoints = Max(biasPoints, 10);
	
	alpha_points_less_than_one = alphaPoints * fraction_alpha_points_less_than_one $ 1;
	bias_points_less_than_one = biasPoints * fraction_bias_points_less_than_one $ 1;
	
	maxGammaStep = (maxAlphaValue-1)^(1/3)/(alphaPoints-alpha_points_less_than_one-1);
	maxBiasStep = (maxBiasValue-1)^(1/3)/(biasPoints-bias_points_less_than_one-1);
	
	index = 0;
	for(i = 0 ; i < alphaPoints ; i += 1)
	{
		for(j = 0 ; j < biasPoints ; j += 1)
		{
			if(i <= alpha_points_less_than_one)
			{
				alphaAndBiasGrid[index][0] = minAlphaValue + (1-minAlphaValue)*(i / alpha_points_less_than_one);
				if(alpha_points_less_than_one == 0)
				{
					alphaAndBiasGrid[index][0] = 1;
				}
			}
			else
			{
				alphaAndBiasGrid[index][0] = 1+(maxGammaStep*(i - alpha_points_less_than_one))^3;
			}
			
			if(j <= bias_points_less_than_one)
			{
				alphaAndBiasGrid[index][1]  = j / bias_points_less_than_one;
				if(bias_points_less_than_one == 0)
				{
					alphaAndBiasGrid[index][1] = 1;
				}
			}
			else
			{
				alphaAndBiasGrid[index][1] = 1+(maxBiasStep*(j - bias_points_less_than_one))^3;

				//make 50% of betas = 1 (so prior is not so unfair)				
				/*if(j < 10)
				{
					alphaAndBiasGrid[index][1] = 1;
				}
				else
				{
					alphaAndBiasGrid[index][1] = 1+((j-10)/50) + (((j-10)^3)/14.87);
				}*/
				
			}
			
			index +=1;
		}
	}


	index = 0;
	for(i = 0 ; i < alphaPoints ; i += 1)
	{
		for(j = 0 ; j < biasPoints ; j += 1)
		{
			//alphaAndBiasGrid[index][0] = 1;
			//alphaAndBiasGrid[index][1] = 1;
			index +=1;
		}
	}
    return alphaAndBiasGrid;   
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
	modelSTDINoverload["2"] = "Fixed Rates";
	//modelSTDINoverload["3"] = "Beta-Gamma";
	//modelSTDINoverload["4"] = "4";
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
