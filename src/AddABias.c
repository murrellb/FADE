pThresh = 0.001;

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

/*
 * Adds a bias to the specified amino acid in the FADE model.
 * @param biasedBase a numeric value (0 - 19) specifying the amino acid to be biased.
 */
function AddABiasFADE (ModelMatrixName&, ModelMatrixName2&, biasedBase)
{
	ModelMatrixName2 = {20,20};
	
	gamma = 1;
	bias = 1;
	_numericRateMatrix = ModelMatrixName;
	
	for (ri = 0; ri < 20; ri = ri+1)
	{
		for (ci = ri+1; ci < 20; ci = ci+1)
		{
			ModelMatrixName2[ri][ci] := _numericRateMatrix__[ri__][ci__] * gamma;
			ModelMatrixName2[ci][ri] := _numericRateMatrix__[ri__][ci__] * gamma;
		
		}
	}

	if (biasedBase < 20)
	{
		global rateBiasTo 	  = bias;
		global rateBiasFrom	 := 1/bias;
			
		for (ri = 0; ri < 20; ri = ri+1)
		{
			if (ri != biasedBase)
			{
				ModelMatrixName2[ri][biasedBase] := _numericRateMatrix__[ri__][biasedBase__] * gamma * rateBiasTo;
				ModelMatrixName2[biasedBase][ri] := _numericRateMatrix__[ri__][biasedBase__] * gamma * rateBiasFrom;
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