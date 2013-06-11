/* define an associative array with key = amino acid and value = integer
	index for AA in alphabetical order */
AAString    = "ACDEFGHIKLMNPQRSTVWY";
             //01234567890123456789


//run_residue = 16; // run only this residue, -1 => run all residues
run_residue = -1;
num_alpha = 15; // number of site-to-site rate variation grid points
num_beta = 15; // number of bias grid points
_cachingOK = 1;
concentration = 0.5;
save_conditionals = 0;  // for testing purposes only

if(run_residue < 0)
{
	fprintf(stdout, "Testing all amino acids...", "\n");
}
else
{
	fprintf(stdout, "Testing only amino acid ",AAString[run_residue],"...", "\n");
}

SKIP_MODEL_PARAMETER_LIST = 0;

#include "FADE_tools.ibf";
#include "GrabBag.bf";
#include "FUBAR_tools_iterative.ibf";
#include "CodonToProtein.bf";
LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("ReadDelimitedFiles");

ChoiceList						(reloadFlag, "Reload/New", 1, SKIP_NONE, "New analysis","Start a new analysis",
																	      "Reload","Reload a baseline protein fit.");
																		  
ACCEPT_ROOTED_TREES 			= 1;

//PHASE 1 NEEDS TO BE MOVED TO A SEPERATE FILE, SO THAT beomap CAN BE USED SO THAT THE TREE DOESN'T OPTIMIZE ON THE MASTER NODE
// Phase 1: Optimize the baseline model or load existing from file
if (reloadFlag == 0) // optimize baseline model
{
//	ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{HYPHY_LIB_DIRECTORY[0][Abs(HYPHY_LIB_DIRECTORY)-2],"TemplateBatchFiles","_MFReader_.ibf"}}));

	/* new analysis, fit baseline model */
	
	DataSet	ds					 = ReadDataFile (PROMPT_FOR_FILE); // must be protein or codon alignment
	basePath 					 = LAST_FILE_PATH;
	DataSetFilter   filteredData = CreateFilter (ds,1);
	GetDataInfo		(checkCharacters, filteredData, "CHARACTERS");
	if (Columns (checkCharacters) != 20) // check if amino acid
	{
		if(Columns (checkCharacters) < 10) // convert to amino acid if codon
		{
			fprintf(stdout, "Codon data needs to be translated to amino acid data.");
			translatedDataset  = translateCodonToAminoAcid(LAST_FILE_PATH, 1);
			DataSetFilter filteredData = CreateFilter (bigDataSet,1);
		}
		else
		{			
			fprintf (stdout, "ERROR: please ensure that the input alignment contains codon or protein sequences");
			return 0;
		}
	}
	
	
	promptModel (0); // prompts user for AA model - no rate variation

	AddABiasFADE					(modelNameString,"backgroundMatrix",21); // 21 = don't add a bias
	Model						FG = (backgroundMatrix, vectorOfFrequencies, 1); 
	Model						BG = (backgroundMatrix, vectorOfFrequencies, 1);
	
	
	// prompts the user for a tree, returns givenTree
	ExecuteAFile 					(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "queryTree.bf");
	
	//DANGER POINT: NEED TO FIGURE OUT WHAT IS GOING ON WITH THE ROOTING. DO BOTH CHOICES IN THIS LIST BEHAVE PROPERLY?
	//DOES THE TREE EVEN GET REROOTED? WHAT SHOULD WE DO ON DATAMONKEY? SAME AS CURRENT DEPS: ROOT ON SELECTED NODE?
	ChoiceList						(fixFB, "Fix Branch", 1, SKIP_NONE, "Unknown root","The character at the root of the tree is drawn from the stationary distribution",
																		"Fix 1st sequence as root","The 1st sequence in the file (assumed to be a direct descendant of the root) is used to populate the root sequences.");


	/* check if the tree is rooted */
	
	treeAVL  = givenTree^0;
	rootNode = treeAVL[(treeAVL[0])["Root"]];
	if (Abs(rootNode["Children"]) != 2)
	{
		fprintf (stdout, "ERROR: please ensure that the tree is rooted");
		return 0;
	}
	root_left  = "givenTree." + (treeAVL[(rootNode["Children"])[0]])["Name"] + ".t";
	root_right = "givenTree." + (treeAVL[(rootNode["Children"])[1]])["Name"] + ".t";
	
	if (fixFB>0)
	{
		/* branch to first sequence is collapsed to length zero;
			enforcing identifiability with root sequence */
			
		ExecuteCommands ("givenTree."+TipName(givenTree,0)+".t:=0");
	}
	else
	{
		if (fixFB < 0)
		{
			return 0;
		}
	}

	//DANGER POINT: If this runs and the "fixFB>0" block above runs, then won't the tree will have two 0 length branches coming off the root?
	/* only the sum of branch lengths leading to two immediate descendants
		of root can be estimated; this measure prevents one of the branches
		from collapsing to zero length */
	ExecuteCommands					(root_left + ":=" + root_right); 
	
	LikelihoodFunction lf 		= 	(filteredData, givenTree);
	fprintf							(stdout, "[PHASE 1.1] Standard model fit\n"); 
	
	
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
	
	baselineLogL				= res0[1][0];
}
else
{
	/* import baseline LF from file */
	modelNameString = "_customAAModelMatrix";
	SetDialogPrompt ("Locate an existing fit:");
	ExecuteAFile (PROMPT_FOR_FILE);
	GetString (lfInfo,lf,-1);
	if ((lfInfo["Models"])[0] == "mtREVModel")
	{
		modelNameString = "mtREVMatrix";	
	}
	bpSplit						 = splitFilePath (LAST_FILE_PATH);
	basePath					 = bpSplit["DIRECTORY"] + bpSplit["FILENAME"];
	outPath						 = basePath + ".base";
	treeString 					 = Format (givenTree,0,0);	/* topology only, no branch lengths */

	treeAVL  = givenTree^0;
	rootNode = treeAVL[(treeAVL[0])["Root"]];
	if (Abs(rootNode["Children"]) != 2)
	{
		fprintf (stdout, "ERROR: please ensure that the tree is rooted");
		return 0;
	}


	ChoiceList						(fixFB, "Fix Branch", 1, SKIP_NONE, "Unknown root","The character at the root of the tree is drawn from the stationary distribution",
																		"Fix 1st sequence as root","The 1st sequence in the file (assumed to be a direct descendant of the root) is used to populate the root sequences.");
	if (fixFB>0)
	{
		ExecuteCommands ("givenTree."+TipName(givenTree,0)+".t:=0");
	}
	else
	{
		if (fixFB < 0)
		{
			return 0;
		}
	}
	
	LFCompute (lf,LF_START_COMPUTE);
	LFCompute (lf,baselineLogL);
	LFCompute (lf,LF_DONE_COMPUTE);
}

treeContainsTags = (treeString||"{FG}")[0][0] != -1;
//WHY DOES THIS GET CALLED AGAIN?
root_left  = "biasedTree." + (treeAVL[(rootNode["Children"])[0]])["Name"] + ".t";
root_right = "biasedTree." + (treeAVL[(rootNode["Children"])[1]])["Name"] + ".t";

/* vector of branch lengths for entire tree */
baselineBL						= BranchLength (givenTree,-1);

referenceL						= (baselineBL * (Transpose(baselineBL)["1"]))[0];

fprintf							(stdout,      "[PHASE 1.2] Standard model fit. Log-L = ",baselineLogL,". Tree length = ",referenceL, " subs/site \n"); 

ExecuteAFile 					(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "GrabBag.bf");
//DO WE NEED THIS?
ExecuteAFile 					(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "AncestralMapper.bf");

fixGlobalParameters ("lf");

byResidueSummary = {};
bySiteSummary	 = {};

/*------------------------------------------------------------------------------*/
//THIS IS ACTUALLY A CONCENTRATION MULTIPLIER NOW, SINCE THE PRIOR IS A VECTOR
_fubarPriorShape = prompt_for_a_value ("The concentration parameter of the Dirichlet prior",0.5,0.001,1,0);    
fprintf (stdout, "[DIAGNOSTIC] FADE will use the Dirichlet prior concentration parameter of ", _fubarPriorShape, "\n"); 

fadeGrid = 					defineFadeGrid (num_alpha, num_beta);

fprintf	(stdout,      "[PHASE 2] Computing conditional likelihoods.\n"); 
for (residue = 0; residue < 20; residue = residue + 1) // Stage 2, 3 and 4 fror each amino acid
{
	if(run_residue == residue || run_residue < 0)
	{ 
		AddABiasFADE2					(backgroundMatrix,"backgroundMatrix2",21); // 21 = don't add a bias
		AddABiasFADE2					(backgroundMatrix,"biasedMatrix",residue);	

		index = 0;
		for(_x = 0 ; _x < num_alpha ; _x += 1)
		{
			for(_y = 0 ; _y < num_beta ; _y += 1)
			{
				//fprintf(stdout,_x,",",_y,"\t",fadeGrid[index][0],"\t",fadeGrid[index][1],"\n");
				index += 1;
			}
		}


		Model				FG = (biasedMatrix, vectorOfFrequencies, 1); // vectorOfFrequencies comes from Custom_AA_empirical.mdl, in turn imported from a file such as "HIVWithin" rate matrix is multiplied by this vector (third argument)				
		if(treeContainsTags) // only run the background model if the tree contains tags
		{
			Model 			BG =  (backgroundMatrix2, vectorOfFrequencies, 1); // need to change this to baseline matrix
		}

		Tree				biasedTree = treeString;
		global				treeScaler = 1;

		/* constrains tree to be congruent to that estimated under baseline model */
		ReplicateConstraint 			("this1.?.?:=treeScaler*this2.?.?__",biasedTree,givenTree);
		//IS THIS SHIT REQUIRED IF WE AREN'T CALLING Optimize()?
		ExecuteCommands					(root_left + "=" + root_left);
		ExecuteCommands					(root_right + "=" + root_right);
		LikelihoodFunction lfb 		= 	(filteredData, biasedTree);
		
		gridInfoFile = basePath +"."+AAString[residue]+".grid_info";

		if(_cachingOK && !gridInfoFile)
		{  
			 fprintf (stdout, "[CACHED] Grid info file found for residue: ",AAString[residue], "\n"); 
		}
		else
		{
			gridInfo = computeLFOnGrid("lfb", fadeGrid, 1);
			points = Rows (gridInfo["conditionals"]);
			sites = Columns (gridInfo["conditionals"]);
			
			if(save_conditionals) // conditionals folder needs to be present in data directory, otherwise this will crash
			{
				conditionalsGrid = gridInfo["conditionals"];
				for(_s = 0 ; _s < sites ; _s += 1)
				{
					conditionals_out = basePath+"/conditionals/"+AAString[residue]+"."+_s+".csv";
					fprintf(conditionals_out,CLEAR_FILE,getGridStringForSite(conditionalsGrid,fadeGrid,num_alpha,num_beta,_s));
				}
			}
			fprintf (stdout," ---- Grid computed for AA: ",AAString[residue],"\n");
			fprintf (gridInfoFile,CLEAR_FILE, fadeGrid, "\n", gridInfo);
		}

	}	
}

// Phase 3
fprintf	(stdout,      "[PHASE 3] Computing Dirichlet weights using iterative algorithm.\n"); 
ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{PATH_TO_CURRENT_BF[0][Abs(PATH_TO_CURRENT_BF)-2],"FADE_PHASE_3_iterative.bf"}}), {"0": "" + basePath, "1": "" +  concentration});

// Phase 4
fprintf	(stdout,      "[PHASE 4] Computing posteriors and tabulating results.\n"); 
ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{PATH_TO_CURRENT_BF[0][Abs(PATH_TO_CURRENT_BF)-2],"FADE_PHASE_4_iterative.bf"}}), {"0": "" + basePath,"1": "" + basePath+".allresults.csv", "2": "" + basePath+".webdata.mat"});

<<<<<<< HEAD
// GENERATING WEB RESULTS PAGE
ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{PATH_TO_CURRENT_BF[0][Abs(PATH_TO_CURRENT_BF)-2],"FADE_Processor.ibf"}}), {"0": "" + basePath, "1": "" + basePath+".webdata.mat"});
=======
	// Phase 4 iterative
	fprintf	(stdout,      "[PHASE 4] Computing posteriors and tabulating results.\n"); 
	ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{PATH_TO_CURRENT_BF[0][Abs(PATH_TO_CURRENT_BF)-2],"FADE_PHASE_4_iterative.bf"}}), {"0": "" + LAST_FILE_PATH, "1": "" + LAST_FILE_PATH+"_allresults.csv", "2": LAST_FILE_PATH+"_summary.txt", "3": "" + LAST_FILE_PATH+"_webdata.mat"});
}
>>>>>>> 950db49cfbfd36939ed4131851b7d3662277c733

//THESE FUNCTIONS MOSTLY HELP WRITE OUT GRID FILES
function vectorToMatrix(columnVector, rows)
{
	npoints = Rows(columnVector);
	if(npoints % rows != 0)
	{
		fprintf(stdout, npoints, " is not divisible by ", rows,".\n");
		return 0;
	}
	
	columns = npoints / rows;
	matrixret = {rows, columns};
	index = 0;
	for(_x = 0 ; _x < rows ; _x += 1)
	{
		for(_y = 0 ; _y < columns ; _y += 1)
		{
			matrixret[_x][_y] = columnVector[index];
			index += 1;
		}
	}

	return matrixret;
}

function vectorToMatrixCSVstring(columnVector, rows)
{
	npoints = Rows(columnVector);
	if(npoints % rows != 0)
	{
		fprintf(stdout, "vectorToMatrix:", npoints, " is not divisible by ", rows,".\n");
		return 0;
	}
	
	columns = npoints / rows;
	matrixret = "";
	index = 0;
	for(_x = 0 ; _x < rows ; _x += 1)
	{
		matrixret = matrixret+(1+_x);
		for(_y = 0 ; _y < columns ; _y += 1)
		{
			matrixret = matrixret + ","+columnVector[index];
			index += 1;
		}
		matrixret = matrixret+"\n";
	}

	return matrixret;
}

function saveGridForSite(gridFileName, conditionalsGrid, num_alpha, num_beta, col)
{	
	fprintf(gridFileName,CLEAR_FILE);			
	fprintf(gridFileName,AAString[residue],",",num_alpha,",",num_beta,"\n");
	index = 0;
	for(_x = 0 ; _x < num_alpha ; _x += 1)
	{
		for(_y = 0 ; _y < num_beta ; _y += 1)
		{
			fprintf(gridFileName,conditionalsGrid[sites*index+col]);
			if(_y < num_beta - 1)
			{
				fprintf(gridFileName,",");
			}
			index += 1;
		}
		fprintf(gridFileName,"\n");
	}
	return 0;
}

function getGridStringForSite(conditionalsGrid, grid, num_alpha, num_beta, col)
{	
	ret_string = "";
	
	ret_string += ",";
	index = 0;
	for(_x = 0 ; _x < num_alpha ; _x += 1)
	{
		for(_y = 0 ; _y < num_beta ; _y += 1)
		{
			if(_x == 0)
			{
				ret_string += grid[index][1];
				if(_y < num_beta - 1)
				{
					ret_string += ",";
				}
			}
			index += 1;
		}
	}

	ret_string += "\n";
	index = 0;
	for(_x = 0 ; _x < num_alpha ; _x += 1)
	{
		ret_string += grid[index][0];
		ret_string += ",";
		for(_y = 0 ; _y < num_beta ; _y += 1)
		{
			ret_string += conditionalsGrid[sites*index+col];
			if(_y < num_beta - 1)
			{
				ret_string += ",";
			}
			index += 1;
		}
		ret_string += "\n";
	}
	return ret_string;
}
