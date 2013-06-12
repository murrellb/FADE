/* define an associative array with key = amino acid and value = integer
	index for AA in alphabetical order */
AAString    = "ACDEFGHIKLMNPQRSTVWY";
             //01234567890123456789


//run_residue = 16; // run only this residue, -1 => run all residues
run_residue = -1;
num_alpha = 20; // number of site-to-site rate variation grid points
num_beta = 20; // number of bias grid points
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

// Phase 1: User input in this section, optimize call in FADE_PHASE_1_iterative.bf
if (reloadFlag == 0) // optimize baseline model
{
	/* new analysis, fit baseline model */
	
	/* this should be a protein alignment */	
	SetDialogPrompt ("\nSelect an alignment file:");
	DataSet			ds 			 = ReadDataFile (PROMPT_FOR_FILE);
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
	fprintf(stdout,"\n");
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

	/* export baseline model LF without optimization*/
	LIKELIHOOD_FUNCTION_OUTPUT = 7;
	outPath = basePath + ".base";
	fprintf (outPath, CLEAR_FILE, lf, CLOSE_FILE); // write out unoptimized base-line model

	fprintf	(stdout,      "[PHASE 1.1] Optimizing the tree branch lengths under the baseline model.\n");
	// optimizes the baseline model and overwrites existing baseline file
	ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{PATH_TO_CURRENT_BF[0][Abs(PATH_TO_CURRENT_BF)-2],"FADE_PHASE_1_iterative.bf"}}), {"0":outPath});	
}
else
{
	fprintf (stdout, "Locate an existing fit (.base):");
	fscanf(stdin, "String", outPath);
	fprintf							(stdout,      "[PHASE 1.1] Loading the baseline model.\n");
}

// Load the optimized baseline
modelNameString = "_customAAModelMatrix";

ExecuteAFile (outPath);
GetString (lfInfo,lf,-1);
if ((lfInfo["Models"])[0] == "mtREVModel")
{
	modelNameString = "mtREVMatrix";	
}
bpSplit						 = splitFilePath (outPath);
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


if(reloadFlag != 0)
{
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
}
LFCompute (lf,LF_START_COMPUTE);
LFCompute (lf,baselineLogL);
LFCompute (lf,LF_DONE_COMPUTE);

root_left  = "biasedTree." + (treeAVL[(rootNode["Children"])[0]])["Name"] + ".t";
root_right = "biasedTree." + (treeAVL[(rootNode["Children"])[1]])["Name"] + ".t";

treeContainsTags = (treeString||"{FG}")[0][0] != -1;

/* vector of branch lengths for entire tree */
baselineBL						= BranchLength (givenTree,-1);
referenceL						= (baselineBL * (Transpose(baselineBL)["1"]))[0];

fprintf	(stdout,      "[PHASE 1.2] Standard model fit. Log-L = ",baselineLogL,". Tree length = ",referenceL, " subs/site \n"); 

ExecuteAFile 					(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "GrabBag.bf");

fixGlobalParameters ("lf");
/*------------------------------------------------------------------------------*/
//THIS IS ACTUALLY A CONCENTRATION MULTIPLIER NOW, SINCE THE PRIOR IS A VECTOR
fadeConcentrationParameter = prompt_for_a_value ("The concentration parameter of the Dirichlet prior",0.5,0.001,1,0);    
fprintf (stdout, "[DIAGNOSTIC] FADE will use the Dirichlet prior concentration parameter of ", fadeConcentrationParameter, "\n"); 


// Phase 2
fprintf	(stdout,      "[PHASE 2] Computing conditional likelihoods.\n"); 
ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{PATH_TO_CURRENT_BF[0][Abs(PATH_TO_CURRENT_BF)-2],"FADE_PHASE_2_iterative.bf"}}), {"0":outPath,"1":""+num_alpha,"2":""+num_beta,"3":""+treeContainsTags});

// Phase 3
fprintf	(stdout,      "[PHASE 3] Computing Dirichlet weights using iterative algorithm.\n"); 
ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{PATH_TO_CURRENT_BF[0][Abs(PATH_TO_CURRENT_BF)-2],"FADE_PHASE_3_iterative.bf"}}), {"0": "" + basePath, "1": "" +  concentration});

// Phase 4
fprintf	(stdout,      "[PHASE 4] Computing posteriors and tabulating results.\n"); 
ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{PATH_TO_CURRENT_BF[0][Abs(PATH_TO_CURRENT_BF)-2],"FADE_PHASE_4_iterative.bf"}}), {"0": "" + basePath,"1": "" + basePath+".allresults.csv", "2": "" + basePath+".webdata.mat", "3": "" + basePath+".summary.txt"});

// GENERATING WEB RESULTS PAGE
ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{PATH_TO_CURRENT_BF[0][Abs(PATH_TO_CURRENT_BF)-2],"FADE_Processor.ibf"}}), {"0": "" + basePath, "1": "" + basePath+".webdata.mat"});

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
