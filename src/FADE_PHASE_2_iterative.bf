fscanf		    (stdin, "String", baseFile);
fscanf              (stdin, "Number", num_alpha);
fscanf              (stdin, "Number", num_beta);
fscanf              (stdin, "Number", treeContainsTags);

ExecuteAFile (baseFile); // so backgroundMatrix gets loaded into memory

fadeGrid = defineFadeGrid (num_alpha, num_beta);
for (residue = 0; residue < 20; residue = residue + 1) 
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
