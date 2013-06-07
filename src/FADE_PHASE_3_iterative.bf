filepath =  LAST_FILE_PATH;
bf_path = PATH_TO_CURRENT_BF;
AAString    = "ACDEFGHIKLMNPQRSTVWY";

fscanf              (stdin, "Number", _concentration);
aminoAcids = 20;

debug = 0;

if (!debug) {
    funcText = "";

    funcsToExport = {"0": "runIterativeDeterministic"};
    for (k = 0; k < Abs (funcsToExport); k+=1) {
	funcText += exportFunctionDefinition (funcsToExport[k]);
    }


    funcText += "\nfor (_chainIndex=start; _chainIndex<end; _chainIndex+=1) {runIterativeDeterministic(filepath, bf_path,_chainIndex,_concentration);} return 0;";

    variablesToExport =  "filepath = \"" + filepath + "\";\n" + 
			 "bf_path = \"" + bf_path + "\";\n" + 
			 "_sampleFile = \"" + _sampleFile + "\";\n" + 
			 "_gridInfo = \"" + _gridInfo + "\";\n" +   
			 "_concentration = " + _concentration + ";\n";

     if (MPI_NODE_COUNT > 1) {
	    per_node    = Max(1,aminoAcids $ MPI_NODE_COUNT);
	    _startPoint = aminoAcids-per_node;
	    leftover    = aminoAcids-per_node*MPI_NODE_COUNT;
	    
	    from          = 0;
	    to            = per_node + (leftover>0);
	    node_ranges   = {MPI_NODE_COUNT,2};
	    
	    for (node_id = 1; node_id < Min(aminoAcids,MPI_NODE_COUNT); node_id += 1) {
	                                
		
	        MPISend				(node_id, variablesToExport + ";start = " +from + ";end=" + to+";" + funcText); 
	        
	        
	        node_ranges [node_id][0]         = from;
	        node_ranges [node_id][1]         = to;
	        
	        from                             = to;
	        to                              += per_node+(node_id<=leftover);  
	    } 
	} else {
	_startPoint = 0;    
    }
	
    for (_r = _startPoint; _r < aminoAcids; _r += 1){
	runIterativeDeterministic(filepath,bf_path,_r,_concentration);
    }
    if (MPI_NODE_COUNT > 1 && points > MPI_NODE_COUNT) {
	for (node_id = 1; node_id < Min(aminoAcids,MPI_NODE_COUNT); node_id += 1) {
	    MPIReceive (-1,fromNode,res);
	}
    }
}


function runIterativeDeterministic(file_path, bf_path, residue, _concentration)
{
	AAString    = "ACDEFGHIKLMNPQRSTVWY";
	
	_gridInfo = file_path +"."+AAString[residue]+".grid_info";
	_sampleFile = file_path+"."+AAString[residue]+".theta";
	
	ExecuteAFile        ( bf_path + "FUBAR_tools_iterative.ibf");

	fscanf (_gridInfo, REWIND, "NMatrix,Raw", grid, gridInfo);
	gridInfo = Eval(gridInfo);

	points             = Rows(grid);
	sites              = Columns(gridInfo["conditionals"]);
	normalize_by_site  = ({1,points}["1"])*(gridInfo["conditionals"]);
	normalized_weights = (gridInfo["conditionals"])*({sites,sites}["1/normalize_by_site[_MATRIX_ELEMENT_ROW_]*(_MATRIX_ELEMENT_ROW_==_MATRIX_ELEMENT_COLUMN_)"]);
	sum_by_site        = normalized_weights * ({sites,1}["1"]);


	countBetaOnes = +({points,1} ["grid[_MATRIX_ELEMENT_ROW_][1]==1"]);
	priorvec = {points,1} ["_concentration*(1+(countBetaOnes-2)*(grid[_MATRIX_ELEMENT_ROW_][1]==1))"]; // points[i] = 1, except where beta = 1, points[i] = 19
	//fprintf(stdout,"priorvec",priorvec,"\n");

	weights = {1,points}["1"];
	weights = weights * (1/(+weights));
	oldweights = Transpose(weights);

	diffSum = 1;
	iters=1;
	t0 = Time (1);
	while (diffSum > 0.0000000001) {
	 phiUN = {points,sites}["normalized_weights[_MATRIX_ELEMENT_ROW_][_MATRIX_ELEMENT_COLUMN_]*weights[_MATRIX_ELEMENT_ROW_]"];
	 phiNormalizers  = ({1,points}["1"])*phiUN;
	 phi = phiUN*({sites,sites}["1/phiNormalizers[_MATRIX_ELEMENT_ROW_]*(_MATRIX_ELEMENT_ROW_==_MATRIX_ELEMENT_COLUMN_)"]);
	weights = (phi * ({sites,1}["1"]))+priorvec;
	// weights = (phi * ({sites,1}["1"]))+_concentration;
	 weights = weights * (1/(+weights));
	 diffVec = weights - oldweights;
	 diffSum = +({points,1}["diffVec[_MATRIX_ELEMENT_ROW_]*diffVec[_MATRIX_ELEMENT_ROW_]"]);
	 SetParameter (STATUS_BAR_STATUS_STRING, "Iteration: "+ iters + " ------ delta:" + diffSum + " ----- Time:" + _formatTimeString(Time(1)-t0),0);
	 oldweights = weights;
	 iters = iters+1;
	 }


	fprintf (stdout, "\nTime taken:" + _formatTimeString(Time(1)-t0), "\n");
	convergefile = _sampleFile;
	fprintf (convergefile,CLEAR_FILE, weights);
	

	return 0;
}
