fscanf              (stdin, "String", nuc_fit_file);
fscanf              (stdin, "String", grid_file);
fscanf              (stdin, "String", sample_base_file);
fscanf              (stdin, "Number", _chainCount);
fscanf              (stdin, "String", results_file);

_fubar_do_simulations = 0;

ExecuteAFile        (PATH_TO_CURRENT_BF + "FUBAR_tools.ibf");
LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("WriteDelimitedFiles");

/* This causes grid info files to be overwritten??
ExecuteAFile        (nuc_fit_file);
sequences = nucData_1.species;
GetInformation (treeCount,"^nuc_tree_[0-9]+$");
fileCount       = Columns (treeCount);
treeLengths     = {fileCount,1};

//getting branch lengths
for (fileID = 1; fileID <= fileCount; fileID += 1)
{
	treeLengths [fileID-1] = + Eval("BranchLength(nuc_tree_"+fileID+",-1)");
}
*/
fscanf (grid_file, REWIND, "NMatrix,Raw", grid, site_probs);
site_probs = Eval (site_probs);
sites   = Columns (site_probs["conditionals"]);


readMCMCSamples (sample_base_file,_chainCount);


pointsWithNoBias = {points,1} ["grid[_MATRIX_ELEMENT_ROW_][1]==1"];
//pointsWithNoBias = {points,1} ["grid[_MATRIX_ELEMENT_ROW_][1]>1"];
pointsWithNoBiasCount     = +pointsWithNoBias;

priorMean            = {1, points};
sampleFromThisDistro = {pointsWithNoBiasCount,2};

tabulateGridResults (points, sites, samples, _chainCount);

from = 0;
for (_point = 0; _point < points; _point += 1) {
    priorMean [_point] = (+jointSamples[-1][_point])/samples;
    if (pointsWithNoBias [_point]) {
        sampleFromThisDistro [from][0] = _point;
        sampleFromThisDistro [from][1] = priorMean [_point];
        from += 1;
    }
}

priorNN = +(sampleFromThisDistro [-1][1]);
fubar_results = reportSiteResults   (sites, 0, priorNN, _fubar_do_simulations);
fubarRowCount     = Rows (fubar_results);


site_counter = {};
for (currentFubarIndex = 0; currentFubarIndex < fubarRowCount; currentFubarIndex += 1) {
    site_counter + (currentFubarIndex+1);
}


WriteSeparatedTable (results_file, {{"Codon","alpha","bias","Prob[bias>1]", "Prob[bias=1]", "BayesFactor","PSRF", "Neff"}}, fubar_results, site_counter, ",");

