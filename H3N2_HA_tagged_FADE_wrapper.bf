inputRedirect = {};
inputRedirect["01"]="New analysis"; // reload base file
inputRedirect["02"]="../datasets/H3N2_HA/H3N2_HA.fas"; // alignment file
inputRedirect["03"]="Flu H5N1"; // influenza amino-acid rate matrix
inputRedirect["04"]="/home/michael/Development/FADE/datasets/H3N2_HA/H3N2_HA_tagged.nwk"; // tagged tree
inputRedirect["05"]="Fix 1st sequence as root"; // use specified root
inputRedirect["06"]="2"; // number of MCMC chains
inputRedirect["07"]="500000"; // number of MCMC iterations
inputRedirect["08"]="250000"; // burn-in
inputRedirect["09"]="100"; //number of samples to draw
inputRedirect["10"]="0.5"; // Dirichlet concentration parameters
ExecuteAFile ("src/FADE.bf", inputRedirect);
