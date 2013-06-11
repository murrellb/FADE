inputRedirect = {};
inputRedirect["01"]="Reload"; // reload base file
inputRedirect["02"]="../datasets/H3N2_HA/H3N2_HA.fas.base"; // base file name
inputRedirect["03"]="Fix 1st sequence as root"; // use specified root
inputRedirect["04"]="2"; // number of MCMC chains
inputRedirect["05"]="500000"; // number of MCMC iterations
inputRedirect["06"]="250000"; // burn-in
inputRedirect["07"]="100"; //number of samples to draw
inputRedirect["08"]="0.5"; // Dirichlet concentration parameters
ExecuteAFile ("src/FADE.bf", inputRedirect);
