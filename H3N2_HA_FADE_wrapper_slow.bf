inputRedirect = {};
inputRedirect["01"]="Reload"; // reload base file
inputRedirect["02"]="../datasets/H3N2_HA/H3N2_HA.fas.base"; // base file name
inputRedirect["03"]="Fix 1st sequence as root"; // use specified root
inputRedirect["04"]="3"; // number of MCMC chains
inputRedirect["05"]="5000000"; // number of MCMC iterations
inputRedirect["06"]="2000000"; // burn-in
inputRedirect["07"]="500"; //number of samples to draw
inputRedirect["08"]="1"; // Dirichlet concentration parameters
ExecuteAFile ("src/FADE.c", inputRedirect);
