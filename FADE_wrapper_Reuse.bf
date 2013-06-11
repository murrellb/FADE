inputRedirect = {};
inputRedirect["01"]="Reload"; // reload base file
inputRedirect["02"]="../datasets/lysin/lysinAA.both.base"; // base file name
inputRedirect["03"]="Fix 1st sequence as root"; // use specified root
inputRedirect["04"]="0.5"; // Dirichlet concentration parameters
ExecuteAFile ("src/FADE.bf", inputRedirect);
