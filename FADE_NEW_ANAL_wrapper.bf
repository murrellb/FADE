inputRedirect = {};
inputRedirect["01"]="New analysis"; // reload base file
inputRedirect["02"]="../datasets/lysin/lysinAA.both"; // alignment file
inputRedirect["03"]="LG"; // amino-acid rate matrix
inputRedirect["04"]="Y"; // use that tree
inputRedirect["05"]="Fix 1st sequence as root"; // use specified root
inputRedirect["06"]="0.5"; // Dirichlet concentration parameters
ExecuteAFile ("src/FADE.bf", inputRedirect);
