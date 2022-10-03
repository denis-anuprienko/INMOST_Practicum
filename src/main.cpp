#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
using namespace std;

// Create three tags with different data types
// and defined on different mesh elements
void main_create_tag(Mesh &m)
{
	Tag tagNodeCellNum;  // for each node stores number of surrounding cells
	Tag tagFaceArea;     // face areas
	Tag tagCellCentroid; // cell centroids

	// 1. Name: "NumberOfCells"
	// 2. Data type: integers
	// 3. Element type: mesh nodes
	// 4. Sparsity: not needed yet, set to NONE
	// 5. Tag size: 1
	tagNodeCellNum = m.CreateTag("NumberOfCells", DATA_INTEGER, NODE, NONE, 1);


	// 1. Name: "FaceArea"
	// 2. Data type: real numbers
	// 3. Element type: mesh faces
	// 4. Sparsity: not needed yet, set to NONE
	// 5. Tag size: 1
	tagFaceArea = m.CreateTag("FaceArea", DATA_REAL, FACE, NONE, 1);

	// 1. Name: "CellCentroid"
	// 2. Data type: real numbers
	// 3. Element type: mesh cells
	// 4. Sparsity: not needed yet, set to NONE
	// 5. Tag size: 3
	tagCellCentroid = m.CreateTag("CellCentroid", DATA_REAL, CELL, NONE, 3);

	// Loops to fill tags

	// Node loop
	for(Mesh::iteratorNode inode = m.BeginNode(); inode != m.EndNode(); inode++){
		Node n = inode->getAsNode();
		int ncells = static_cast<int>(n.getCells().size());
		n.Integer(tagNodeCellNum) = ncells;
	}

	// Face loop
	for(Mesh::iteratorFace iface = m.BeginFace(); iface != m.EndFace(); iface++){
		Face f = iface->getAsFace();
		f.Real(tagFaceArea) = f.Area();
	}

	// Cell loop
	for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
		Cell c = icell->getAsCell();
		double x[3];
		c.Centroid(x);
		c.RealArray(tagCellCentroid)[0] = x[0];
		c.RealArray(tagCellCentroid)[1] = x[1];
		c.RealArray(tagCellCentroid)[2] = x[2];
	}
}

// Check if mesh has a tag
void main_get_tag(Mesh &m)
{
	string tname = "Solution";
	Tag tag;
	if(m.HaveTag(tname)){
		tag = m.GetTag(tname);
	}
	else{
		cout << "Mesh has no '" << tname << "' tag!" << endl;
		exit(1);
	}
}

// Compute distance between two nodes
double nodeDist(const Node &n1, const Node &n2)
{
	double x1[3], x2[3];
	n1.Centroid(x1);
	n2.Centroid(x2);
	return sqrt(
				(  x1[0] - x2[0])*(x1[0] - x2[0])
				+ (x1[1] - x2[1])*(x1[1] - x2[1])
				+ (x1[2] - x2[2])*(x1[2] - x2[2])
				);
}

// Compute cell diameter
double cellDiam(const Cell &c)
{
	ElementArray<Node> nodes = c.getNodes();
	unsigned m = static_cast<unsigned>(nodes.size());
	double diam = 0.;
	for(unsigned i = 0; i < m; i++){
		for(unsigned j = 0; j < m; j++){
			diam = max(diam, nodeDist(nodes[i], nodes[j]));
		}
	}
	return diam;
}

// Compute and print mesh diameter
void main_mesh_diam(Mesh &m)
{
	double diam = 0.0;
	for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
		diam = max(diam, cellDiam(icell->self()));
	}
	cout << "Mesh diameter is " << diam << endl;
}

// find difference between two functions
// which are cellwise-constant
void main_diff_cell_funcs(Mesh &m)
{
	Tag tagU, tagUapprox;
	string tname = "U";
	if(m.HaveTag(tname))
		tagU = m.GetTag(tname);
	else{
		cout << "Mesh has no '" << tname << "' tag!" << endl;
		exit(1);
	}
	tname = "U_approx";
	if(m.HaveTag(tname))
		tagUapprox = m.GetTag(tname);
	else{
		cout << "Mesh has no '" << tname << "' tag!" << endl;
		exit(1);
	}

	double normC  = 0.0;
	double normL2 = 0.0;
	for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
		Cell c = icell->getAsCell();
		double diff = fabs(c.Real(tagU) - c.Real(tagUapprox));
		normC = max(normC, diff);
		normL2 += diff * c.Volume(); // functions are constant on cell
	}
	normL2 = sqrt(normL2);
	cout << "|u - u_approx|_C  = " << normC  << endl;
	cout << "|u - u_approx|_L2 = " << normL2 << endl;
}

// Create linear system for nodal values
// Get solver and solve the system
// Write values to nodal tag
void main_linear_solver(Mesh &m)
{
	// Get number of nodes
	unsigned N = static_cast<unsigned>(m.NumberOfNodes());

	// Create sparse matrix, RHS vector and solution vector
	Sparse::Matrix A;
	Sparse::Vector b;
	Sparse::Vector sol;
	// Set their size - number of nodes
	A.SetInterval(0, N);
	b.SetInterval(0, N);
	sol.SetInterval(0, N);
	// Make A identity matrix
	// Make b: b_i = i
	for(Mesh::iteratorNode inode = m.BeginNode(); inode != m.EndNode(); inode++){
		Node n = inode->getAsNode();
		// Use 'local id' as index for matrix and vector
		unsigned i = static_cast<unsigned>(n.LocalID());
		A[i][i] = 1.0;
		b[i] = i;
	}

	// Get solver
	// All inner INMOST solvers are BiCGStab
	// with different preconditioners, let's use ILU2
	Solver S("inner_ilu2");
	S.SetParameter("absolute_tolerance", "1e-10");
	S.SetParameter("relative_tolerance", "1e-6");

	// Set matrix in the solver;
	// this also computes preconditioner
	S.SetMatrix(A);
	S.SetParameter("absolute_tolerance", "1e-10");
	S.SetParameter("relative_tolerance", "1e-");

	// Solve
	bool solved = S.Solve(b, sol);
	cout << "lin.it.: " << S.Iterations() << endl;
	if(!solved){
		cout << "Linear solver failure!" << endl;
		cout << "Reason: " << S.ReturnReason() << endl;
	}
	// Now we have solution in 'sol'
	// Create tag and save it here
	Tag tagSol = m.CreateTag("Sol", DATA_REAL, NODE, NONE, 1);
	for(Mesh::iteratorNode inode = m.BeginNode(); inode != m.EndNode(); inode++){
		Node n = inode->getAsNode();
		n.Real(tagSol) = sol[static_cast<unsigned>(n.LocalID())];
	}
}

int main(int argc, char ** argv)
{
	if( argc < 2 )
	{
		printf("Usage: %s mesh_file\n",argv[0]);
		return -1;
	}

	Mesh m;
	m.Load(argv[1]);
	main_create_tag(m);
	//main_get_tag(m);
	//main_mesh_diam(m);
	//main_diff_cell_funcs(m);
	//main_linear_solver(m);
	m.Save("res.vtk");
	return 0;
}
