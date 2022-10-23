#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
using namespace std;

const double dx = 1.0;
const double dy = 1.0;
const double dxy = 0.0;
const double pi = 3.1415926535898;
const double a = 1;

double C(double x, double y)
{
	return sin(a*x) * sin(a*y);
}

double source(double x, double y)
{
	return -a*a * (2.*dxy * cos(a*x)*cos(a*y) - (dx+dy) * sin(a*x)*sin(a*y));
}

enum BoundCondType
{
	BC_DIR = 1,
	BC_NEUM = 2
};

// Class including everything needed
class Problem
{
private:
	/// Mesh
	Mesh &m;
	// =========== Tags =============
	/// Solution tag: 1 real value per node
	Tag tagConc;
	/// Diffusion tensor tag: 3 real values (Dx, Dy, Dxy) per cell
	Tag tagD;
	/// Boundary condition type tag: 1 integer value per node
	Tag tagBCtype;
	/// Boundary condition value tag: 1 real value per node, sparse on nodes
	Tag tagBCval;
	/// Right-hand side tag: 1 real value per node, sparse on nodes
	Tag tagSource;
	/// Analytical solution tag: 1 real value per node
	Tag tagConcAn;
	/// Global index tag: 1 integer value per node
	Tag tagGlobInd;

	// =========== Tag names ===========
	const string tagNameConc = "Concentration";
	const string tagNameD = "Diffusion_tensor";
	const string tagNameBCtype = "BC_type";
	const string tagNameBCval = "BC_value";
	const string tagNameSource = "Source";
	const string tagNameConcAn = "Concentration_analytical";
	const string tagNameGlobInd = "Global_Index";

	// =========== Markers
	/// Marker for Dirichlet nodes
	MarkerType mrkDirNode;
	/// Number of Dirichlet nodes
	unsigned numDirNodes;
public:
	Problem(Mesh &m_);
	~Problem();
	void initProblem();
	void assembleGlobalSystem(Sparse::Matrix &A, Sparse::Vector &rhs);
	void assembleLocalSystem(const Cell &c, rMatrix &A_loc, rMatrix &rhs_loc);
	void run();
};

Problem::Problem(Mesh &m_) : m(m_)
{
}

Problem::~Problem()
{

}

void Problem::initProblem()
{
	// Init tags
	tagConc = m.CreateTag(tagNameConc, DATA_REAL, NODE, NONE, 1);
	tagD = m.CreateTag(tagNameD, DATA_REAL, CELL, NONE, 3);
	tagBCtype = m.CreateTag(tagNameBCtype, DATA_INTEGER, NODE, NODE, 1);
	tagBCval = m.CreateTag(tagNameBCval, DATA_REAL, NODE, NODE, 1);
	tagSource =  m.CreateTag(tagNameSource, DATA_REAL, NODE, NONE, 1);
	tagConcAn =  m.CreateTag(tagNameConcAn, DATA_REAL, NODE, NONE, 1);
	tagGlobInd = m.CreateTag(tagNameGlobInd, DATA_INTEGER, NODE, NONE, 1);

	// Cell loop
	// 1. Check that cell is a triangle
	// 2. Set diffusion tensor values
	for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
		Cell c = icell->getAsCell();
		ElementArray<Node> nodes = c.getNodes();
		if (nodes.size() != 3){
			cout<<"Cell is not a triangle!!!";
			exit(1);
		}
		
		c.RealArray(tagD)[0] = dx;
		c.RealArray(tagD)[1] = dy;
		c.RealArray(tagD)[2] = dxy;
	}

	// Node loop
	// 1. Write analytical solution and source
	// 2. Check if node is boundary, set BC type and value, mark it
	// 3. Assign global indices
	mrkDirNode = m.CreateMarker();
	numDirNodes = 0;
	int glob_ind = 0;
	for(Mesh::iteratorNode inode = m.BeginNode(); inode != m.EndNode(); inode++){
		Node n = inode->getAsNode();
		double xn[3];
		n.Centroid(xn);
		n.Real(tagConcAn) = C(xn[0], xn[1]);
		n.Real(tagSource) = source(xn[0], xn[1]);

		if(n.Boundary()){
			n.SetMarker(mrkDirNode);
			n.Integer(tagBCtype) = BC_DIR;
			numDirNodes++;
			n.Real(tagBCval) = n.Real(tagConcAn);

		}
		else{
			n.Integer(tagGlobInd) = glob_ind;
			glob_ind++;
		}
	}
}

void Problem::assembleLocalSystem(const Cell &c, rMatrix &A_loc, rMatrix &rhs_loc)
{

}

void Problem::assembleGlobalSystem(Sparse::Matrix &A, Sparse::Vector &rhs)
{
	// Cell loop
	// For each cell assemble local system
	// and incorporate it into global
	for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
		Cell c = icell->getAsCell();
		rMatrix A_loc, rhs_loc;
		assembleLocalSystem(c, A_loc, rhs_loc);
		// Now A_loc is 3x3, rhs_loc is 3x1
		//  
		//
		//
		//

		ElementArray<Node> nodes = c.getNodes();
		unsigned glob_ind[3];
		for(unsigned loc_ind = 0; loc_ind < 3; loc_ind++){
			glob_ind[loc_ind] = static_cast<unsigned>(nodes[loc_ind].Integer(tagGlobInd));
		}
		
		for(unsigned loc_ind = 0; loc_ind < 3; loc_ind++){
			// Consider node with local index 'loc_ind'
			for(unsigned j = 0; j < 3; j++)
				A[glob_ind[loc_ind]][glob_ind[j]] += A_loc(loc_ind,j);
		}
	}
}

void Problem::run()
{
	// Matrix size
	unsigned N = static_cast<unsigned>(m.NumberOfNodes()) - numDirNodes;
	// Global matrix called 'stiffness matrix'
	Sparse::Matrix A;
	// Solution vector
	Sparse::Vector sol;
	// Right-hand side vector
	Sparse::Vector rhs;

	A.SetInterval(0, N);
	sol.SetInterval(0, N);
	rhs.SetInterval(0, N);

	assembleGlobalSystem(A, rhs);

	A.Save("A.mtx");
	rhs.Save("rhs.mtx");

	string solver_name = "inner_mptiluc";
	Solver S(solver_name);

	S.SetMatrix(A);
	bool solved = S.Solve(rhs, sol);
	if(!solved){
		printf("Linear solver failed: %s\n", S.GetReason().c_str());
		printf("Number of iterations: %d\n", S.Iterations());
		printf("Residual:             %e\n", S.Residual());
		exit(1);
	}

	for(Mesh::iteratorNode inode = m.BeginNode(); inode != m.EndNode(); inode++){
		Node n = inode->getAsNode();
		if(n.GetMarker(mrkDirNode))
			continue;
		unsigned ind = static_cast<unsigned>(n.Integer(tagGlobInd));
		n.Real(tagConc) = sol[ind];
	}

	m.Save("res.vtk");
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
	Problem P(m);
	P.initProblem();
	P.run();
	printf("Success\n");
	return 0;
}
