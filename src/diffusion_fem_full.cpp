#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
using namespace std;

//const double dx = 1.0;
//const double dy = 1.0;
//const double dxy = 0.0;
// Corresponds to tensor
// [ 1  0 ]
// [ 0 10 ]
// rotated by M_PI/6
const double dx = 10;//3.25;
const double dy = 6;//-0.433013;
const double dxy = 2*sqrt(3);//0.25;
const double pi = 3.1415926535898;
const double a = 1;

double C(double x, double y)
{
        return 0;//sin(a*x) * sin(a*y);
}

double source(double x, double y)
{
        double r = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));
	if(r < 0.1)
	        return -10;
	return 0;//-a*a * (2.*dxy * cos(a*x)*cos(a*y) - (dx+dy) * sin(a*x)*sin(a*y));
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
	tagBCval = m.CreateTag(tagNameBCval, DATA_REAL, NODE, NONE, 1);
	tagSource = m.CreateTag(tagNameSource, DATA_REAL, NODE, NONE, 1);
	tagConcAn = m.CreateTag(tagNameConcAn, DATA_REAL, NODE, NONE, 1);
	tagGlobInd = m.CreateTag(tagNameGlobInd, DATA_INTEGER, NODE, NONE, 1);

	// Cell loop
	// 1. Check that cell is a triangle
	// 2. Set diffusion tensor values
	for(Mesh::iteratorCell icell = m.BeginCell(); icell != m.EndCell(); icell++){
	        Cell c = icell->getAsCell();
		ElementArray<Node> nodes = c.getNodes();
		if(nodes.size() != 3){
		        printf("Cell %d is not a triangle, has %llu nodes!",
			           c.LocalID(), nodes.size());
		}

		c.RealArray(tagD)[0] = dx; // Dx
		c.RealArray(tagD)[1] = dy; // Dy
		c.RealArray(tagD)[2] = dxy; // Dxy
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
			printf("BC val = %lf\n", n.Real(tagBCval));
			n.Real(tagConc) = n.Real(tagConcAn);
			n.Integer(tagGlobInd) = -1;
		}
		else{
		        n.Integer(tagGlobInd) = glob_ind;
			glob_ind++;
		}
	}
	printf("Total nodes: %d, Dir: %d, gi = %d\n", m.NumberOfNodes(), numDirNodes, glob_ind);
}

void Problem::assembleLocalSystem(const Cell &c, rMatrix &A_loc, rMatrix &rhs_loc)
{;
        ElementArray<Node> nodes = c.getNodes();

	double x0[2], x1[2], x2[2];
	nodes[0].Barycenter(x0);
	nodes[1].Barycenter(x1);
	nodes[2].Barycenter(x2);

	rMatrix D(2,2); // Diffusion tensor
	D(0,0) = c.RealArray(tagD)[0];
	D(1,1) = c.RealArray(tagD)[1];
	D(1,0) = c.RealArray(tagD)[2];
	D(0,1) = c.RealArray(tagD)[2];

	rMatrix Bk(2,2);
	Bk(0,0) = x1[0] - x0[0]; //x2 - x1;
	Bk(0,1) = x2[0] - x0[0]; //x3 - x1;
	Bk(1,0) = x1[1] - x0[1]; //y2 - y1;
	Bk(1,1) = x2[1] - x0[1]; //y3 - y1;

	rMatrix Ck = Bk.Invert() * D * Bk.Invert().Transpose();
	//Ck = Dk * Ck;

	double detBk = Bk(0,0)*Bk(1,1) - Bk(0,1)*Bk(1,0);

	rMatrix Kee(3,3), Knn(3,3), Ken(3,3);
	Kee.Zero();
	Knn.Zero();
	Ken.Zero();

	Kee(0,0) = Kee(1,1) =  1.;
	Kee(0,1) = Kee(1,0) = -1.;
	Knn(0,0) = Knn(2,2) =  1.;
	Knn(0,2) = Knn(2,0) = -1.;
	Ken(0,0) = Ken(1,2) =  1.;
	Ken(1,0) = Ken(0,2) = -1.;

	Kee *= 0.5;
	Knn *= 0.5;
	Ken *= 0.5;

	A_loc = rMatrix(3,3);
	A_loc = Ck(0,0)*Kee + Ck(1,1)*Knn + Ck(0,1)*(Ken + Ken.Transpose());
	A_loc *= fabs(detBk);


	// RHS

	rMatrix res(3,1);

	res.Zero();
	res(0,0) += nodes[0].Real(tagSource) + nodes[1].Real(tagSource) + nodes[2].Real(tagSource);
	res(1,0) = res(0,0);
	res(2,0) = res(0,0);

	rhs_loc = res * fabs(detBk) / 18.;
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

		ElementArray<Node> nodes = c.getNodes();
		// This loop goes over 3 nodes
		// using 'i' as local index in A_loc
		// In A, indexing is done through LocalID()
		// In LocalID(), 'Local' means local in the
		// sense of MPI-parallel computations

		unsigned ind[3];
		for(unsigned i = 0; i < 3; i++)
		        ind[i] = static_cast<unsigned>(nodes[i].Integer(tagGlobInd));

		for(unsigned i = 0; i < 3; i++){
		        if(nodes[i].GetMarker(mrkDirNode)){
			        // There's no row corresponding to nodes[0]
			        double bcVal = nodes[i].Real(tagBCval);
				for(unsigned j = 0; j < 3; j++){
				        if(j == i)
					        continue;
					if(!nodes[j].GetMarker(mrkDirNode))
					        rhs[ind[j]] -= bcVal * A_loc(i,j);
				}
			}
			else{
			        A[ind[i]][ind[0]] += A_loc(i,0);
				A[ind[i]][ind[1]] += A_loc(i,1);
				A[ind[i]][ind[2]] += A_loc(i,2);

				rhs[ind[i]] += rhs_loc(i,0);
			}
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
	Mesh::GeomParam table;
	table[BARYCENTER] = CELL;
	table[CENTROID] = CELL | NODE;
	table[MEASURE] = CELL;
	m.PrepareGeometricData(table);
	Problem P(m);
	P.initProblem();
	P.run();
	printf("Success\n");
	return 0;
}
