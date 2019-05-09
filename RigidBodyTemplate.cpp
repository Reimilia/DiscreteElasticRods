#include "RigidBodyTemplate.h"
#include <iostream>
#include <igl/readOBJ.h>
#include <Eigen/Dense>
#include <fstream>
#include <map>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <Eigen/Sparse>
#include "Distance.h"

using namespace std;
using namespace Eigen;

RigidBodyTemplate::RigidBodyTemplate(const std::string &meshFilename, double scale) : volume_(0)
{
    inertiaTensor_.setZero();
    Eigen::MatrixXd mV;
    Eigen::MatrixXi mF;
    igl::readOBJ(meshFilename, mV, mF);

    mV *= scale;

    igl::copyleft::tetgen::tetrahedralize(mV, mF, "pq1.414a0.01", V, T, F);
    computeFaces();
    initialize();
}

RigidBodyTemplate::RigidBodyTemplate(const Eigen::MatrixX3d &verts, const Eigen::MatrixX4i &tets) : volume_(0)
{
    V = verts;
    T = tets;
    computeFaces();
    initialize();
}

RigidBodyTemplate::~RigidBodyTemplate()
{    
}

void RigidBodyTemplate::initialize()
{
    computeVolume();
    Vector3d cm = computeCenterOfMass();
    for(int i=0; i<V.rows(); i++)
        V.row(i) -= cm;

    computeInertiaTensor();    
    computeDistances();
}

void RigidBodyTemplate::computeFaces()
{

    struct triple
    {
        triple(int aa, int bb, int cc) : a(aa), b(bb), c(cc)
        {
            if(a < b)
                std::swap(a,b);
            if(a < c)
                std::swap(a,c);
            if(b < c)
                std::swap(b,c);
        }

        int a, b, c;
        bool operator<(const triple &other) const
        {
            if(a < other.a)
                return true;
            else if(a > other.a)
                return false;
            if(b < other.b)
                return true;
            else if(b > other.b)
                return false;
            return c < other.c;
        }
    };

    int ntets = (int)T.rows();
    MatrixX3i allfaces(4*ntets, 3);
    Matrix<int, 4, 3> faceidx;
    faceidx << 0, 1, 3,
        3, 1, 2,
        3, 2, 0,
        0, 2, 1;

    for(int i=0; i<ntets; i++)
    {
        Vector4i tet = T.row(i);
        for(int face=0; face<4; face++)
        {
            for(int k=0; k<3; k++)
                allfaces(4*i+face, k) = tet[faceidx(face,k)];
        }
    }

    map<triple, vector<int> > faces;
    for(int i=0; i<4*ntets; i++)
    {
        triple t(allfaces(i,0), allfaces(i,1), allfaces(i,2));
        faces[t].push_back(i/4);
    }

    int nfaces=0;
    for(map<triple, vector<int> >::iterator it = faces.begin(); it != faces.end(); ++it)
        if(it->second.size() == 1)
            nfaces++;

    F.resize(nfaces,3);
    int idx=0;

    for(int i=0; i<4*ntets; i++)
    {
        triple t(allfaces(i,0), allfaces(i,1), allfaces(i,2));
        if(faces[t].size() == 1)
        {
            F.row(idx) = allfaces.row(i);
            idx++;
        }        
    }    
}

void RigidBodyTemplate::computeVolume()
{
    volume_ = 0;
    for (int i = 0; i < F.rows(); i++)
    {
        Vector3d pts[3];
        Vector3d centroid(0, 0, 0);
        for (int j = 0; j < 3; j++)
        {
            pts[j] = V.row(F(i, j));
            centroid += pts[j];
        }
        Vector3d normal = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
        double area = 0.5 * normal.norm();
        normal /= normal.norm();

        centroid /= 3.0;
        volume_ += centroid.dot(normal) * area / 3.0;
    } 
}

Vector3d RigidBodyTemplate::computeCenterOfMass()
{
    Vector3d cm(0, 0, 0);
    for (int i = 0; i < F.rows(); i++)
    {
        Vector3d pts[3];
        for (int j = 0; j < 3; j++)
        {
            pts[j] = V.row(F(i, j));
        }
        Vector3d normal = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
        double area = 0.5 * normal.norm();
        normal /= normal.norm();

        Vector3d term(0, 0, 0);
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                for (int l = k; l < 3; l++)
                    term[j] += pts[k][j] * pts[l][j];
            }
            term[j] *= area*normal[j] / 12.0;
        }

        cm += term;
    }

    return cm / volume_;
}

void RigidBodyTemplate::computeInertiaTensor()
{
    Vector3d quads(0, 0, 0);
    Vector3d mixed(0, 0, 0);
    for (int i = 0; i < F.rows(); i++)
    {
        Vector3d pts[3];
        for (int j = 0; j < 3; j++)
        {
            pts[j] = V.row(F(i, j));
        }
        Vector3d normal = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
        double area = 0.5 * normal.norm();
        normal /= normal.norm();


        Vector3d term(0, 0, 0);
        Vector3d mixterm(0, 0, 0);
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                for (int l = k; l < 3; l++)
                {
                    for (int m = l; m < 3; m++)
                        term[j] += pts[k][j] * pts[l][j] * pts[m][j];
                }
            }
            term[j] *= area*normal[j] / 30.0;
        }
        double mix = 0;
        for (int j = 0; j < 3; j++)
        {
            mix += 6.0*pts[j][0] * pts[j][1] * pts[j][2];
            for (int k = 0; k < 3; k++)
            {
                mix += 2.0*pts[j][k] * pts[j][(k + 1) % 3] * pts[(j + 1) % 3][(k + 2) % 3];
                mix += 2.0*pts[j][k] * pts[j][(k + 1) % 3] * pts[(j + 2) % 3][(k + 2) % 3];
            }
            mix += pts[j][0] * pts[(j + 1) % 3][1] * pts[(j + 2) % 3][2];
            mix += pts[j][2] * pts[(j + 1) % 3][1] * pts[(j + 2) % 3][0];
        }
        for (int j = 0; j < 3; j++)
            mixterm[j] = mix*area*normal[j] / 60.0;

        quads += term;
        mixed += mixterm;
    }

    inertiaTensor_ << quads[1] + quads[2], -mixed[2], -mixed[1],
        -mixed[2], quads[0] + quads[2], -mixed[0],
        -mixed[1], -mixed[0], quads[0] + quads[1];
}

void RigidBodyTemplate::computeDistances()
{
    int nverts = (int)V.rows();
    int nfaces = (int)F.rows();
    distances.resize(nverts);
    for (int i = 0; i < nverts; i++)
    {
        double dist = numeric_limits<double>::infinity();
        for (int j = 0; j < nfaces; j++)
        {
            double dummy;
            if (Distance::vertexPlaneDistanceLessThan(V.row(i), V.row(F(j, 0)), V.row(F(j, 1)), V.row(F(j, 2)), dist))
            {
                Vector3d distvec = Distance::vertexFaceDistance(V.row(i), V.row(F(j, 0)), V.row(F(j, 1)), V.row(F(j, 2)), dummy, dummy, dummy);
                dist = min(dist, distvec.norm());
            }
        }
        distances[i] = dist;        
    }
}

double RigidBodyTemplate::distance(Vector3d p, int tet) const
{
	// TODO compute distance from point to object boundary
	Matrix3d A;
	A.row(0) = V.row(T(tet, 0)) - V.row(T(tet, 3));
	A.row(1) = V.row(T(tet, 1)) - V.row(T(tet, 3));
	A.row(2) = V.row(T(tet, 2)) - V.row(T(tet, 3));
	A.transposeInPlace();
  
	Vector3d coef = A.colPivHouseholderQr().solve(p - V.row(T(tet, 3)).transpose());
  
	//std::cout << coef.transpose() << std::endl;
  
	/*if (coef[0] < 0 || coef[1] < 0 || coef[2] < 0 || coef.sum()>1)
	{
		std::cout << "Warning: the point does not lie in the tetrahedral!\n";
	}*/
	
  
	return coef[0] * distances[T(tet, 0)] + coef[1] * distances[T(tet, 1)] +
		coef[2] * distances[T(tet, 2)] + (1-coef.sum()) * distances[T(tet, 3)];
}
  
Vector3d RigidBodyTemplate::Ddistance(int tet) const
{
	//TODO: compute derivative of distance from point to boundary
	Vector3d result(0, 0, 0);
  
	Matrix3d A;
	A.row(0) = V.row(T(tet, 0)) - V.row(T(tet, 3));
	A.row(1) = V.row(T(tet, 1)) - V.row(T(tet, 3));
	A.row(2) = V.row(T(tet, 2)) - V.row(T(tet, 3));
	//A.transposeInPlace();
  
	Vector3d d(
		distances[T(tet, 0)]- distances[T(tet, 3)],
		distances[T(tet, 1)] - distances[T(tet, 3)],
		distances[T(tet, 2)] - distances[T(tet, 3)]
		);
  
	result = A.inverse() * d;  // £¨d.T * A.(-T)).T = A.(-1) * d
  
	return result;
}

double RigidBodyTemplate::signedDistance(Eigen::Vector3d p, int tet) const
{
	// Resolve ambiguity
	return -distance(p, tet);
}

Eigen::Vector3d RigidBodyTemplate::signedDdistance(int tet) const
{
	// Resolve ambiguity
	return -Ddistance(tet);
}
