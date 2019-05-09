#ifndef RIGIDBODYTEMPLATE_H
#define RIGIDBODYTEMPLATE_H

#include <string>
#include <Eigen/Core>
#include <set>
#include <vector>

class SignedDistanceField;

class RigidBodyTemplate
{
public:
    RigidBodyTemplate(const std::string &meshFilename, double scale);
    RigidBodyTemplate(const Eigen::MatrixX3d &verts, const Eigen::MatrixX4i &tets);
    ~RigidBodyTemplate();

    double getVolume() const {return volume_;}
    const Eigen::Matrix3d getInertiaTensor() const {return inertiaTensor_;}    
    
    const Eigen::MatrixX3d &getVerts() const {return V;}
    const Eigen::MatrixX3i &getFaces() const {return F;}      
    const Eigen::MatrixX4i &getTets() const { return T; }

    double distance(Eigen::Vector3d p, int tet) const;
	double signedDistance(Eigen::Vector3d p, int tet) const;
    Eigen::Vector3d Ddistance(int tet) const;
	Eigen::Vector3d signedDdistance(int tet) const;

private:
    RigidBodyTemplate(const RigidBodyTemplate &other) = delete;
    RigidBodyTemplate &operator=(const RigidBodyTemplate &other) = delete;

    void initialize();
    
    void computeFaces();
    void computeVolume();
    Eigen::Vector3d computeCenterOfMass();
    void computeInertiaTensor();
    void computeDistances();

    
    Eigen::MatrixX3d V;
    Eigen::MatrixX3i F;
    Eigen::MatrixX4i T;

    std::vector<double> distances;
    
    double volume_;
    Eigen::Matrix3d inertiaTensor_;    
};

#endif // RIGIDBODYTEMPLATE_H
