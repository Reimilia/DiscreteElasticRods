#ifndef RIGIDBODYINSTANCE_H
#define RIGIDBODYINSTANCE_H

#include <Eigen/Core>
#include <list>
#include <vector>

class RigidBodyTemplate;
struct AABBNode;

class RigidBodyInstance
{
public:
    RigidBodyInstance(const RigidBodyTemplate &rbtemplate, const Eigen::Vector3d &c, const Eigen::Vector3d &theta, const Eigen::Vector3d &cvel, const Eigen::Vector3d &w, double density);
    ~RigidBodyInstance();

    Eigen::Vector3d c;
    Eigen::Vector3d theta;

    Eigen::Vector3d cvel;
    Eigen::Vector3d w;

    double density;

    AABBNode *AABB;
    
    const RigidBodyTemplate &getTemplate() const {return rbtemplate_;}
    
private:
    const RigidBodyTemplate &rbtemplate_;
};

#endif // RIGIDBODYINSTANCE_H
