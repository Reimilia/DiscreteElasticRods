#include "CollisionDetection.h"
#include "RigidBodyInstance.h"
#include "RigidBodyTemplate.h"
#include <Eigen/Core>
#include "VectorMath.h"
#include <Eigen/Geometry>
#include <iostream>

using namespace std;

struct BBox
{
    double mins[3];
    double maxs[3];
};

bool intersects(const BBox &b1, const BBox &b2)
{
    for(int i=0; i<3; i++)
    {
        if( (b1.maxs[i] < b2.mins[i]) || (b2.maxs[i] < b1.mins[i]))
            return false;
    }
    return true;
}

struct AABBNode
{
    AABBNode() : left(NULL), right(NULL), childtet(-1) {}
    ~AABBNode() {delete left; delete right;}

    AABBNode *left;
    AABBNode *right;
    BBox box;
    int childtet;
};

class NodeComparator
{
public:
    NodeComparator(int axis) : axis(axis) {}

    int axis;

    bool operator()(const AABBNode *left, const AABBNode *right) const
    {
        return left->box.mins[axis] < right->box.mins[axis];
    }
};

AABBNode *buildAABB(vector<AABBNode *> nodes)
{
    if(nodes.size() == 0)
        return NULL;
    else if(nodes.size() == 1)
        return nodes[0];

    double axismins[3];
    double axismaxs[3];
    for(int i=0; i<3; i++)
    {
        axismins[i] = numeric_limits<double>::infinity();
        axismaxs[i] = -numeric_limits<double>::infinity();
    }
    int nnodes = (int)nodes.size();
    for(int i=0; i<nnodes; i++)
    {
        for(int j=0; j<3; j++)
        {
            axismins[j] = min(axismins[j], nodes[i]->box.mins[j]);
            axismaxs[j] = max(axismaxs[j], nodes[i]->box.maxs[j]);
        }
    }
    double widths[3];
    for(int i=0; i<3; i++)
        widths[i] = axismaxs[i] - axismins[i];
    int splitaxis = -1;
    if(widths[0] >= widths[1] && widths[0] >= widths[2])
        splitaxis = 0;
    else if(widths[1] >= widths[0] && widths[1] >= widths[2])
        splitaxis = 1;
    else
        splitaxis = 2;
    std::sort(nodes.begin(), nodes.end(), NodeComparator(splitaxis));
    vector<AABBNode *> left(nnodes/2);
    vector<AABBNode *> right(nnodes - nnodes/2);
    for(int i=0; i<nnodes/2; i++)
    {
        left[i] = nodes[i];
    }
    for(int i=nnodes/2; i<nnodes; i++)
    {
        right[i-nnodes/2] = nodes[i];
    }
    AABBNode *node = new AABBNode;
    node->left = buildAABB(left);
    node->right = buildAABB(right);
    for(int i=0; i<3; i++)
    {
        node->box.mins[i] = min(node->left->box.mins[i], node->right->box.mins[i]);
        node->box.maxs[i] = max(node->left->box.maxs[i], node->right->box.maxs[i]);
    }
    return node;
}

void refitAABB(const RigidBodyInstance * instance, AABBNode *node)
{
    if(node->childtet != -1)
    {
        Eigen::Vector4i tet = instance->getTemplate().getTets().row(node->childtet);
        for(int k=0; k<3; k++)
        {
            node->box.mins[k] = numeric_limits<double>::infinity();
            node->box.maxs[k] = -numeric_limits<double>::infinity();
        }
        for(int k=0; k<4; k++)
        {
            Eigen::Vector3d point = instance->c + VectorMath::rotationMatrix(instance->theta)*instance->getTemplate().getVerts().row(tet[k]).transpose();
            for(int l=0; l<3; l++)
            {
                node->box.mins[l] = min(node->box.mins[l], point[l]);
                node->box.maxs[l] = max(node->box.maxs[l], point[l]);
            }
        }
    }
    else if(node->left && node->right)
    {
        refitAABB(instance, node->left);
        refitAABB(instance, node->right);
        for(int i=0; i<3; i++)
        {
            node->box.mins[i] = min(node->left->box.mins[i], node->right->box.mins[i]);
            node->box.maxs[i] = max(node->left->box.maxs[i], node->right->box.maxs[i]);
        }
    }
}

AABBNode *buildAABB(const RigidBodyInstance * instance)
{    
    int ntets = (int)instance->getTemplate().getTets().rows();
    vector<AABBNode *> leaves(ntets);
    for(int j=0; j<ntets; j++)
    {
        AABBNode *leaf = new AABBNode;
        leaf->childtet = j;
        Eigen::Vector4i tet = instance->getTemplate().getTets().row(j);
        BBox box;
        for(int k=0; k<3; k++)
        {
            box.mins[k] = numeric_limits<double>::infinity();
            box.maxs[k] = -numeric_limits<double>::infinity();
        }
        for(int k=0; k<4; k++)
        {
            Eigen::Vector3d point = instance->c + VectorMath::rotationMatrix(instance->theta)*instance->getTemplate().getVerts().row(tet[k]).transpose();
            for(int l=0; l<3; l++)
            {
                box.mins[l] = min(box.mins[l], point[l]);
                box.maxs[l] = max(box.maxs[l], point[l]);
            }
        }
        leaf->box = box;
        leaves[j] = leaf;
    }
    return buildAABB(leaves);
}

bool vertInTet(const Eigen::Vector3d &p, const Eigen::Vector3d &q1, const Eigen::Vector3d &q2, const Eigen::Vector3d &q3, const Eigen::Vector3d &q4)
{
    if( (q2-p).cross(q3-p).dot(q4-p) < 0)
        return false;
    if( (p-q1).cross(q3-q1).dot(q4-q1) < 0)
        return false;
    if( (q2-q1).cross(p-q1).dot(q4-q1) < 0)
        return false;
    if( (q2-q1).cross(q3-q1).dot(p-q1) < 0)
        return false;
    return true;
}

void tetTetIntersect(const AABBNode *node1, const AABBNode *node2, int body1, int body2, const std::vector<RigidBodyInstance *> instances, std::set<Collision> &collisions)
{
    if(body1==body2)
        return;

    Eigen::Vector4i tet1 = instances[body1]->getTemplate().getTets().row(node1->childtet);
    Eigen::Vector4i tet2 = instances[body2]->getTemplate().getTets().row(node2->childtet);
    Eigen::Vector3d verts1[4];
    Eigen::Vector3d verts2[4];
    for(int i=0; i<4; i++)
    {
        verts1[i] = instances[body1]->c + VectorMath::rotationMatrix(instances[body1]->theta)*instances[body1]->getTemplate().getVerts().row(tet1[i]).transpose();
        verts2[i] = instances[body2]->c + VectorMath::rotationMatrix(instances[body2]->theta)*instances[body2]->getTemplate().getVerts().row(tet2[i]).transpose();
    }
    for(int i=0; i<4; i++)
    {
        if(vertInTet(verts1[i], verts2[0], verts2[1], verts2[2], verts2[3]))
        {
            Collision c;
            c.body1 = body1;
            c.body2 = body2;
            c.collidingVertex = tet1[i];
            c.collidingTet = node2->childtet;
            collisions.insert(c);
        }
        if(vertInTet(verts2[i], verts1[0], verts1[1], verts1[2], verts1[3]))
        {
            Collision c;
            c.body1 = body2;
            c.body2 = body1;
            c.collidingVertex = tet2[i];
            c.collidingTet = node1->childtet;
            collisions.insert(c);
        }
    }
}

void intersect(const AABBNode *node1, const AABBNode *node2, int body1, int body2, const std::vector<RigidBodyInstance *> instances, std::set<Collision> &collisions)
{
    if(!node1 || !node2)
        return;

    if(!intersects(node1->box, node2->box))
        return;

    if(node1->childtet != -1)
    {
        if(node2->childtet != -1)
        {
            tetTetIntersect(node1, node2, body1, body2, instances, collisions);
        }
        else
        {
            intersect(node1, node2->left, body1, body2, instances, collisions);
            intersect(node1, node2->right, body1, body2, instances, collisions);
        }
    }
    else
    {
        if(node2->childtet != -1)
        {
            intersect(node1->left, node2, body1, body2, instances, collisions);
            intersect(node1->right, node2, body1, body2, instances, collisions);
        }
        else
        {
            intersect(node1->left, node2->left, body1, body2, instances, collisions);
            intersect(node1->left, node2->right, body1, body2, instances, collisions);
            intersect(node1->right, node2->left, body1, body2, instances, collisions);
            intersect(node1->right, node2->right, body1, body2, instances, collisions);
        }
    }
}

void collisionDetection(const std::vector<RigidBodyInstance *> instances, std::set<Collision> &collisions)
{
    collisions.clear();
    int nbodies = (int)instances.size();

    for(int i=0; i<nbodies; i++)
    {
        refitAABB(instances[i], instances[i]->AABB);
    }

    for(int i=0; i<nbodies; i++)
    {
        for(int j=i+1; j<nbodies; j++)
        {
            intersect(instances[i]->AABB, instances[j]->AABB, i, j, instances, collisions);
        }
    }

    // floor

    for (int i=0; i<nbodies; i++)
    {
        int nverts = instances[i]->getTemplate().getVerts().rows();
        for (int j = 0; j < nverts; j++)
        {
            Eigen::Vector3d point = VectorMath::rotationMatrix(instances[i]->theta) * instances[i]->getTemplate().getVerts().row(j).transpose() + instances[i]->c;
            if (point[1] <= -1.0)
            {
                Collision c;
                c.body1 = i;
                c.body2 = -1;
                c.collidingTet = -1;
                c.collidingVertex = j;
                collisions.insert(c);
            }
        }
    }
}
