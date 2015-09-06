#include <vector>
#include <cassert>
#include <initializer_list>
#include <iostream>

#include "kdTree.h"

template<size_t N>
class ndPoint
{
public:
    static const size_t dimension = N;

    ndPoint()
    :mPnt(dimension)
    {
    }

    ndPoint(std::initializer_list<double> l)
    :mPnt(l)
    {
        assert(l.size() == dimension);
    }

    void set(std::vector<double> pnt)
    {
        assert(pnt.size() == dimension);
        mPnt=pnt;
    }

    double getDimension(size_t dim) const
    {
        assert(dim < dimension);
        return mPnt[dim];
    }

    friend std::ostream & operator << ( std::ostream & output, const ndPoint & pnt)
    {
        for(size_t i=0;i<dimension-1;++i)
        {
            output<<pnt.mPnt[i]<<",";
        }
        output<<pnt.mPnt[dimension-1]<<std::endl;
        return output;
    }

private:
    std::vector<double> mPnt;
};

int main(void)
{
    typedef ndPoint<4> pointType;
    KDTree<pointType> tree;

    std::vector<pointType> pts{{0,0,0,0},
        {1,1,1,1},
        {-1, 3, 4,3},
        {5, 6, 7,5},
        {2, -6, 8,6},
        {4, 5, -4,-3},
        {2, 3, 4,-1},
        {2, 5, 6,-2}
    };

    tree.buildTree(pts);

    tree.dumpTreeInorder();

    std::cout << "searching near 0,0,0.1" << std::endl;
    auto closeNodes = tree.getPointsWithinCube({0, 0,0, 0.1}, 0.2);
    std::cout << "found" << std::endl;
    for(auto n : closeNodes){
      tree.dumpNode(n);
    }

    std::cout << "searching near 0.5,0.5,0.5" << std::endl;
    closeNodes = tree.getPointsWithinCube({0.5, 0.5, 0.5,0.5}, 0.7);
    std::cout << "found" << std::endl;
    for(auto n : closeNodes){
      tree.dumpNode(n);
    }

    std::cout << "min, x: " << std::endl;
    tree.dumpNode(tree.findMin(0));

    std::cout << "min, y: " << std::endl;
    tree.dumpNode(tree.findMin(1));
    std::cout << "min, z: " << std::endl;
    tree.dumpNode(tree.findMin(2));

    for(auto i = 0; i < pts.size() -1; ++i){
      std::cout << "deleting node " << i << std::endl;
      tree.deletePoint(i);
      tree.dumpTreeInorder();
    }

    std::cout << "inserting 2 points" << std::endl;
    tree.insertPoint({1, 4, 5,3});
    tree.dumpTreeInorder();
    tree.insertPoint({3, 8, 6,-4});
    tree.dumpTreeInorder();

    return 0;
}
