#ifndef _BJ_KD_TREE_H
#define _BJ_KD_TREE_H

#include <vector>
#include <memory>
#include <numeric>
#include <iostream>

template <typename PointType, int SplitDimension>
class KDNode {
public:
    typedef KDNode<PointType, (SplitDimension + 1)%PointType::dimension> ChildType;

    KDNode(size_t ind) : treeIndex(ind) {}
    size_t treeIndex; //particle's position in tree's list
    std::unique_ptr<ChildType> leftChild, rightChild;
};


template <typename PointType, typename PointArray=std::vector<PointType> >
class KDTree {
public:
    KDTree() {}
    KDTree(const PointArray& pointsIn){
        buildTree(pointsIn);
    }

    void buildTree(const PointArray& pointsIn);
    void dumpTreeInorder();

    template <typename F>
    void inorderTraversal(F func);

    std::vector<size_t> getPointsWithinCube(PointType testPoint, double radius);
    size_t findMin(int dimension);

    void dumpNode(size_t i){
        std::cout << mPoints[i] << std::endl;
    }

    void deletePoint(size_t nodeIndex);
    PointType getPoint(size_t nodeIndex){return mPoints[nodeIndex];}
    void insertPoint(const PointType& p);


    //END PUBLIC API
private:
    std::unique_ptr<KDNode<PointType, 0> > mRoot;
    PointArray mPoints;
    std::vector<size_t> mPointIndeces;

    template <int SplitDimension>
    std::unique_ptr< KDNode<PointType, SplitDimension> >
    buildSubtree(std::vector<size_t>::iterator begin, std::vector<size_t>::iterator end );

    template<int SplitDimension>
    void dumpSubtree(std::unique_ptr< KDNode<PointType, SplitDimension> >& node);

    template<int SplitDimension>
    void getPointsWithinCubeSubtree( PointType testPoint, double queryRange[2*PointType::dimension],
        std::unique_ptr< KDNode<PointType, SplitDimension> >& node, std::vector<size_t>& ret);

    static bool pointInRange(PointType testPoint, double queryRange[2*PointType::dimension]){
        for(int i = 0; i < PointType::dimension; ++i){
            if(testPoint.getDimension(i) < queryRange[2*i] ||
                testPoint.getDimension(i) > queryRange[2*i +1]){
                //std::cout << "point not in range: " << testPoint << " dim: " << i << std::endl;
                return false;
            }
        }
        return true;
    }

    template<int SplitDimension>
    size_t findMinSubtree(int dimension, std::unique_ptr<KDNode<PointType, SplitDimension> >& node);

    template<int SplitDimension>
    std::unique_ptr<KDNode<PointType, SplitDimension> >
    deleteFromSubtree(size_t nodeIndex,std::unique_ptr<KDNode<PointType, SplitDimension> >& node);

    template<int SplitDimension, typename F>
    void inorderTraversalSubtree(F func, std::unique_ptr<KDNode<PointType, SplitDimension> >& node);


    template<int SplitDimension>
    std::unique_ptr<KDNode<PointType, SplitDimension> >
    insertPointSubtree(std::unique_ptr<KDNode<PointType, SplitDimension> >& node,size_t pointIndex);

};


template<typename PointType, typename PointArray>
void KDTree<PointType, PointArray>::buildTree(const PointArray& pointsIn){
    mPoints = pointsIn;
    mPointIndeces.resize( mPoints.size() );
    std::iota(begin(mPointIndeces), end(mPointIndeces), 0);
    mRoot = buildSubtree<0>(begin(mPointIndeces), end(mPointIndeces));
}

template<typename PointType, typename PointArray>
template<int SplitDimension>
std::unique_ptr<KDNode<PointType, SplitDimension> >
KDTree<PointType, PointArray>::buildSubtree( std::vector<size_t>::iterator begin,
    std::vector<size_t>::iterator end){

    auto rangeSize = std::distance(begin, end);

    if(rangeSize == 0){
        return std::unique_ptr<KDNode<PointType, SplitDimension> >(nullptr);
    } else {
        std::sort(begin, end,
            [this]( size_t a, size_t b){
                return mPoints[a].getDimension(SplitDimension) < mPoints[b].getDimension(SplitDimension);
            });
        auto median = begin + rangeSize/2;
        while(median != begin  &&
            mPoints[*(median)].getDimension(SplitDimension) ==
                mPoints[*(median - 1)].getDimension(SplitDimension)){
            --median;
            //put all the nodes with equal coord value in the right subtree
        }
        auto ret = std::unique_ptr<KDNode<PointType, SplitDimension> >
            ( new KDNode<PointType, SplitDimension>(*median));
        ret->leftChild  = buildSubtree<(SplitDimension +1)%PointType::dimension>(begin, median);
        ret->rightChild = buildSubtree<(SplitDimension +1)%PointType::dimension>(median + 1, end);
        return ret;
    }
}


template<typename PointType, typename PointArray>
void KDTree<PointType, PointArray>::dumpTreeInorder(){
    dumpSubtree<0>(mRoot);
}

template<typename PointType, typename PointArray>
template<int SplitDimension>
void KDTree<PointType, PointArray>::dumpSubtree(std::unique_ptr<KDNode<PointType, SplitDimension> >& node){
    if(node->leftChild){
        std::cout << "dumping left: " << std::endl;
        dumpSubtree<(SplitDimension +1)%PointType::dimension %PointType::dimension>(node->leftChild);
    }
    std::cout << "dumping this: " << std::endl;
    std::cout << node->treeIndex << ": " << mPoints[node->treeIndex] << std::endl;
    if(node->rightChild){
        std::cout << "dumping right: " << std::endl;
        dumpSubtree<(SplitDimension +1)%PointType::dimension %PointType::dimension>(node->rightChild);
    }
}

template<typename PointType, typename PointArray>
std::vector<size_t> KDTree<PointType, PointArray>::getPointsWithinCube(PointType testPoint, double radius){

    double queryRange[2*PointType::dimension];
    for(auto i = 0; i < PointType::dimension; ++i){
        queryRange[2*i] = testPoint.getDimension(i) - radius;
        queryRange[2*i +1] = testPoint.getDimension(i) + radius;
    }

    std::vector<size_t> ret;
    getPointsWithinCubeSubtree<0>(testPoint, queryRange, mRoot, ret);

    return ret;
}

template<typename PointType, typename PointArray>
template<int SplitDimension>
void KDTree<PointType, PointArray>::getPointsWithinCubeSubtree(PointType testPoint,
    double queryRange[2*PointType::dimension],std::unique_ptr<KDNode<PointType,
    SplitDimension> >& node,std::vector<size_t>& ret){

    if(node == nullptr){
        return;
    }

    //std::cout << "query range: " << std::endl;
    //for(int i= 0; i < PointType::dimension; ++i){
    //  std::cout << queryRange[2*i] << ' ' << queryRange[2*i +1] << std::endl;
    //}

    auto nodePoint = mPoints[node->treeIndex];
    if(pointInRange(nodePoint, queryRange)){
        ret.push_back(node->treeIndex);
    }
    if(nodePoint.getDimension(SplitDimension) >= queryRange[2*SplitDimension]){
        //query range goes into the left subtree
        //std::cout << "recurse left" << std::endl;
        getPointsWithinCubeSubtree<(SplitDimension +1)%PointType::dimension>(testPoint,
            queryRange,node->leftChild,ret);
    }
    if(nodePoint.getDimension(SplitDimension) <= queryRange[2*SplitDimension + 1]){
        //query range goes into the right subtree
        //std::cout << "recurse right" << std::endl;
        getPointsWithinCubeSubtree<(SplitDimension +1)%PointType::dimension>(testPoint,
            queryRange,node->rightChild,ret);
    }
}


template<typename PointType, typename PointArray>
size_t KDTree<PointType, PointArray>::findMin(int dimension){
    return findMinSubtree<0>(dimension, mRoot);
}


template<typename PointType, typename PointArray>
template<int SplitDimension>
size_t KDTree<PointType, PointArray>::findMinSubtree(int dimension,
    std::unique_ptr<KDNode<PointType,SplitDimension> >& node){

    if(SplitDimension == dimension){
        if(node->leftChild == nullptr){
            return node->treeIndex;
        } else {
            return findMinSubtree<(SplitDimension+1)%PointType::dimension>(dimension,
                node->leftChild);
        }
    } else {
        size_t leftMin = 123456, rightMin= 123456;
        if(node->leftChild){
            leftMin = findMinSubtree<(SplitDimension+1)%PointType::dimension>(dimension,
                node->leftChild);
        }
        if(node->rightChild){
            rightMin = findMinSubtree<(SplitDimension+1)%PointType::dimension>(dimension,
                node->rightChild);
        }
        auto nodeValue = mPoints[node->treeIndex].getDimension(dimension);
        if(node->leftChild && mPoints[leftMin].getDimension(dimension) < nodeValue){
            if(node->rightChild){
                return ( mPoints[leftMin].getDimension(dimension) <
                    mPoints[rightMin].getDimension(dimension) ) ? leftMin : rightMin;
            } else {
                return leftMin;
            }
        } else if(node->rightChild && mPoints[rightMin].getDimension(dimension) < nodeValue){
            return rightMin;
        } else {
            return node->treeIndex;
        }
    }
}

template<typename PointType, typename PointArray>
  void KDTree<PointType, PointArray>::deletePoint(size_t nodeIndex){

  mRoot = deleteFromSubtree<0>(nodeIndex, mRoot);
}

template<typename PointType, typename PointArray>
  template<int SplitDimension>
  std::unique_ptr<KDNode<PointType, SplitDimension> >
  KDTree<PointType, PointArray>::deleteFromSubtree(size_t nodeIndex,
                           std::unique_ptr<KDNode<PointType, SplitDimension> >& node){

  constexpr size_t nextDimension = (SplitDimension +1)%PointType::dimension;

  if(node->treeIndex == nodeIndex){
    if(node->rightChild){
      auto rightMin = findMinSubtree<nextDimension>(SplitDimension, node->rightChild);
      node->treeIndex = rightMin;
      node->rightChild = deleteFromSubtree<nextDimension>(rightMin,
                       node->rightChild);
    } else if(node->leftChild){
      auto leftMin = findMinSubtree<nextDimension>(SplitDimension, node->leftChild);
      node->treeIndex = leftMin;
      node->rightChild = deleteFromSubtree<nextDimension>(leftMin,
                              node->leftChild);
      node->leftChild = nullptr;
    } else {
      return nullptr;
    }
  } else if(mPoints[nodeIndex].getDimension(SplitDimension) <
        mPoints[node->treeIndex].getDimension(SplitDimension)){

    node->leftChild = deleteFromSubtree<nextDimension>(nodeIndex,
                    node->leftChild);
  } else {
    node->rightChild = deleteFromSubtree<nextDimension>(nodeIndex,
                            node->rightChild);
  }
  return std::move(node);

}

template<typename PointType, typename PointArray>
template <typename F>
void KDTree<PointType, PointArray>::inorderTraversal(F func){

    inorderTraversalSubtree<0, F>(func, mRoot);
}

template<typename PointType, typename PointArray>
template<int SplitDimension, typename F>
void KDTree<PointType, PointArray>::inorderTraversalSubtree(F func,
    std::unique_ptr<KDNode<PointType,SplitDimension> >& node){

    auto constexpr nextDimension = (SplitDimension +1)%PointType::dimension;
    if(node->leftChild){
        inorderTraversalSubtree<nextDimension, F>(func, node->leftChild);
    }
    func(mPoints[node->treeIndex]);
    if(node->rightChild){
        inorderTraversalSubtree<nextDimension, F>(func, node->rightChild);
    }
}

template<typename PointType, typename PointArray>
void KDTree<PointType, PointArray>::insertPoint(const PointType& point){
    mPoints.push_back(point);
    mRoot = insertPointSubtree<0>(mRoot, mPoints.size() -1);
}

template<typename PointType, typename PointArray>
template<int SplitDimension>
std::unique_ptr<KDNode<PointType, SplitDimension> >
KDTree<PointType, PointArray>::insertPointSubtree(std::unique_ptr<KDNode<PointType,
    SplitDimension> >& node,size_t pointIndex){

    auto constexpr nextDimension = (SplitDimension +1)%PointType::dimension;

    if(node == nullptr){
        std::cout << "new node" << std::endl;
        return std::unique_ptr<KDNode<PointType, SplitDimension> > (new KDNode<PointType, SplitDimension>(pointIndex));
    } else if (mPoints[pointIndex].getDimension(SplitDimension) <
        mPoints[node->treeIndex].getDimension(SplitDimension)){

        std::cout << "adding left" << std::endl;
        node->leftChild = insertPointSubtree<nextDimension>(node->leftChild, pointIndex);
        std::cout << "added left " << std::endl;
    } else {
        std::cout << "adding right " << std::endl;
        node->rightChild = insertPointSubtree<nextDimension>(node->rightChild, pointIndex);
        std::cout << "added right" << std::endl;
    }
    return std::move(node);
}
#endif //_BJ_KD_TREE_H
