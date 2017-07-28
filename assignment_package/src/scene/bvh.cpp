#include "bvh.h"
#include <QTime>




// Memory Allocation Functions
//void *AllocAligned(size_t size) {

//#if defined (PBRT_IS_WINDOWS)
//    return _aligned_malloc(size, PBRT_L1_CACHE_LINE_SIZE);
//#elif defined(PBRT_IS_OPENBSD) || defined(PBRT_IS_OSX) || defined(PBRT_IS_FREEBSD)
//    void *ptr;
//    if (posix_memalign(&ptr, PBRT_L1_CACHE_LINE_SIZE, size) != 0) ptr = nullptr;
//    return ptr;
////#else
////    return memalign(PBRT_L1_CACHE_LINE_SIZE, size);
//#endif
//}

//void FreeAligned(void *ptr) {
//    if (!ptr) return;
//#if defined(PBRT_IS_WINDOWS)
//    _aligned_free(ptr);
//#else
//    free(ptr);
//#endif
//}




// Feel free to ignore these structs entirely!
// They are here if you want to implement any of PBRT's
// methods for BVH construction.

struct BVHPrimitiveInfo {
    BVHPrimitiveInfo() {}
    BVHPrimitiveInfo(size_t primitiveNumber, const Bounds3f &bounds)
        : primitiveNumber(primitiveNumber),
          bounds(bounds),
          centroid(.5f * bounds.min + .5f * bounds.max) {}
    int primitiveNumber;
    Bounds3f bounds;
    Point3f centroid;
};

struct BVHBuildNode {
    // BVHBuildNode Public Methods
    void InitLeaf(int first, int n, const Bounds3f &b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
        children[0] = children[1] = nullptr;
    }
    void InitInterior(int axis, BVHBuildNode *c0, BVHBuildNode *c1) {
        children[0] = c0;
        children[1] = c1;
        bounds = Union(c0->bounds, c1->bounds);
        splitAxis = axis;
        nPrimitives = 0;
    }
    Bounds3f bounds;
    BVHBuildNode *children[2];
    int splitAxis, firstPrimOffset, nPrimitives;
};

struct MortonPrimitive {
    int primitiveIndex;
    unsigned int mortonCode;
};

struct LBVHTreelet {
    int startIndex, nPrimitives;
    BVHBuildNode *buildNodes;
};

struct LinearBVHNode {
    Bounds3f bounds;
    union {
        int primitivesOffset;   // leaf
        int secondChildOffset;  // interior
    };
    unsigned short nPrimitives;  // 0 -> interior node, 16 bytes
    unsigned char axis;          // interior node: xyz, 8 bytes
    unsigned char pad[1];        // ensure 32 byte total size
};


struct BucketInfo{
    int count = 0;
    Bounds3f bounds;
};


BVHAccel::~BVHAccel()
{
    delete [] nodes;
}

// Constructs an array of BVHPrimitiveInfos, recursively builds a node-based BVH
// from the information, then optimizes the memory of the BVH
BVHAccel::BVHAccel(const std::vector<std::shared_ptr<Primitive> > &p, int maxPrimsInNode)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)), primitives(p)
{

    //TODO

    // Test whether primitives is empty or not
    if(primitives.size() == 0){
        return;
    }


    QTime myTimer;
    std::cout << "BVH construction begins!" << std::endl;
    myTimer.start();


    // Build BVH from primitives

    // 1. Initialize permitiveInfo array for primitives
    std::vector<BVHPrimitiveInfo> primitiveInfo(primitives.size());
    for(size_t i = 0; i < primitives.size(); i++){
        primitiveInfo[i] = BVHPrimitiveInfo(i, primitives[i]->WorldBound());

    }



    // 2. Build BVH tree for primitives using primitiveInfo
    // MemoryArena arena(1024 * 1024);
    int totalNodes = 0;
    std::vector<std::shared_ptr<Primitive>> orderedPrims;
    orderedPrims.reserve(primitives.size());
    BVHBuildNode *root;

    root = recursiveBuild(primitiveInfo, 0, primitives.size(), &totalNodes, orderedPrims);

    primitives.swap(orderedPrims);


    // 3. Compute representation of depth-first traversal of BVH tree
    nodes = new LinearBVHNode[totalNodes];
    int offset = 0;
    flattenBVHTree(root, &offset);



    // 4. delete BVHBuildNode root
    deleteBuildNode(root);


    std::cout << "BVH construction ends, and total time used: "
              << myTimer.elapsed() << " milliseconds." << std::endl;

}

bool BVHAccel::Intersect(const Ray &ray, Intersection *isect) const
{
    //TODO

    if(!nodes) return false;

    bool hit = false;

//    bool SameLevelBoundingBoxIntersect = false;

    Vector3f invDir(1.0f / ray.direction.x, 1.0f / ray.direction.y, 1.0f / ray.direction.z);
    int dirIsNeg[3] = {invDir.x < 0.f, invDir.y < 0.f, invDir.z < 0.f};


    // Follow ray through BVH nodes to find primitive intersections
    int toVisitOffset = 0, currentNodeIndex = 0;
    int nodesToVisit[64];

    while(true){

        const LinearBVHNode *node = &nodes[currentNodeIndex];

        // Check ray against BVH node
        float temp_t = 0.f;
            // For the very first bounding box, we don't know whether the ray
            // interects (root node) or not, so wo do initial test here
            // but for the following bounding box added, we
            // make sure they are SURE to intersect!
        bool node_isect = (currentNodeIndex == 0) ? node->bounds.Intersect(ray, &temp_t) : true;

//      if(node_isect){
//      if(node->bounds.Intersect(ray, &temp_t)){
      if(node->bounds.IntersectP(ray, invDir, dirIsNeg)){

            // If this is a leaf node
            if(node->nPrimitives > 0){

                // Intersect ray with EVERY primitive in leaf BVH node
                for (int i = 0; i < node->nPrimitives; i++){
                    Intersection temp_isect;
                    if(primitives[node->primitivesOffset + i]->Intersect(ray, &temp_isect)){

                        hit = true;

                        // if iscet is still null,
                        // we need initialize it
                        if(isect->t == -1.0f && isect->objectHit == nullptr){
                            (*isect) = temp_isect;
                        }
                        else{
                            if(temp_isect.t < isect->t){
                                (*isect) = temp_isect;
                            }
                        }
                    }
                }

                // If it hits one primitive, loop break;
                // since bounding box is put in ToVisit list
                // according to t value(from small to large)

//                if(toVisitOffset == 0 || hit) break;

                // no.. actually, it's incorrect!
                // can't stop loop, there may be intersection with smaller t value!
                if(toVisitOffset == 0) break;


//                if(SameLevelBoundingBoxIntersect){
//                    if(toVisitOffset == 0) break;
//                }
//                else{
//                    if(toVisitOffset == 0 || hit) break;
//                }

                currentNodeIndex = nodesToVisit[--toVisitOffset];

            }

            // If this is a interior node
            else{

//                // Since this is a binary tree
//                // we do bounding box intersection test
//                // to its two children
//                float firstChild_t = 0.f;
//                bool firstChild_isect = nodes[currentNodeIndex + 1].bounds.Intersect(ray, &firstChild_t);
//                float secondChild_t = 0.f;
//                bool secondChild_isect = nodes[node->secondChildOffset].bounds.Intersect(ray, &secondChild_t);


//                // If it hits both children, we need sort them,
//                // put the one with smaller t value(closer) first!
//                // then the other with larger t value
//                if(firstChild_isect && secondChild_isect){

//                    if(!SameLevelBoundingBoxIntersect){
//                        SameLevelBoundingBoxIntersect = nodes[currentNodeIndex + 1].bounds.IntersectBoundingBox(
//                                    nodes[node->secondChildOffset].bounds.min,
//                                    nodes[node->secondChildOffset].bounds.max);
//                    }


//                    if(firstChild_t < secondChild_t){
//                        nodesToVisit[toVisitOffset++] = node->secondChildOffset;
//                        currentNodeIndex = currentNodeIndex + 1;
//                    }
//                    else{
//                        nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
//                        currentNodeIndex = node->secondChildOffset;
//                    }

//                }
//                // If we only intersects with one child,
//                // this child node will be next one to process
//                else if(firstChild_isect){
//                    currentNodeIndex = currentNodeIndex + 1;
//                }
//                else if(secondChild_isect){
//                    currentNodeIndex = node->secondChildOffset;
//                }
//                else{
//                    if(toVisitOffset == 0) break;
//                    currentNodeIndex = nodesToVisit[--toVisitOffset];
//                }


                if(dirIsNeg[node->axis]){
                    nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
                    currentNodeIndex = node->secondChildOffset;
                }
                else{
                    nodesToVisit[toVisitOffset++] = node->secondChildOffset;
                    currentNodeIndex = currentNodeIndex + 1;
                }



            }

        }

        // If the root node hits nothing
        else{
            if (toVisitOffset == 0) break;
            currentNodeIndex = nodesToVisit[--toVisitOffset];
//            break;
        }

    }

    return hit;
}



int BVHAccel::flattenBVHTree(BVHBuildNode *node, int *offset){

    LinearBVHNode *linearNode = &nodes[*offset];
    linearNode->bounds = node->bounds;
    int myOffset = (*offset)++;
    if(node->nPrimitives > 0){
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    }
    else{
        // Create interior flattened BVH node
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
        flattenBVHTree(node->children[0], offset);
        linearNode->secondChildOffset = flattenBVHTree(node->children[1], offset);
    }

    return myOffset;
}

BVHBuildNode *BVHAccel::recursiveBuild(std::vector<BVHPrimitiveInfo> &primitiveInfo,
    int start, int end, int *totalNodes,
    std::vector<std::shared_ptr<Primitive>> &orderedPrims){

    BVHBuildNode *node = new BVHBuildNode();
    //BVHBuildNode *node = arena.Alloc<BVHBuildNode>();

    (*totalNodes)++;

    // Compute bounds of all primitives in BVH node
    Bounds3f bounds = primitiveInfo[start].bounds;

    for(int i = start; i < end; i++){
        bounds = Union(bounds, primitiveInfo[i].bounds);
    }

    int nPrimitives = end - start;

    if(nPrimitives == 1){
        // Create leaf
        int firstPrimOffset = orderedPrims.size();
        for(int i = start; i < end; i++){
            int primNum = primitiveInfo[i].primitiveNumber;
            orderedPrims.push_back(primitives[primNum]);
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
        return node;
    }
    else{
        // Compute bound of primitive centroids, choose split dimension dim
        Bounds3f centroidBounds = Bounds3f(primitiveInfo[start].centroid);
        for(int i = start; i< end; i++){
            centroidBounds = Union(centroidBounds, primitiveInfo[i].centroid);
        }
        int dim = centroidBounds.MaximumExtent();


        // Partition primitives into two set and build children
        float mid = ((float)start + (float)end) / 2.0f;
        if(centroidBounds.max[dim] == centroidBounds.min[dim]){
            // Create leaf BVHBuildNode
            int firstPrimOffset = orderedPrims.size();
            for(int i = start; i < end; i++){
                int primNum = primitiveInfo[i].primitiveNumber;
                orderedPrims.push_back(primitives[primNum]);
            }

            node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
            return node;
        }
        else{
            // ------------------------------------------
            // Partition primitives using approximate SAH
            // if nPrimitives is smaller than 4,
            // we partition them into equally sized subsets
            // ------------------------------------------

            if(nPrimitives <= 2){
                //Partition primitives into equally sized subsets
                //mid = (start + end) / 2;
                mid = ((float)start + (float)end) / 2.0f;
                std::nth_element(&primitiveInfo[start], &primitiveInfo[(int)mid],
                                 &primitiveInfo[end - 1]+1,
                        [dim](const BVHPrimitiveInfo &a, const BVHPrimitiveInfo &b){
                                return a.centroid[dim] < b.centroid[dim];
                                });
            }

            else{
                // Allocate BucketInfo for SAH partition buckets
                constexpr int nBuckets = 12;


                BucketInfo buckets[nBuckets];

                // Initialize BucketInfo for SAH partition buckets
                for(int i = start; i < end; i++){
                    int b = nBuckets * centroidBounds.Offset(primitiveInfo[i].centroid)[dim];
                    if(b == nBuckets) b = nBuckets - 1;
                    buckets[b].count++;
                    buckets[b].bounds = Union(buckets[b].bounds, primitiveInfo[i].bounds);
                }

                // Compute costs for splitting after each bucket
                float cost[nBuckets - 1];
                for (int i = 0; i < nBuckets - 1; i++){
                    Bounds3f b0;
                    Bounds3f b1;
//                    Bounds3f b0, b1;
                    int count0 = 0;
                    int count1 = 0;
                    for(int j = 0; j <= i; j++){
                        b0 = Union(b0, buckets[j].bounds);
                        count0 += buckets[j].count;
                    }
                    for(int j = i + 1; j < nBuckets; j++){
                        b1 = Union(b1, buckets[j].bounds);
                        count1 += buckets[j].count;
                    }
                    cost[i] = 1.f + (count0 * b0.SurfaceArea() + count1 * b1.SurfaceArea()) / bounds.SurfaceArea();
                }

                // Find bucket to split at that minimizes SAH metric
                float minCost = cost[0];
                int minCostSplitBucket = 0;
                for (int i = 1; i < nBuckets - 1; i++){
                    if(cost[i] < minCost){
                        minCost = cost[i];
                        minCostSplitBucket = i;
                    }
                }

                // Either create leaf or split primitives at selected SAH bucket
                float leafCost = nPrimitives;
                if(nPrimitives > maxPrimsInNode || minCost < leafCost){
                    BVHPrimitiveInfo *pmid = std::partition(&primitiveInfo[start],
                                                            &primitiveInfo[end - 1] + 1,
                                                            [=](const BVHPrimitiveInfo &pi){
                                                                    int b = nBuckets * centroidBounds.Offset(pi.centroid)[dim];
                                                                    if(b == nBuckets) b = nBuckets - 1;
                                                                    return b <= minCostSplitBucket;
                                                            });
                    mid = pmid - &primitiveInfo[0];
                }
                else{
                    // Create leaf BVHBuildNode
                    int firstPrimOffset = orderedPrims.size();
                    for(int i = start; i < end; i++){
                        int primNum = primitiveInfo[i].primitiveNumber;
                        orderedPrims.push_back(primitives[primNum]);
                    }

                    node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
                    return node;
                }
            }
            // ------------------------------------------
            // --------- partition part end -------------
            // ------------------------------------------


            node->InitInterior(dim,
                               recursiveBuild(primitiveInfo, start, (int)mid, totalNodes, orderedPrims),
                               recursiveBuild(primitiveInfo, (int)mid, end, totalNodes, orderedPrims));


        }

    }
    return node;
}


void BVHAccel::deleteBuildNode(BVHBuildNode *root){

    // If this is a leaf node
    if(root->children[0] == nullptr &&
       root->children[1] == nullptr){
        delete root;
        return;
    }

    deleteBuildNode(root->children[0]);

    deleteBuildNode(root->children[1]);

    delete root;

    return;
}
