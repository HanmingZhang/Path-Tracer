#pragma once
#include "geometry/primitive.h"


// Forward declarations of structs used by our BVH tree
// They are defined in the .cpp file
struct BVHBuildNode;
struct BVHPrimitiveInfo;
struct MortonPrimitive;
struct LinearBVHNode;

class MemoryArena;

class BVHAccel : public Primitive
{
    //Functions
public:

    BVHAccel(const std::vector<std::shared_ptr<Primitive>> &p,
             int maxPrimsInNode = 1);
    Bounds3f WorldBound() const;
    ~BVHAccel();
    bool Intersect(const Ray &ray, Intersection *isect) const;

private:
    BVHBuildNode *recursiveBuild(
        std::vector<BVHPrimitiveInfo> &primitiveInfo,
        int start, int end, int *totalNodes,
        std::vector<std::shared_ptr<Primitive>> &orderedPrims);

    BVHBuildNode *buildUpperSAH(std::vector<BVHBuildNode *> &treeletRoots,
                                int start, int end, int *totalNodes) const;

    int flattenBVHTree(BVHBuildNode *node, int *offset);


    void deleteBuildNode(BVHBuildNode *root);


    //Members
    const int maxPrimsInNode;
    std::vector<std::shared_ptr<Primitive>> primitives;
    LinearBVHNode *nodes = nullptr;


    std::vector<BVHBuildNode> buildNodesVector;
};




//#ifndef PBRT_L1_CACHE_LINE_SIZE
//#define PBRT_L1_CACHE_LINE_SIZE 64
//#endif


//#define PBRT_IS_WINDOWS

//void *AllocAligned(size_t size);

//template <typename T>
//T *AllocAligned(size_t count) {
//    return (T *)AllocAligned(count * sizeof(T));
//}

//void FreeAligned(void *);



//class MemoryArena {
//public:
//    // MemoryArena Public Methods
//    MemoryArena(size_t blockSize = 262144) : blockSize(blockSize) {}
//    ~MemoryArena() {
//        FreeAligned(currentBlock);
//        for (auto &block : usedBlocks) FreeAligned(block.second);
//        for (auto &block : availableBlocks) FreeAligned(block.second);
//    }
//    void *Alloc(size_t nBytes) {
//        // Round up _nBytes_ to minimum machine alignment
//        nBytes = ((nBytes + 15) & (~15));
//        if (currentBlockPos + nBytes > currentAllocSize) {
//            // Add current block to _usedBlocks_ list
//            if (currentBlock) {
//                usedBlocks.push_back(
//                    std::make_pair(currentAllocSize, currentBlock));
//                currentBlock = nullptr;
//                currentAllocSize = 0;
//            }

//            // Get new block of memory for _MemoryArena_

//            // Try to get memory block from _availableBlocks_
//            for (auto iter = availableBlocks.begin();
//                 iter != availableBlocks.end(); ++iter) {
//                if (iter->first >= nBytes) {
//                    currentAllocSize = iter->first;
//                    currentBlock = iter->second;
//                    availableBlocks.erase(iter);
//                    break;
//                }
//            }
//            if (!currentBlock) {
//                currentAllocSize = std::max(nBytes, blockSize);
//                currentBlock = AllocAligned<uint8_t>(currentAllocSize);
//            }
//            currentBlockPos = 0;
//        }
//        void *ret = currentBlock + currentBlockPos;
//        currentBlockPos += nBytes;
//        return ret;
//    }
//    template <typename T>
//    T *Alloc(size_t n = 1, bool runConstructor = true) {
//        T *ret = (T *)Alloc(n * sizeof(T));
//        if (runConstructor)
//            for (size_t i = 0; i < n; ++i) new (&ret[i]) T();
//        return ret;
//    }
//    void Reset() {
//        currentBlockPos = 0;
//        availableBlocks.splice(availableBlocks.begin(), usedBlocks);
//    }
//    size_t TotalAllocated() const {
//        size_t total = currentAllocSize;
//        for (const auto &alloc : usedBlocks) total += alloc.first;
//        for (const auto &alloc : availableBlocks) total += alloc.first;
//        return total;
//    }

//private:
//    MemoryArena(const MemoryArena &) = delete;
//    MemoryArena &operator=(const MemoryArena &) = delete;
//    // MemoryArena Private Data
//    const size_t blockSize;
//    size_t currentBlockPos = 0, currentAllocSize = 0;
//    uint8_t *currentBlock = nullptr;
//    std::list<std::pair<size_t, uint8_t *>> usedBlocks, availableBlocks;
//};
