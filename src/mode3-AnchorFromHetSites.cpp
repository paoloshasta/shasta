// Shasta
#include "mode3-Anchor.hpp"
#include "dset64-gccAtomic.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
#include "Marker.hpp"
#include "MarkerKmers.hpp"
#include "Reads.hpp"
#include "LongBaseSequence.hpp"
#include "extractKmer.hpp"

// Standard library
#include <algorithm>
#include <array>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <limits>

using namespace shasta;
using namespace mode3;


// Constructor for creating Anchors from heterozygous sites.
Anchors::Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const std::vector<uint64_t>& clusterRepresentatives,
    DisjointSets& disjointSets,
    const MemoryMapped::Vector< std::pair<OrientedReadId, uint32_t> >& positionPairs,
    const MemoryMapped::Vector<uint8_t>& positionPairAlleles,
    const MemoryMapped::Vector<VariantPositionContext>& positionPairContexts,
    uint64_t minClusterCoverage,
    uint64_t minAlleleCoverage,
    uint64_t threadCount) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    k(k),
    markers(markers)
{
    kHalf = k / 2;

    cout << timestamp << "Anchor creation from het sites begins over " << clusterRepresentatives.size() << " clusters." << endl;

    // Store arguments so all threads can see them.
    ConstructFromHetSitesData& data = constructFromHetSitesData;
    data.clusterRepresentatives = &clusterRepresentatives;
    data.disjointSets = &disjointSets;
    data.positionPairs = &positionPairs;
    data.positionPairAlleles = &positionPairAlleles;
    data.positionPairContexts = &positionPairContexts;
    data.minClusterCoverage = minClusterCoverage;
    data.minAlleleCoverage = minAlleleCoverage;

    // Create the membersByRepIdx structure, a data structure that maps
    // each cluster (by its index) to all its member position pair IDs.

    // Pre-bucket position pair ids by their representative
    const uint64_t N = data.positionPairs->size();

    // Build sorted lookup array: O(C log C)
    std::vector<std::pair<uint64_t, size_t>> sortedReps;
    sortedReps.reserve(clusterRepresentatives.size());
    for(size_t i = 0; i < clusterRepresentatives.size(); i++) {
        sortedReps.push_back({clusterRepresentatives[i], i});
    }
    std::sort(sortedReps.begin(), sortedReps.end());
    
    // Allocate result storage
    data.membersByRepIdx.resize(clusterRepresentatives.size());
    
    // Single pass with binary search: O(N log C)
    for(uint64_t id = 0; id < N; id++) {
        const uint64_t rep = disjointSets.find(id);
        
        auto it = std::lower_bound(sortedReps.begin(), sortedReps.end(), 
                                    std::make_pair(rep, size_t(0)));
        
        if(it != sortedReps.end() && it->first == rep) {
            data.membersByRepIdx[it->second].push_back(id);
        }
    }

    // Each thread stores the anchors it finds in a separate vector.
    data.threadAnchors.clear();
    data.threadAnchors.resize(threadCount);
    const uint64_t batchSize = 1;
    setupLoadBalancing(clusterRepresentatives.size(), batchSize);
    runThreads(&Anchors::constructFromHetSitesThreadFunction, threadCount);

    // Gather the anchors found by all threads.
    performanceLog << timestamp << "Gathering the anchors found by all threads." << endl;
    anchorMarkerIntervals.createNew(
            largeDataName("AnchorMarkerIntervals"),
            largeDataPageSize);
    anchorSequences.createNew(
        largeDataName("AnchorSequences"), largeDataPageSize);
    anchorInfos.createNew(largeDataName("AnchorInfos"), largeDataPageSize);
    for(uint64_t threadId=0; threadId<threadCount; threadId++) {
        auto& threadAnchorsPointer = data.threadAnchors[threadId];
        auto& threadAnchors = *threadAnchorsPointer;

        // Loop over the anchors found by this thread.
        for(uint64_t i=0; i<threadAnchors.size(); i++) {
            const auto threadAnchor = threadAnchors[i];
            anchorMarkerIntervals.appendVector();
            for(const auto& markerInfo: threadAnchor) {
                anchorMarkerIntervals.append(AnchorMarkerInterval(markerInfo.orientedReadId, markerInfo.ordinal));
            }
        }
        threadAnchors.remove();
        threadAnchorsPointer = 0;
    }

    // The anchor sequences are all empty.
    // The ordinal offsets are all 0 (ordinal0 == ordinal1 since we store a single marker ordinal).
    const uint64_t anchorCount = anchorMarkerIntervals.size();
    anchorInfos.resize(anchorCount);
    for(AnchorId anchorId=0; anchorId<anchorCount; anchorId++) {
        anchorSequences.appendVector();
        anchorInfos[anchorId].ordinalOffset = 0;
    }

    cout << timestamp << "Generated " << anchorCount << " anchors from het sites." << endl;
    cout << timestamp << "Anchor creation from het sites ends." << endl;
}



void Anchors::constructFromHetSitesThreadFunction(uint64_t threadId)
{
    using MarkerInfo = MarkerKmers::MarkerInfo;

    ConstructFromHetSitesData& data = constructFromHetSitesData;
    const uint64_t minClusterCoverage = data.minClusterCoverage;
    const uint64_t minAlleleCoverage = data.minAlleleCoverage;

    const MemoryMapped::Vector< std::pair<OrientedReadId, uint32_t> >& positionPairs = *(data.positionPairs);
    const MemoryMapped::Vector<uint8_t>& positionPairAlleles = *(data.positionPairAlleles);
    const MemoryMapped::Vector<VariantPositionContext>& positionPairContexts = *(data.positionPairContexts);
    

    // Initialize the anchors that will be found by this thread
    auto& threadAnchorsPointer = data.threadAnchors[threadId];
    threadAnchorsPointer = make_shared<MemoryMapped::VectorOfVectors<ConstructFromHetSitesData::MarkerInfo, uint64_t> >();
    auto& threadAnchors = *threadAnchorsPointer;
    threadAnchors.createNew(largeDataName("tmp-threadAnchors-" + to_string(threadId)), largeDataPageSize);

    // Process batches of clusters (each thread processes different clusters)
    uint64_t begin, end;
    uint64_t clustersProcessed = 0;
    uint64_t clustersSkipped = 0;
    uint64_t clustersSkippedDuplicate = 0;
    uint64_t clustersSkippedStrand1Lowest = 0;
    
    while(getNextBatch(begin, end)) {
        // Loop over clusters assigned to this thread.
        for(uint64_t clusterIdx = begin; clusterIdx < end; clusterIdx++) {
            clustersProcessed++;
            
            // Get the member ids of this cluster
            const std::vector<uint64_t>& memberIds = data.membersByRepIdx[clusterIdx];
            if(memberIds.empty()) {
                continue;
            }

            // Gather reads and allele support for this cluster
            std::unordered_set<uint32_t> allReads;
            std::array<std::unordered_set<uint32_t>, 5> readsByAllele{};

            for(uint64_t id: memberIds) {
                const auto& pp = positionPairs[id];
                const OrientedReadId rid = pp.first;
                const uint32_t ridVal = rid.getValue();
                allReads.insert(ridVal);

                if(id < positionPairAlleles.size()) {
                    const uint8_t a = positionPairAlleles[id];
                    if(a < 5) {
                        readsByAllele[a].insert(ridVal);
                    }
                }
            }

            // Eligibility check: at least 2 alleles with >= minReadsPerAllele coverage
            if(allReads.size() < minClusterCoverage) {
                continue;
            }
            uint64_t eligibleAlleles = 0;
            for(uint8_t a=0; a<5; a++) {
                if(readsByAllele[a].size() >= minAlleleCoverage) {
                    eligibleAlleles++;
                }
            }
            if(eligibleAlleles < 2) {
                continue;
            }

            // Step 1: Check for duplicate orientedReadIds
            std::set<uint32_t> uniqueOrientedReadIds;
            bool hasDuplicates = false;
            for(uint64_t id: memberIds) {
                const auto& pp = positionPairs[id];
                const uint32_t orientedReadIdValue = pp.first.getValue();
                if(uniqueOrientedReadIds.count(orientedReadIdValue)) {
                    hasDuplicates = true;
                    break;
                }
                uniqueOrientedReadIds.insert(orientedReadIdValue);
            }
            
            if(hasDuplicates) {
                clustersSkipped++;
                clustersSkippedDuplicate++;
                continue;
            }

            // Step 2: Find the lowest number orientedReadId
            uint32_t lowestOrientedReadIdValue = std::numeric_limits<uint32_t>::max();
            OrientedReadId lowestOrientedReadId;
            for(uint64_t id: memberIds) {
                const auto& pp = positionPairs[id];
                const uint32_t orientedReadIdValue = pp.first.getValue();
                if(orientedReadIdValue < lowestOrientedReadIdValue) {
                    lowestOrientedReadIdValue = orientedReadIdValue;
                    lowestOrientedReadId = pp.first;
                }
            }
            
            // Check if lowest is on strand 1
            if(lowestOrientedReadId.getStrand() == 1) {
                clustersSkipped++;
                clustersSkippedStrand1Lowest++;
                continue;
            }

            // Step 3: If lowest is on strand 0, do nothing for now
            // (Future processing will go here)
            
        }
    }
    
    // // DEBUG: Print summary for this thread
    // cout << "Thread " << threadId << " summary:" << endl;
    // cout << "  Clusters processed: " << clustersProcessed << endl;
    // cout << "  Clusters skipped (duplicates): " << clustersSkippedDuplicate << endl;
    // cout << "  Clusters skipped (strand-1 lowest): " << clustersSkippedStrand1Lowest << endl;
    // cout << "  Total clusters skipped: " << clustersSkipped << endl;
}