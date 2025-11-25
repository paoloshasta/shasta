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
#include "markerAccessFunctions.hpp"

// Standard library
#include <algorithm>
#include <array>
#include <cmath>
#include <unordered_set>
#include <unordered_map>
#include <map>
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
    const MemoryMapped::Vector<uint8_t>& variantClusteringValidClustersCompatible,
    const MemoryMapped::Vector<uint8_t>& variantClusteringMemberStatus,
    uint64_t minClusterCoverage,
    uint64_t minAlleleCoverage,
    double minCommonKmerFraction,
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
    data.variantClusteringValidClustersCompatible = &variantClusteringValidClustersCompatible;
    data.variantClusteringMemberStatus = &variantClusteringMemberStatus;
    data.minClusterCoverage = minClusterCoverage;
    data.minAlleleCoverage = minAlleleCoverage;
    data.minCommonKmerFraction = minCommonKmerFraction;
    data.anchorsSkippedNoCommonKmer = 0;
        
    // Initialize the lock-free marker tracking array
    // Allocate it to the total number of markers across all oriented reads
    data.markerCount = markers.totalSize();
    data.markerUsed.reset(new std::atomic<bool>[data.markerCount]);
    for(uint64_t i = 0; i < data.markerCount; i++) {
        data.markerUsed[i].store(false, std::memory_order_relaxed);
    }
    data.anchorsSkippedDuplicateMarkers = 0;

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
    
    // Report how many anchors couldn't be generated due to no common kmer
    const uint64_t skippedNoKmer = data.anchorsSkippedNoCommonKmer;
    if(skippedNoKmer > 0) {
        cout << timestamp << "Warning: " << skippedNoKmer 
             << " anchor(s) could not be generated (no common kmer found)." << endl;
    }
    
    // Report how many anchors were skipped due to duplicate markers
    const uint64_t skippedDuplicates = data.anchorsSkippedDuplicateMarkers;
    if(skippedDuplicates > 0) {
        cout << timestamp << "Info: " << skippedDuplicates 
             << " anchor(s) skipped (markers already used by other anchors)." << endl;
    }
    
    cout << timestamp << "Anchor creation from het sites ends." << endl;
}




// ============================================================================
// Helper functions for anchor generation from het sites
// ============================================================================

namespace {

    // Helper function to get KmerId for a marker (inline for performance)
    inline KmerId getMarkerKmerId(
        const MarkerKmers::MarkerInfo& markerInfo,
        uint64_t k,
        const Reads& reads,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers)
    {
        return getOrientedReadMarkerKmerId(
            markerInfo.orientedReadId,
            markerInfo.ordinal,
            k, reads, markers);
    }
    
    // Helper to get marker at a specific offset from a reference marker
    inline bool getMarkerAtOffset(
        const MarkerKmers::MarkerInfo& refMarker,
        int offset,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        MarkerKmers::MarkerInfo& result)
    {
        const auto orientedReadMarkers = markers[refMarker.orientedReadId.getValue()];
        const int32_t targetOrdinal = static_cast<int32_t>(refMarker.ordinal) + offset;
        
        if(targetOrdinal < 0 || targetOrdinal >= static_cast<int32_t>(orientedReadMarkers.size())) {
            return false;
        }
        
        result.orientedReadId = refMarker.orientedReadId;
        result.ordinal = static_cast<uint32_t>(targetOrdinal);
        return true;
    }
    
    // Cache structure for pre-extracted kmers
    struct ReadKmers {
        std::array<KmerId, 6> kmers;      // prev-2, prev-1, prev, next, next+1, next+2
        std::array<MarkerKmers::MarkerInfo, 6> markers;
        uint8_t validMask;                 // Bitmask: bit i set if position i is valid
        
        ReadKmers() : validMask(0) {}
    };
    
    // Step 1: Check if enough reads (>=minFraction) have same kmer at a fixed position
    inline bool tryFixedPositionCached(
        const std::vector<ReadKmers>& readKmersCache,
        const size_t cacheIndex,
        std::vector<MarkerKmers::MarkerInfo>& markerInfos,
        const double minFraction = 0.8)  // Default 80%
    {
        const uint8_t requiredBit = (1 << cacheIndex);
        const size_t minReadsRequired = static_cast<size_t>(std::ceil(readKmersCache.size() * minFraction));
        
        // Count kmers that appear at this position
        std::map<KmerId, std::vector<size_t>> kmerToReadIndices;
        
        for(size_t i = 0; i < readKmersCache.size(); i++) {
            const auto& rk = readKmersCache[i];
            if(rk.validMask & requiredBit) {
                kmerToReadIndices[rk.kmers[cacheIndex]].push_back(i);
            }
        }
        
        // Find the most common kmer that meets the threshold
        for(const auto& [kmerId, readIndices] : kmerToReadIndices) {
            if(readIndices.size() >= minReadsRequired) {
                // Found a kmer shared by enough reads!
                markerInfos.clear();
                markerInfos.reserve(readIndices.size());
                for(size_t readIdx : readIndices) {
                    markerInfos.push_back(readKmersCache[readIdx].markers[cacheIndex]);
                }
                return true;
            }
        }
        
        return false;
    }
    
    // Step 2/3: Find common kmer across multiple positions (all vs all)
    // Now supports partial matches: kmer must appear in >= minFraction of reads
    inline bool tryMultipleOffsetsCached(
        const std::vector<ReadKmers>& readKmersCache,
        const std::vector<size_t>& cacheIndices,
        std::vector<MarkerKmers::MarkerInfo>& markerInfos,
        const double minFraction = 0.8)  // Default 80%
    {
        const size_t minReadsRequired = static_cast<size_t>(std::ceil(readKmersCache.size() * minFraction));
        
        // Map: KmerId -> [(readIdx, cacheIdx), ...]
        std::map<KmerId, std::vector<std::pair<size_t, size_t>>> kmerToReads;
        
        // Collect all kmers from all reads at specified positions
        for(size_t readIdx = 0; readIdx < readKmersCache.size(); readIdx++) {
            const ReadKmers& rk = readKmersCache[readIdx];
            
            for(size_t cacheIdx : cacheIndices) {
                if(rk.validMask & (1 << cacheIdx)) {
                    kmerToReads[rk.kmers[cacheIdx]].push_back({readIdx, cacheIdx});
                }
            }
        }
        
        // Find a kmer that appears in >= minReadsRequired reads (once per read)
        for(const auto& [kmerId, readPairs] : kmerToReads) {
            // Check if this kmer appears in enough reads
            if(readPairs.size() < minReadsRequired) {
                continue;
            }
            
            // Verify each read appears exactly once (no duplicates from same read)
            std::vector<uint8_t> readSeen(readKmersCache.size(), 0);
            std::vector<std::pair<size_t, size_t>> uniqueReadPairs;
            uniqueReadPairs.reserve(readPairs.size());
            
            for(const auto& [readIdx, cacheIdx] : readPairs) {
                if(!readSeen[readIdx]) {
                    readSeen[readIdx] = 1;
                    uniqueReadPairs.push_back({readIdx, cacheIdx});
                }
            }
            
            // Check if we still have enough unique reads after deduplication
            if(uniqueReadPairs.size() < minReadsRequired) {
                continue;
            }
            
            // Found a common kmer shared by enough reads!
            markerInfos.clear();
            markerInfos.reserve(uniqueReadPairs.size());
            for(const auto& [readIdx, cacheIdx] : uniqueReadPairs) {
                markerInfos.push_back(readKmersCache[readIdx].markers[cacheIdx]);
            }
            return true;
        }
        
        return false;
    }
    
    // Helper to compute MarkerId from MarkerInfo (OrientedReadId + ordinal)
    inline MarkerId getMarkerId(
        const MarkerKmers::MarkerInfo& markerInfo,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers)
    {
        return (markers.begin(markerInfo.orientedReadId.getValue()) - markers.begin()) 
               + markerInfo.ordinal;
    }
    
    // Lock-free helper function to atomically check and claim markers
    // Returns true if all markers were available and have been marked as used
    // Returns false if any marker was already used
    // Uses compare_exchange for lock-free atomic operations
    inline bool tryClaimMarkers(
        const std::vector<MarkerKmers::MarkerInfo>& markerInfos,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        std::atomic<bool>* markerUsed)
    {
        // Try to atomically claim each marker
        std::vector<MarkerId> claimedMarkers;
        claimedMarkers.reserve(markerInfos.size());
        
        for(const auto& markerInfo : markerInfos) {
            const MarkerId markerId = getMarkerId(markerInfo, markers);
            
            // Atomically try to change from false (available) to true (used)
            bool expected = false;
            if(markerUsed[markerId].compare_exchange_strong(expected, true, 
                                                            std::memory_order_acquire,
                                                            std::memory_order_relaxed)) {
                // Successfully claimed this marker
                claimedMarkers.push_back(markerId);
            } else {
                // This marker was already used - need to release the ones we claimed
                for(MarkerId releasedId : claimedMarkers) {
                    markerUsed[releasedId].store(false, std::memory_order_release);
                }
                return false;
            }
        }
        
        // All markers successfully claimed
        return true;
    }
    
}  // anonymous namespace


void Anchors::constructFromHetSitesThreadFunction(uint64_t threadId)
{

    ofstream debugOut("Debug-" + to_string(threadId) + ".txt");
    debugOut << "Thread ID: " << threadId << " is starting anchor generation from het sites." << endl;

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
    uint64_t clustersSkippedNotCompatible = 0;
    
    while(getNextBatch(begin, end)) {

        // Loop over clusters assigned to this thread.
        for(uint64_t clusterIdx = begin; clusterIdx < end; clusterIdx++) {

            // Get the actual cluster ID (representative ID)
            const uint64_t clusterId = (*data.clusterRepresentatives)[clusterIdx];

            // FILTER: Only process clusters that were marked compatible in ReadGraph5
            // Use clusterId, NOT clusterIdx
            if (!(*data.variantClusteringValidClustersCompatible)[clusterId]) {
                clustersSkipped++;
                clustersSkippedNotCompatible++;
                continue;
            }
            
            clustersProcessed++;
            
            // Get the member ids of this cluster
            // This is a list of the position pair ids that are members of this cluster.
            std::vector<uint64_t> memberIds = data.membersByRepIdx[clusterIdx];
            if(memberIds.empty()) {
                continue;
            }

            // If cluster coverage is less than minClusterCoverage, skip it.
            if(memberIds.size() < minClusterCoverage) {
                clustersSkipped++;
                continue;
            }

            // Sort memberIds by their position pairs to enable efficient duplicate checking
            std::sort(memberIds.begin(), memberIds.end(),
            [&positionPairs](uint64_t a, uint64_t b) {
                return positionPairs[a] < positionPairs[b];
            });

            
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


            // Check if there are at least 2 alleles with >= minAlleleCoverage coverage
            uint64_t eligibleAlleles = 0;
            for(uint8_t a=0; a<5; a++) {
                if(readsByAllele[a].size() >= minAlleleCoverage) {
                    eligibleAlleles++;
                }
            }
            if(eligibleAlleles < 2) {
                clustersSkipped++;
                continue;
            }



            // Check for duplicate OrientedReadIds and same ReadId on both strands
            // Since the memberIds are sorted, we only need to check consecutive pairs
            bool hasDuplicates = false;
            bool hasSameReadIdBothStrands = false;

            if(memberIds.size() >= 2) {
                for(uint64_t i = 1; i < memberIds.size(); i++) {
                    const auto& pp0 = positionPairs[memberIds[i-1]];
                    const auto& pp1 = positionPairs[memberIds[i]];
                    
                    const OrientedReadId orientedReadId0 = pp0.first;
                    const OrientedReadId orientedReadId1 = pp1.first;
                    
                    // Check for duplicate OrientedReadIds
                    if(orientedReadId0 == orientedReadId1) {
                        hasDuplicates = true;
                        break;
                    }
                    
                    // Check if same ReadId appears on both strands
                    if(orientedReadId0.getReadId() == orientedReadId1.getReadId()) {
                        hasSameReadIdBothStrands = true;
                        break;
                    }
                }
            }

            if(hasDuplicates) {
                clustersSkipped++;
                clustersSkippedDuplicate++;
                continue;
            }

            if(hasSameReadIdBothStrands) {
                clustersSkipped++;
                continue;
            }



            // Find the lowest number orientedReadId.
            // If the lowest orientedReadId is on strand 1, skip it.
            // Since for each genomic position we have created, by design, two
            // clusters (one set of OrientedReadIds and one set of OrientedReadIds on their opposite strand),
            // we want to create anchors from the processing of one of them.
            // Since memberIds is sorted, the first element has the lowest OrientedReadId.
            const auto& firstPair = positionPairs[memberIds.front()];
            const OrientedReadId lowestOrientedReadId = firstPair.first;

            // Check if lowest is on strand 1
            if(lowestOrientedReadId.getStrand() == 1) {
                clustersSkipped++;
                clustersSkippedStrand1Lowest++;
                continue;
            }


            // If getting here, we will generate a pair of Anchors corresponding
            // to this cluster, for each allele that has >= minAlleleCoverage coverage.

            // ============================================================================
            // MAIN ANCHOR GENERATION LOOP
            // ============================================================================

            // Thread-local storage to avoid repeated allocations
            std::vector<MarkerInfo> anchorMarkerInfos;
            std::vector<uint64_t> alleleMemberIds;

            // For each allele that has >= minAlleleCoverage coverage
            for(uint8_t allele = 0; allele < 5; allele++) {

                // If the allele has less than minAlleleCoverage coverage, skip it.
                if(readsByAllele[allele].size() < minAlleleCoverage) {
                    continue;
                }
                
                // Collect the memberIds for this allele.
                // This is a list of the position pair ids that are members of this allele.
                alleleMemberIds.clear();
                alleleMemberIds.reserve(readsByAllele[allele].size());
                
                for(uint64_t id: memberIds) {
                    // Skip members marked as stray by graph refinement
                    if(id < data.variantClusteringMemberStatus->size() && 
                       (*data.variantClusteringMemberStatus)[id] == 1) {
                        continue;  // Skip stray/filtered reads
                    }
                    
                    if(id < positionPairAlleles.size() && positionPairAlleles[id] == allele) {
                        alleleMemberIds.push_back(id);
                    }
                }
                
                if(alleleMemberIds.empty()) continue;
                
                // ========================================================================
                // OPTIMIZATION: Pre-extract all kmers once
                // ========================================================================
                std::vector<ReadKmers> readKmersCache(alleleMemberIds.size());
                
                for(size_t readIdx = 0; readIdx < alleleMemberIds.size(); readIdx++) {
                    const auto& context = positionPairContexts[alleleMemberIds[readIdx]];
                    ReadKmers& rk = readKmersCache[readIdx];
                    
                    // Extract prev-2, prev-1, prev (cache indices 0, 1, 2)
                    MarkerInfo marker;
                    for(int offset = -2; offset <= 0; offset++) {
                        if(getMarkerAtOffset(context.prevMarkerInfo, offset, markers, marker)) {
                            const size_t idx = offset + 2;
                            rk.markers[idx] = marker;
                            rk.kmers[idx] = getMarkerKmerId(marker, k, reads, markers);
                            rk.validMask |= (1 << idx);
                        }
                    }
                    
                    // Extract next, next+1, next+2 (cache indices 3, 4, 5)
                    for(int offset = 0; offset <= 2; offset++) {
                        if(getMarkerAtOffset(context.nextMarkerInfo, offset, markers, marker)) {
                            const size_t idx = offset + 3;
                            // Index 0: prev-2 marker
                            // Index 1: prev-1 marker
                            // Index 2: prev marker (the one before the variant)
                            // Index 3: next marker (the one after the variant)
                            // Index 4: next+1 marker
                            // Index 5: next+2 marker
                            rk.markers[idx] = marker;
                            rk.kmers[idx] = getMarkerKmerId(marker, k, reads, markers);
                            rk.validMask |= (1 << idx);
                        }
                    }

                    // Example:
                    // readKmersCache[0] = {
                    //     kmers:     [123456, 234567, 345678, 456789, 567890, 678901],
                    //     markers:   [{read5-0,13}, {read5-0,14}, {read5-0,15}, {read5-0,16}, {read5-0,17}, {read5-0,18}],
                    //     validMask: 0b111111  // All 6 positions are valid
                    // }
                }
                
                // ========================================================================
                // Progressive search for common kmer
                // ========================================================================
                bool foundCommonKmer = false;
                
                // Step 1: Try prev (cache index 2)
                if(tryFixedPositionCached(readKmersCache, 2, anchorMarkerInfos, data.minCommonKmerFraction)) {
                    foundCommonKmer = true;
                }
                // Step 1b: Try next (cache index 3)
                else if(tryFixedPositionCached(readKmersCache, 3, anchorMarkerInfos, data.minCommonKmerFraction)) {
                    foundCommonKmer = true;
                }
                // Step 2: Try prev and prev-1 (cache indices 2, 1)
                else if(tryMultipleOffsetsCached(readKmersCache, {2, 1}, anchorMarkerInfos, data.minCommonKmerFraction)) {
                    foundCommonKmer = true;
                }
                // Step 2b: Try next and next+1 (cache indices 3, 4)
                else if(tryMultipleOffsetsCached(readKmersCache, {3, 4}, anchorMarkerInfos, data.minCommonKmerFraction)) {
                    foundCommonKmer = true;
                }
                // Step 3: Try prev, prev-1, prev-2 (cache indices 2, 1, 0)
                else if(tryMultipleOffsetsCached(readKmersCache, {2, 1, 0}, anchorMarkerInfos, data.minCommonKmerFraction)) {
                    foundCommonKmer = true;
                }
                // Step 3b: Try next, next+1, next+2 (cache indices 3, 4, 5)
                else if(tryMultipleOffsetsCached(readKmersCache, {3, 4, 5}, anchorMarkerInfos, data.minCommonKmerFraction)) {
                    foundCommonKmer = true;
                }
                
                // Step 4: If no common kmer found, print warning and skip
                if(!foundCommonKmer) {
                    // Get cluster representative info for debugging
                    const uint64_t representativeId = (*(data.clusterRepresentatives))[clusterIdx];
                    const auto& repPositionPair = positionPairs[representativeId];
                    const OrientedReadId repReadId = repPositionPair.first;
                    const uint32_t repPosition = repPositionPair.second;
                    
                    debugOut << "\n========================================" << endl;
                    debugOut << "Warning: Could not find common kmer for allele " << int(allele) << endl;
                    debugOut << "Cluster " << clusterIdx 
                        << " - Representative: " << repReadId << " at position " << repPosition << endl;
                    debugOut << "Total cluster size: " << memberIds.size() << " reads" << endl;
                    debugOut << "Allele " << int(allele) << " size: " << alleleMemberIds.size() << " reads" << endl;
                    
                    // Print allele distribution summary
                    debugOut << "Allele coverage: ";
                    for(uint8_t a = 0; a < 5; a++) {
                        if(readsByAllele[a].size() > 0) {
                            debugOut << "A" << int(a) << "=" << readsByAllele[a].size() << " ";
                        }
                    }
                    debugOut << endl;
                    
                    // Print full cluster details
                    debugOut << "Full cluster details:" << endl;
                    for(uint64_t id : memberIds) {
                        const auto& pp = positionPairs[id];
                        const OrientedReadId rid = pp.first;
                        const uint32_t pos = pp.second;
                        
                        debugOut << "  Member " << id << ": " << rid << " pos:" << pos;
                        
                        // Print allele if available
                        if(id < positionPairAlleles.size()) {
                            const uint8_t a = positionPairAlleles[id];
                            debugOut << " allele:" << int(a);
                            if(a == allele) {
                                debugOut << " <-- FAILED ALLELE";
                            }
                        } else {
                            debugOut << " allele:N/A";
                        }
                        debugOut << endl;
                    }
                    
                    debugOut << "Skipping anchor generation for this allele." << endl;
                    debugOut << "========================================\n" << endl;
                    
                    data.anchorsSkippedNoCommonKmer++;
                    continue;
                }
                
                // ========================================================================
                // Generate anchor pair (forward + reverse complement)
                // ========================================================================
                
                // Log if we're using a subset of reads (less than 100%)
                if(anchorMarkerInfos.size() < alleleMemberIds.size()) {
                    debugOut << "Note: Found common kmer for " << anchorMarkerInfos.size() 
                        << " out of " << alleleMemberIds.size() << " reads (" 
                        << (100.0 * anchorMarkerInfos.size() / alleleMemberIds.size()) 
                        << "%) in allele " << int(allele) << " of cluster " << clusterIdx << endl;
                }
                
                // Prepare reverse complement marker infos first
                std::vector<MarkerInfo> rcAnchorMarkerInfos = anchorMarkerInfos;
                for(auto& markerInfo: rcAnchorMarkerInfos) {
                    markerInfo.orientedReadId.flipStrand();
                    const uint64_t markerCount = markers.size(markerInfo.orientedReadId.getValue());
                    markerInfo.ordinal = uint32_t(markerCount) - 1 - markerInfo.ordinal;
                }
                
                // Try to claim BOTH forward and RC markers (all-or-nothing)
                bool forwardClaimed = tryClaimMarkers(anchorMarkerInfos, markers, data.markerUsed.get());
                bool rcClaimed = false;
                
                if(forwardClaimed) {
                    // Forward succeeded, now try RC
                    rcClaimed = tryClaimMarkers(rcAnchorMarkerInfos, markers, data.markerUsed.get());
                    
                    if(!rcClaimed) {
                        // RC failed - rollback the forward markers we just claimed
                        for(const auto& markerInfo : anchorMarkerInfos) {
                            const MarkerId markerId = getMarkerId(markerInfo, markers);
                            data.markerUsed[markerId].store(false, std::memory_order_release);
                        }
                    }
                }
                
                // Only create anchors if BOTH forward and RC were claimed successfully
                if(forwardClaimed && rcClaimed) {
                    // Both succeeded - create the anchor pair
                    threadAnchors.appendVector(anchorMarkerInfos);
                    threadAnchors.appendVector(rcAnchorMarkerInfos);
                } else {
                    // At least one failed - skip this anchor entirely
                    debugOut << "Skipping anchor for allele " << int(allele) 
                        << " in cluster " << clusterIdx 
                        << " - markers already used by another anchor (forward=" 
                        << forwardClaimed << ", rc=" << rcClaimed << ")." << endl;
                    data.anchorsSkippedDuplicateMarkers++;
                    continue;
                }
            }
            





        }
    }
    
    // DEBUG: Print summary for this thread (Uncomment to see results)

    debugOut << "Thread " << threadId << " summary:" << endl;
    debugOut << "  Clusters processed: " << clustersProcessed << endl;
    debugOut << "  Clusters skipped (not compatible): " << clustersSkippedNotCompatible << endl;
    debugOut << "  Clusters skipped (duplicates): " << clustersSkippedDuplicate << endl;
    debugOut << "  Clusters skipped (strand-1 lowest): " << clustersSkippedStrand1Lowest << endl;
    debugOut << "  Total clusters skipped: " << clustersSkipped << endl;

}