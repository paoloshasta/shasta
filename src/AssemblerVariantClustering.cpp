// Shasta.
#include "Assembler.hpp"
#include "ProjectedAlignment.hpp"
#include "compressAlignment.hpp"
#include "dset64-gccAtomic.hpp"
#include "SHASTA_ASSERT.hpp"
#include "iostream.hpp"
#include "Reads.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
#include "chrono.hpp"

#include "MarkerKmers.hpp"


using namespace shasta;

// Helper function called during alignment computation in AssemblerAlign.cpp to collect position pairs
// that have SNP differences
void Assembler::collectVariantClusteringPositionPairs(
    const ProjectedAlignment& projectedAlignment,
    const array<OrientedReadId, 2>& orientedReadIds,
    MemoryMapped::Vector< pair<OrientedReadId, uint32_t> >& positionPairs)
{
    const OrientedReadId currentOrientedReadId0 = orientedReadIds[0];
    const OrientedReadId currentOrientedReadId1 = orientedReadIds[1];
    
    // Process each segment to find differences
    for (const ProjectedAlignmentSegment& segment : projectedAlignment.segments) {
        
        const vector<Base>& sequence0 = segment.sequences[0];
        const vector<Base>& sequence1 = segment.sequences[1];
        
        // Track positions in each sequence (no gaps)
        uint64_t position0 = 0;
        uint64_t position1 = 0;
        
        // Iterate through the alignment directly
        for (const auto& basePair : segment.alignment) {
            const bool hasBase0 = basePair.first;
            const bool hasBase1 = basePair.second;
            
            // Check for mismatch (both are bases, not gaps, and different)
            if (hasBase0 && hasBase1) {
                // Both reads have a base at this alignment position
                SHASTA_ASSERT(position0 < sequence0.size());
                SHASTA_ASSERT(position1 < sequence1.size());
                
                const Base& base0 = sequence0[position0];
                const Base& base1 = sequence1[position1];

                if (base0 != base1) {
                    // Mismatch found!
                    uint64_t positionInRead0 = segment.positionsA[0] + position0;
                    uint64_t positionInRead1 = segment.positionsA[1] + position1;

                    positionPairs.push_back(make_pair(currentOrientedReadId0, uint32_t(positionInRead0)));
                    positionPairs.push_back(make_pair(currentOrientedReadId1, uint32_t(positionInRead1)));

                    // Create reverse complement position pairs
                    OrientedReadId currentOrientedReadId0Rc(currentOrientedReadId0.getReadId(), 1 - currentOrientedReadId0.getStrand());
                    OrientedReadId currentOrientedReadId1Rc(currentOrientedReadId1.getReadId(), 1 - currentOrientedReadId1.getStrand());
                    const uint64_t readLength0Rc = getReads().getReadRawSequenceLength(currentOrientedReadId0Rc.getReadId());
                    const uint64_t readLength1Rc = getReads().getReadRawSequenceLength(currentOrientedReadId1Rc.getReadId());
                    positionPairs.push_back(make_pair(currentOrientedReadId0Rc, uint32_t(readLength0Rc - 1 - positionInRead0)));
                    positionPairs.push_back(make_pair(currentOrientedReadId1Rc, uint32_t(readLength1Rc - 1 - positionInRead1)));
                }
                
                // Increment position counters for both sequences
                position0++;
                position1++;
                
            } else if (hasBase0) {
                // Only read 0 has a base (read 1 has a gap)
                position0++;
            } else if (hasBase1) {
                // Only read 1 has a base (read 0 has a gap)
                position1++;
            }
            // else: both have gaps (shouldn't happen in a valid alignment, but handle gracefully)
        }
    }
}

void Assembler::accessVariantClusteringPositionPairsReadOnly()
{
    variantClusteringPositionPairs.accessExistingReadOnly(
        largeDataName("VariantClusteringPositionPairs"));
    cout << "Accessed " << variantClusteringPositionPairs.size() 
         << " position pairs for variant clustering." << endl;
}

void Assembler::accessVariantClusteringPositionPairsReadWrite()
{
    variantClusteringPositionPairs.accessExistingReadWrite(
        largeDataName("VariantClusteringPositionPairs"));
    cout << "Accessed " << variantClusteringPositionPairs.size() 
         << " position pairs for variant clustering." << endl;
}

void Assembler::checkVariantClusteringPositionPairsIsOpen() const
{
    if(!variantClusteringPositionPairs.isOpen) {
        throw runtime_error("Variant clustering position pairs are not accessible.");
    }
}



void Assembler::storeVariantClusteringPositionPairs(
    size_t threadCount,
    ComputeAlignmentsData& data)
{
    // Store position pairs collected for variant clustering
    performanceLog << timestamp << "Storing position pairs for variant clustering." << endl;
    
    // First, compute total number of pairs needed
    size_t totalPairs = 0;
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        auto& threadPairsPointer = data.threadVariantClusteringPositionPairs[threadId];
        if(threadPairsPointer) {
            totalPairs += threadPairsPointer->size();
        }
    }

    // Create vector with appropriate capacity
    variantClusteringPositionPairs.createNew(
        largeDataName("VariantClusteringPositionPairs"), largeDataPageSize, 0, totalPairs);

    // Append all pairs from each thread
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        auto& threadPairsPointer = data.threadVariantClusteringPositionPairs[threadId];
        if(threadPairsPointer) {
            auto& threadPairs = *threadPairsPointer;
            // Append all pairs from this thread in bulk
            const size_t oldSize = variantClusteringPositionPairs.size();
            const size_t newSize = oldSize + threadPairs.size();
            variantClusteringPositionPairs.resize(newSize);
            std::copy(threadPairs.begin(), threadPairs.end(), 
                        variantClusteringPositionPairs.begin() + oldSize);
            // Free the thread-local storage immediately
            threadPairs.remove();
        }
    }
    variantClusteringPositionPairs.unreserve();
    performanceLog << timestamp << "Stored " << variantClusteringPositionPairs.size() << " position pair entries for variant clustering." << endl;
    cout << timestamp << "Stored " << variantClusteringPositionPairs.size() << " position pair entries for variant clustering." << endl;
}





// Phase 2 thread function: Link position pairs using disjoint sets
void Assembler::linkVariantClustersThreadFunction(uint64_t threadId)
{
    DisjointSets& disjointSets = *variantClusteringDisjointSets;
    
    // Pre-compute pointers for fast binary search access
    const auto pairsBegin = variantClusteringPositionPairs.begin();
    const auto pairsEnd = variantClusteringPositionPairs.end();
    const uint64_t maxId = variantClusteringPositionPairs.size();

    // Cache frequently accessed members
    const auto& alignmentDataRef = alignmentData;
    const auto& markersRef = markers;
    const auto& compressedAlignmentsRef = compressedAlignments;
        
    // Get batches of alignment IDs to process
    uint64_t alignmentIdBegin, alignmentIdEnd;
    while (getNextBatch(alignmentIdBegin, alignmentIdEnd)) {
        
        // Process alignments directly
        for (uint64_t alignmentId = alignmentIdBegin; alignmentId < alignmentIdEnd; alignmentId++) {
            
            const AlignmentData& alignmentData = alignmentDataRef[alignmentId];
            
            // Get oriented read IDs from alignment data
            OrientedReadId currentOrientedReadId0(alignmentData.readIds[0], 0);
            OrientedReadId currentOrientedReadId1(alignmentData.readIds[1], alignmentData.isSameStrand ? 0 : 1);
            AlignmentInfo alignmentInfo = alignmentData.info;

            // Decompress alignment
            Alignment alignment;
            const span<const char> compressedAlignment = compressedAlignmentsRef[alignmentId];
            shasta::decompress(compressedAlignment, alignment);
                                
            // Project alignment to base space
            const ProjectedAlignment projectedAlignment(
                *this,
                {currentOrientedReadId0, currentOrientedReadId1},
                alignment,
                ProjectedAlignment::Method::QuickRaw);

            // Process each segment to find differences
            for (const ProjectedAlignmentSegment& segment : projectedAlignment.segments) {
                
                const vector<Base>& sequence0 = segment.sequences[0];
                const vector<Base>& sequence1 = segment.sequences[1];
                
                // Track positions in each sequence (no gaps)
                uint64_t position0 = 0;
                uint64_t position1 = 0;
                
                // Iterate through the alignment directly
                for (const auto& basePair : segment.alignment) {
                    const bool hasBase0 = basePair.first;
                    const bool hasBase1 = basePair.second;
                    
                    // Check for mismatch (both are bases, not gaps, and different)
                    if (hasBase0 && hasBase1) {
                        // Both reads have a base at this alignment position
                        SHASTA_ASSERT(position0 < sequence0.size());
                        SHASTA_ASSERT(position1 < sequence1.size());
                        
                        const Base& base0 = sequence0[position0];
                        const Base& base1 = sequence1[position1];

                        if (base0 != base1) {
                            // Mismatch found!
                            uint64_t positionInRead0 = segment.positionsA[0] + position0;
                            uint64_t positionInRead1 = segment.positionsA[1] + position1;

                            // Create position pairs (in their original coordinate system)
                            pair<OrientedReadId, uint32_t> pair0(currentOrientedReadId0, uint32_t(positionInRead0));
                            pair<OrientedReadId, uint32_t> pair1(currentOrientedReadId1, uint32_t(positionInRead1));
                            
                            // Binary search to find IDs in sorted vector (O(log n) lookup)
                            auto it0 = std::lower_bound(pairsBegin, pairsEnd, pair0);
                            auto it1 = std::lower_bound(pairsBegin, pairsEnd, pair1);

                            // Track if forward pairs were found
                            bool found0 = (it0 != pairsEnd && *it0 == pair0);
                            bool found1 = (it1 != pairsEnd && *it1 == pair1);
                            
                            // Verify we found exact matches and link them
                            if (found0 && found1) {
                                uint64_t id0 = it0 - pairsBegin;
                                SHASTA_ASSERT(id0 < maxId);
                                uint64_t id1 = it1 - pairsBegin;
                                SHASTA_ASSERT(id1 < maxId);
                                
                                // Link these two position pairs in the disjoint sets
                                disjointSets.unite(id0, id1);
                                
                                // Store allele information
                                variantClusteringPositionPairAlleles[id0] = base0.value;
                                variantClusteringPositionPairAlleles[id1] = base1.value;

                                // Store marker context
                                auto& markerInfoContext0 = variantClusteringPositionPairContexts[id0];
                                markerInfoContext0.prevMarkerInfo = MarkerKmers::MarkerInfo(currentOrientedReadId0, segment.ordinalsA[0]);
                                markerInfoContext0.nextMarkerInfo = MarkerKmers::MarkerInfo(currentOrientedReadId0, segment.ordinalsB[0]);

                                auto& markerInfoContext1 = variantClusteringPositionPairContexts[id1];
                                markerInfoContext1.prevMarkerInfo = MarkerKmers::MarkerInfo(currentOrientedReadId1, segment.ordinalsA[1]);
                                markerInfoContext1.nextMarkerInfo = MarkerKmers::MarkerInfo(currentOrientedReadId1, segment.ordinalsB[1]);
                            }
                            
                            
                            // Create reverse complement position pairs
                            OrientedReadId currentOrientedReadId0Rc(currentOrientedReadId0.getReadId(), 1 - currentOrientedReadId0.getStrand());
                            OrientedReadId currentOrientedReadId1Rc(currentOrientedReadId1.getReadId(), 1 - currentOrientedReadId1.getStrand());
                            const uint64_t readLength0 = getReads().getReadRawSequenceLength(currentOrientedReadId0.getReadId());
                            const uint64_t readLength1 = getReads().getReadRawSequenceLength(currentOrientedReadId1.getReadId());
                            pair<OrientedReadId, uint32_t> pair0Rc(currentOrientedReadId0Rc, uint32_t(readLength0 - 1 - positionInRead0));
                            pair<OrientedReadId, uint32_t> pair1Rc(currentOrientedReadId1Rc, uint32_t(readLength1 - 1 - positionInRead1));

                            // Binary search to find IDs in sorted vector (O(log n) lookup)
                            auto it0Rc = std::lower_bound(variantClusteringPositionPairs.begin(), variantClusteringPositionPairs.end(), pair0Rc);
                            auto it1Rc = std::lower_bound(variantClusteringPositionPairs.begin(), variantClusteringPositionPairs.end(), pair1Rc);


                            // Track if reverse complement pairs were found
                            bool found0Rc = (it0Rc != pairsEnd && *it0Rc == pair0Rc);
                            bool found1Rc = (it1Rc != pairsEnd && *it1Rc == pair1Rc);
                            

                            // Verify we found exact matches and link them
                            if (found0Rc && found1Rc) {
                                uint64_t id0Rc = it0Rc - pairsBegin;
                                SHASTA_ASSERT(id0Rc < maxId);
                                uint64_t id1Rc = it1Rc - pairsBegin;
                                SHASTA_ASSERT(id1Rc < maxId);

                                // Link these two position pairs in the disjoint sets
                                disjointSets.unite(id0Rc, id1Rc);

                                // Find the base these position pairs represent in the reverse strand reads
                                Base base0Rc = getReads().getOrientedReadBase(currentOrientedReadId0Rc, readLength0 - 1 - positionInRead0);
                                Base base1Rc = getReads().getOrientedReadBase(currentOrientedReadId1Rc, readLength1 - 1 - positionInRead1);
                                variantClusteringPositionPairAlleles[id0Rc] = base0Rc.value;
                                variantClusteringPositionPairAlleles[id1Rc] = base1Rc.value;

                                // Store marker context for both positions (prev/next MarkerInfo) in Assembler member variantClusteringPositionPairContexts
                                auto& markerInfoContext0Rc = variantClusteringPositionPairContexts[id0Rc];
                                const uint32_t markerCount0Rc = uint32_t(markersRef[currentOrientedReadId0Rc.getValue()].size());
                                // Convert ordinals: ordinalsA becomes ordinalsB_rc, ordinalsB becomes ordinalsA_rc
                                markerInfoContext0Rc.prevMarkerInfo = MarkerKmers::MarkerInfo(currentOrientedReadId0Rc, markerCount0Rc - 1 - segment.ordinalsB[0]);
                                markerInfoContext0Rc.nextMarkerInfo = MarkerKmers::MarkerInfo(currentOrientedReadId0Rc, markerCount0Rc - 1 - segment.ordinalsA[0]);

                                auto& markerInfoContext1Rc = variantClusteringPositionPairContexts[id1Rc];
                                const uint32_t markerCount1Rc = uint32_t(markersRef[currentOrientedReadId1Rc.getValue()].size());
                                markerInfoContext1Rc.prevMarkerInfo = MarkerKmers::MarkerInfo(currentOrientedReadId1Rc, markerCount1Rc - 1 - segment.ordinalsB[1]);
                                markerInfoContext1Rc.nextMarkerInfo = MarkerKmers::MarkerInfo(currentOrientedReadId1Rc, markerCount1Rc - 1 - segment.ordinalsA[1]);
                            }

                        }
                        
                        // Increment position counters for both sequences
                        position0++;
                        position1++;
                        
                    } else if (hasBase0) {
                        // Only read 0 has a base (read 1 has a gap)
                        position0++;
                    } else if (hasBase1) {
                        // Only read 1 has a base (read 0 has a gap)
                        position1++;
                    }
                    // else: both have gaps (shouldn't happen in a valid alignment, but handle gracefully)
                }
            }
        }
    }
}





// Main function to create variant clusters
void Assembler::performGlobalVariantClustering(
    uint64_t minCoverage,
    uint64_t maxCoverage,
    size_t threadCount)
{
    
    const auto tTotalStart = steady_clock::now();
    
    performanceLog << timestamp << "Starting Global Variant Clustering" << endl;
    cout << timestamp << "Starting Global Variant Clustering" << endl;
    
    cout << "\n============================================" << endl;
    cout << "VARIANT CLUSTERING PROFILING" << endl;
    cout << "============================================\n" << endl;
    
    // Check prerequisites
    const auto tCheckStart = steady_clock::now();
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkAlignmentDataAreOpen();
    SHASTA_ASSERT(compressedAlignments.isOpen());
    const auto tCheckEnd = steady_clock::now();
    const double tCheck = seconds(tCheckEnd - tCheckStart);
    
    // Access position pairs that were collected during alignment computation
    const auto tAccessStart = steady_clock::now();
    performanceLog << timestamp << "Accessing position pairs from alignment computation" << endl;
    cout << "Accessing position pairs collected during alignment computation..." << endl;
    
    // The position pairs were already collected in variantClusteringPositionPairs during computeAlignments
    // We need to access them (they're already stored)
    if (!variantClusteringPositionPairs.isOpen) {
        cout << "ERROR: variantClusteringPositionPairs is not open!" << endl;
        return;
    }
    
    cout << "Found " << variantClusteringPositionPairs.size() << " position pairs from alignment computation" << endl;
    const auto tAccessEnd = steady_clock::now();
    const double tAccess = seconds(tAccessEnd - tAccessStart);

    



    // XXX
    // --- START OF: MINIMUM OCCURRENCES FILTER
    //     Sort, count occurrences, and filter position pairs
    //     This creates our "perfect hash" where the index is the OrientedReadId

    const auto tDeduplicateStart = steady_clock::now();
    performanceLog << timestamp << "Sorting, counting, and filtering position pairs" << endl;
    cout << "Sorting and counting position pair occurrences..." << endl;
    const auto tOccurrenceStart = steady_clock::now();
    
    // Sort using std::sort (works because MemoryMapped::Vector::begin/end return T*)
    std::sort(variantClusteringPositionPairs.begin(), variantClusteringPositionPairs.end());

    // Count occurrences of genomic positions (strand-agnostic) and filter
    const uint64_t minOccurrences = 3;
    const uint64_t totalOccurrences = variantClusteringPositionPairs.size();

    MemoryMapped::Vector<std::pair<OrientedReadId, uint32_t>> variantClusteringFilteredPositionPairs;
    variantClusteringFilteredPositionPairs.createNew(
        largeDataName("tmp-VariantClusteringFilteredPositionPairs"),
        largeDataPageSize);
    variantClusteringFilteredPositionPairs.reserve(variantClusteringPositionPairs.size());

    if (!variantClusteringPositionPairs.empty()) {
        uint64_t count = 1;
        auto current = variantClusteringPositionPairs[0];
    
        for (uint64_t i = 1; i < totalOccurrences; ++i) {
            if (variantClusteringPositionPairs[i] == current) {
                ++count;
            } else {
                if (count >= minOccurrences) {
                    variantClusteringFilteredPositionPairs.push_back(current);
                }
                current = variantClusteringPositionPairs[i];
                count = 1;
            }
        }
        if (count >= minOccurrences) {
            variantClusteringFilteredPositionPairs.push_back(current);
        }
    }


    // Sanity check: ensure each filtered pair still has its reverse-complement view
    {
        std::set<std::pair<OrientedReadId, uint32_t>> filteredSet(
            variantClusteringFilteredPositionPairs.begin(),
            variantClusteringFilteredPositionPairs.end());

        uint64_t missingRcCount = 0;
        const uint64_t maxPrint = 20;

        for (const auto& p : variantClusteringFilteredPositionPairs) {
            const ReadId readId = p.first.getReadId();
            const Strand strand = p.first.getStrand();
            const uint32_t position = p.second;
            const uint64_t readLength = getReads().getReadRawSequenceLength(readId);

            const Strand rcStrand = Strand(1 - strand);
            const uint32_t rcPosition = readLength - 1 - position;
            const OrientedReadId rcId(readId, rcStrand);
            const auto rcPair = std::make_pair(rcId, rcPosition);

            if (filteredSet.find(rcPair) == filteredSet.end()) {
                ++missingRcCount;
                if (missingRcCount <= maxPrint) {
                    cout << "Missing RC after filtering: "
                         << p.first << ":" << position
                         << " -> expected " << rcId << ":" << rcPosition << endl;
                }
            }
        }

        if (missingRcCount == 0) {
            cout << "✓ Filtered positions still include both strand views." << endl;
        } else {
            cout << "✗ WARNING: " << missingRcCount
                 << " filtered positions lost their reverse-complement view." << endl;
        }

        // count how many pairs are in 0th strand and how many are in 1st strand
        uint64_t count0 = 0;
        uint64_t count1 = 0;
        for (const auto& pair : variantClusteringFilteredPositionPairs) {
            if (pair.first.getStrand() == 0) {
                ++count0;
            } else {
                ++count1;
            }
        }
        cout << "Number of pairs in 0th strand: " << count0 << endl;
        cout << "Number of pairs in 1st strand: " << count1 << endl;
    }

    const auto tOccurrenceEnd = steady_clock::now();
    const double tOccurrence = seconds(tOccurrenceEnd - tOccurrenceStart);

    cout << "After filtering (min occurrences=" << minOccurrences << "): " 
         << variantClusteringFilteredPositionPairs.size() << " position pairs (from " << totalOccurrences << " total occurrences)" << endl;
    cout << "  Filtered out " << totalOccurrences - variantClusteringFilteredPositionPairs.size() << " position pairs with < " << minOccurrences << " occurrences" << endl;
    

    // XXX
    // --- END OF: MINIMUM OCCURRENCES FILTER
    // 

    // Replace the memory-mapped vector with filtered results
    variantClusteringPositionPairs.clear();
    variantClusteringPositionPairs.reserve(variantClusteringFilteredPositionPairs.size());
    for (const auto& pair : variantClusteringFilteredPositionPairs) {
        variantClusteringPositionPairs.push_back(pair);
    }
    variantClusteringPositionPairs.unreserve();







    // // XXX
    // // --- START OF: WELL-SEPARATED FILTER
    // //     Filter out clusters of nearby SNPs (well-separated filter)
    // //     Sequencing error and artifacts often appear as clusters of nearby SNPs.
    // //     To avoid clusters of errors, the informative SNPs need to be well-separated.
    // //     Only SNPs at least 32bp apart are considered.
    // //     Since variantClusteringFilteredPositionPairs is already sorted by (OrientedReadId, position),
    // //     positions from the same read are grouped together - we can do a single pass!
    
    // const auto tSeparationStart = steady_clock::now();

    // // Remember how many survived the occurrence filter.
    // const uint64_t filteredCountBeforeSeparation = variantClusteringFilteredPositionPairs.size();

    // // --- Filter out clusters of nearby SNPs (well-separated filter) ---
    // const uint64_t minSeparation = 32;

    // MemoryMapped::Vector<pair<OrientedReadId, uint32_t>> wellSeparatedPositionPairs;
    // wellSeparatedPositionPairs.createNew(
    //     largeDataName("tmp-VariantClusteringWellSeparatedPositionPairs"),
    //     largeDataPageSize);
    // wellSeparatedPositionPairs.reserve(variantClusteringFilteredPositionPairs.size());


    // if (!variantClusteringFilteredPositionPairs.empty()) {
    //     // Track last kept canonical (strand 0) position per read.
    //     std::unordered_map<ReadId, uint32_t> lastKeptStrand0Pos;

    //     for (const auto& pair : variantClusteringFilteredPositionPairs) {
    //         const OrientedReadId orientedReadId = pair.first;
    //         const ReadId readId = orientedReadId.getReadId();
    //         const Strand strand = orientedReadId.getStrand();

    //         if (strand != 0) {
    //             continue;   // decisions are made on strand-0 entries only
    //         }

    //         const uint32_t position = pair.second;
    //         const uint64_t readLength = getReads().getReadRawSequenceLength(readId);

    //         auto [it, inserted] = lastKeptStrand0Pos.try_emplace(readId, position);
    //         bool keep = false;

    //         if (inserted) {
    //             keep = true;           // first canonical position for this read
    //         } else {
    //             const uint32_t lastPos = it->second;
    //             const uint32_t distance =
    //                 position > lastPos ? position - lastPos : lastPos - position;
    //             if (distance == 0 || distance >= minSeparation) {
    //                 keep = true;
    //             }
    //         }

    //         if (keep) {
    //             it->second = position;
    //             wellSeparatedPositionPairs.push_back(pair);

    //             // Also push the reverse-complement view for this genomic position.
    //             const OrientedReadId rcId(readId, Strand(1));
    //             const uint32_t rcPosition = uint32_t(readLength - 1 - position);
    //             wellSeparatedPositionPairs.push_back(std::make_pair(rcId, rcPosition));
    //         }
    //     }
    // }


    // const uint64_t wellSeparatedCount = wellSeparatedPositionPairs.size();
    // const uint64_t wellSeparatedFilteredOut = filteredCountBeforeSeparation - wellSeparatedCount;

    // std::cout << "After well-separated filter (min separation="
    //             << minSeparation << "bp): " << wellSeparatedCount
    //             << " position pairs" << std::endl;
    // std::cout << "  Filtered out " << wellSeparatedFilteredOut
    //             << " position pairs within " << minSeparation
    //             << "bp of adjacent positions" << std::endl;

    // const auto tSeparationEnd = steady_clock::now();
    // const double tSeparation = seconds(tSeparationEnd - tSeparationStart);


    // // Sanity check: each canonical position should have both strand views
    // {
    //     struct PairHash {
    //         size_t operator()(const std::pair<ReadId, uint32_t>& k) const noexcept {
    //             return std::hash<ReadId>()(k.first) ^ (std::hash<uint32_t>()(k.second) << 1);
    //         }
    //     };

    //     std::unordered_map<std::pair<ReadId, uint32_t>, uint8_t, PairHash> viewMask;
    //     viewMask.reserve(wellSeparatedPositionPairs.size());

    //     for (const auto& p : wellSeparatedPositionPairs) {
    //         const ReadId readId = p.first.getReadId();
    //         const Strand strand = p.first.getStrand();
    //         const uint32_t position = p.second;
    //         const uint64_t readLength = getReads().getReadRawSequenceLength(readId);
    //         const uint32_t strand0Pos = (strand == 0) ? position : uint32_t(readLength - 1 - position);

    //         auto& mask = viewMask[{readId, strand0Pos}];
    //         mask |= (1u << strand);   // bit 0 for strand 0, bit 1 for strand 1
    //     }

    //     uint64_t missingViews = 0;
    //     for (const auto& [key, mask] : viewMask) {
    //         if (mask != 0b11) {
    //             ++missingViews;
    //             if (missingViews <= 20) {
    //                 const ReadId readId = key.first;
    //                 const uint32_t strand0Pos = key.second;
    //                 const uint64_t readLength = getReads().getReadRawSequenceLength(readId);
    //                 cout << "Missing reverse-complement view for read " << readId
    //                     << " canonical position " << strand0Pos
    //                     << " → expected strand 0:" << strand0Pos
    //                     << " and strand 1:" << (readLength - 1 - strand0Pos) << endl;
    //             }
    //         }
    //     }

    //     if (missingViews == 0) {
    //         cout << "✓ Sanity check: every canonical position has both strand views. Checked " << viewMask.size() << " canonical positions." << endl;
    //     } else {
    //         cout << "✗ WARNING: " << missingViews
    //             << " canonical positions lost a strand view after filtering. Checked " << viewMask.size() << " canonical positions." << endl;
    //     }

    //     // count how many pairs are in 0th strand and how many are in 1st strand
    //     uint64_t count0 = 0;
    //     uint64_t count1 = 0;
    //     for (const auto& pair : wellSeparatedPositionPairs) {
    //         if (pair.first.getStrand() == 0) {
    //             ++count0;
    //         } else {
    //             ++count1;
    //         }
    //     }
    //     cout << "Number of pairs in 0th strand: " << count0 << endl;
    //     cout << "Number of pairs in 1st strand: " << count1 << endl;

    //     // // Print all pairs that contain readId 0
    //     // for (const auto& pair : wellSeparatedPositionPairs) {
    //     //     if (pair.first.getReadId() == 0) {
    //     //         cout << pair.first << ":" << pair.second << endl;
    //     //     }
    //     // }
    // }


    // const auto tDeduplicateEnd = steady_clock::now();
    // const double tDeduplicate = seconds(tDeduplicateEnd - tDeduplicateStart);
    // cout << "  Occurrence filter time: " << tOccurrence << " s" << endl;
    // cout << "  Well-separated filter time: " << tSeparation << " s" << endl;
    // performanceLog << timestamp << "Occurrence filter time " << tOccurrence << " s" << endl;
    // performanceLog << timestamp << "Well-separated filter time " << tSeparation << " s" << endl;


    // // Replace the memory-mapped vector with filtered results
    // variantClusteringPositionPairs.clear();
    // variantClusteringPositionPairs.reserve(wellSeparatedPositionPairs.size());
    // for (const auto& pair : wellSeparatedPositionPairs) {
    //     variantClusteringPositionPairs.push_back(pair);
    // }
    // variantClusteringPositionPairs.unreserve();

    
    // // XXX
    // // --- END OF: WELL-SEPARATED FILTER
    // // 
















    if (variantClusteringPositionPairs.empty()) {
        cout << "No position pairs found. Exiting." << endl;
        return;
    }
    
    // Initialize disjoint sets for linking position pairs
    const auto tDisjointSetInitStart = steady_clock::now();
    performanceLog << timestamp << "Initializing disjoint sets" << endl;
    cout << "Initializing disjoint sets with " << variantClusteringPositionPairs.size() << " elements..." << endl;
    
    const uint64_t disjointSetCount = variantClusteringPositionPairs.size();
    
    // Store clustering data in Assembler members so it persists for mode3Assembly
    variantClusteringDisjointSetTable.createNew(
        largeDataName("tmp-VariantClusteringDisjointSets"), 
        largeDataPageSize);
    variantClusteringDisjointSetTable.resize(disjointSetCount);
    variantClusteringDisjointSets = std::make_shared<DisjointSets>(variantClusteringDisjointSetTable.begin(), disjointSetCount);
    
    // Initialize allele storage for position pairs
    variantClusteringPositionPairAlleles.createNew(
        largeDataName("tmp-VariantClusteringPositionPairAlleles"),
        largeDataPageSize);
    variantClusteringPositionPairAlleles.resize(disjointSetCount);
    
    // Initialize position context storage for position pairs
    variantClusteringPositionPairContexts.createNew(
        largeDataName("tmp-VariantClusteringPositionPairContexts"),
        largeDataPageSize);
    variantClusteringPositionPairContexts.resize(disjointSetCount);
    
    // // Update local pointer to use Assembler member
    // variantClusteringData.disjointSetsPointer = variantClusteringDisjointSets;
    
    cout << "Disjoint sets, allele and context storage initialized" << endl;
    const auto tDisjointSetInitEnd = steady_clock::now();
    const double tDisjointSetInit = seconds(tDisjointSetInitEnd - tDisjointSetInitStart);

    // Phase 2: Re-process alignments to link position pairs using disjoint sets
    const auto tPass2Start = steady_clock::now();
    performanceLog << timestamp << "Phase 2: Linking position pairs with disjoint sets" << endl;
    cout << "\nPhase 2: Linking position pairs with disjoint sets..." << endl;
    
    // Pick the batch size for load balancing alignments
    const uint64_t requestedBatchSize = 1;  // Process alignments in batches
    setupLoadBalancing(alignmentData.size(), requestedBatchSize);
    runThreads(&Assembler::linkVariantClustersThreadFunction, threadCount);

    const auto tPass2End = steady_clock::now();
    const double tPass2 = seconds(tPass2End - tPass2Start);
    cout << "Phase 2 complete" << endl;
    



    {
        // Verify reverse complement consistency
        const auto tVerifyStart = steady_clock::now();
        performanceLog << timestamp << "Verifying reverse complement consistency" << endl;
        cout << "\nVerifying reverse complement consistency..." << endl;

        // Build reverse complement index map: (readId, strand, position) -> index
        std::map<pair<ReadId, pair<Strand, uint32_t>>, uint64_t> pairToIndex;
        for (uint64_t i = 0; i < variantClusteringPositionPairs.size(); i++) {
            const auto& p = variantClusteringPositionPairs[i];
            pairToIndex[{p.first.getReadId(), {p.first.getStrand(), p.second}}] = i;
        }

        uint64_t inconsistencies = 0;
        uint64_t checkedPairs = 0;
        const uint64_t sampleSize = std::min(variantClusteringPositionPairs.size(), uint64_t(10000));

        // Sample check: verify that if (readId, 0, pos) is linked with others,
        // then (readId, 1, readLength-1-pos) is linked with corresponding RC pairs
        for (uint64_t i = 0; i < sampleSize; i++) {
            const auto& pair = variantClusteringPositionPairs[i];
            const ReadId readId = pair.first.getReadId();
            const Strand strand = pair.first.getStrand();
            const uint32_t position = pair.second;
            
            // Only check strand 0 pairs
            if (strand != 0) continue;
            
            checkedPairs++;
            
            // Find reverse complement pair
            const uint64_t readLength = getReads().getReadRawSequenceLength(readId);
            const uint32_t positionRc = readLength - 1 - position;
            const auto rcKey = make_pair(readId, make_pair(Strand(1), positionRc));
            
            auto rcIt = pairToIndex.find(rcKey);
            if (rcIt == pairToIndex.end()) {
                continue; // RC pair was filtered out, skip
            }
            
            const uint64_t rcIndex = rcIt->second;
            const uint64_t setId0 = variantClusteringDisjointSets->find(i);
            const uint64_t setId0Rc = variantClusteringDisjointSets->find(rcIndex);
            
            // Find all pairs in the same set as i
            std::vector<uint64_t> linkedPairs;
            for (uint64_t j = 0; j < variantClusteringPositionPairs.size(); j++) {
                if (variantClusteringDisjointSets->find(j) == setId0) {
                    linkedPairs.push_back(j);
                }
            }
            
            // For each linked pair j (strand 0), check if rc(j) is linked with rc(i)
            for (uint64_t j : linkedPairs) {
                if (i == j) continue;
                
                const auto& pairJ = variantClusteringPositionPairs[j];
                if (pairJ.first.getStrand() != 0) continue; // Only check strand 0 pairs
                
                const ReadId readIdJ = pairJ.first.getReadId();
                const uint64_t readLengthJ = getReads().getReadRawSequenceLength(readIdJ);
                const uint32_t positionJRc = readLengthJ - 1 - pairJ.second;
                const auto rcJKey = make_pair(readIdJ, make_pair(Strand(1), positionJRc));
                
                auto rcJIt = pairToIndex.find(rcJKey);
                if (rcJIt == pairToIndex.end()) continue;
                
                const uint64_t rcJIndex = rcJIt->second;
                const uint64_t setIdJRc = variantClusteringDisjointSets->find(rcJIndex);
                
                // Check consistency: rc(i) and rc(j) should be in the same set
                if (setIdJRc != setId0Rc) {
                    inconsistencies++;
                    if (inconsistencies <= 5) {
                        cout << "INCONSISTENCY: " << pair.first << ":" << position 
                            << " linked with " << pairJ.first << ":" << pairJ.second
                            << ", but RC pairs are NOT linked" << endl;
                    }
                }
            }
        }

        const auto tVerifyEnd = steady_clock::now();
        const double tVerify = seconds(tVerifyEnd - tVerifyStart);

        cout << "Verification: checked " << checkedPairs << " pairs, found " 
            << inconsistencies << " inconsistencies" << endl;
        if (inconsistencies == 0) {
            cout << "✓ Reverse complement consistency verified!" << endl;
        } else {
            cout << "✗ WARNING: Found inconsistencies - check Phase 2 linking logic!" << endl;
        }
        performanceLog << timestamp << "Verification: " << checkedPairs << " pairs checked, " 
                    << inconsistencies << " inconsistencies" << endl;
    }
    




    
    // Identify cluster representatives
    const auto tIdentifyClustersStart = steady_clock::now();
    performanceLog << timestamp << "Identifying cluster representatives" << endl;
    cout << "Identifying cluster representatives..." << endl;
    
    // Store cluster representatives in Assembler member so they persist for mode3Assembly
    variantClusteringClusterRepresentatives.clear();
    for (uint64_t id = 0; id < disjointSetCount; id++) {
        if (variantClusteringDisjointSets->find(id) == id) {
            variantClusteringClusterRepresentatives.push_back(id);
        }
    }
    
    cout << "Found " << variantClusteringClusterRepresentatives.size() << " clusters from " 
         << disjointSetCount << " unique position pairs" << endl;

    const auto tIdentifyClustersEnd = steady_clock::now();
    const double tIdentifyClusters = seconds(tIdentifyClustersEnd - tIdentifyClustersStart);
    
    // Print timing summary
    const auto tTotalEnd = steady_clock::now();
    const double tTotal = seconds(tTotalEnd - tTotalStart);
    
    cout << "\n============================================" << endl;
    cout << "VARIANT CLUSTERING TIMING SUMMARY" << endl;
    cout << "============================================" << endl;
    cout << std::left << std::setw(40) << "Step" << std::right << std::setw(12) << "Time (s)" << std::setw(12) << "Percent" << endl;
    cout << std::string(64, '-') << endl;
    cout << std::left << std::setw(40) << "Prerequisites check" << std::right << std::setw(12) << tCheck << std::setw(12) << (100.0 * tCheck / tTotal) << endl;
    cout << std::left << std::setw(40) << "Access position pairs" << std::right << std::setw(12) << tAccess << std::setw(12) << (100.0 * tAccess / tTotal) << endl;
    cout << std::left << std::setw(40) << "Occurrence filter" << std::right << std::setw(12) << tOccurrence << std::setw(12) << (100.0 * tOccurrence / tTotal) << endl;
    // cout << std::left << std::setw(40) << "Well-separated filter" << std::right << std::setw(12) << tSeparation << std::setw(12) << (100.0 * tSeparation / tTotal) << endl;
    cout << std::left << std::setw(40) << "Initialize disjoint sets" << std::right << std::setw(12) << tDisjointSetInit << std::setw(12) << (100.0 * tDisjointSetInit / tTotal) << endl;
    cout << std::left << std::setw(40) << "Phase 2: Link pairs" << std::right << std::setw(12) << tPass2 << std::setw(12) << (100.0 * tPass2 / tTotal) << endl;
    cout << std::left << std::setw(40) << "Identify clusters" << std::right << std::setw(12) << tIdentifyClusters << std::setw(12) << (100.0 * tIdentifyClusters / tTotal) << endl;
    cout << std::string(64, '-') << endl;
    cout << std::left << std::setw(40) << "Total Time" << std::right << std::setw(12) << tTotal << std::setw(12) << "100.0" << endl;
    cout << "============================================\n" << endl;
    
    performanceLog << timestamp << "createVariantClusters ends" << endl;
}