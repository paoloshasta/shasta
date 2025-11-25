#include "Assembler.hpp"
#include "timestamp.hpp"
#include "deduplicate.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <mutex>
#include "Reads.hpp"

using namespace shasta;

void Assembler::createReadGraph5()
{
    cout << timestamp << "createReadGraph5 begins." << endl;

    // Check that we have the necessary data.
    checkMarkersAreOpen();
    checkAlignmentDataAreOpen();
    checkVariantClusteringPositionPairsIsOpen();

    if(!variantClusteringDisjointSets) {
        throw runtime_error("Variant clustering disjoint sets not available.");
    }

    // 1. Build membersByRepIdx (Cluster -> PositionPairs map)
    // We use the 2-pass approach for MemoryMapped::VectorOfVectors to be memory efficient
    // and avoid creating a huge temporary vector<vector> in RAM.
    cout << timestamp << "Building cluster to reads map." << endl;
    
    const uint64_t positionPairCount = variantClusteringPositionPairs.size();
    const uint64_t clusterCount = variantClusteringDisjointSets->size();

    variantClusteringMembersByRepIdx.createNew(
        largeDataName("VariantClusteringMembersByRepIdx"),
        largeDataPageSize
    );

    // Pass 1: Count the number of members in each cluster (representative).
    // This prepares the "Table of Contents" (TOC).
    variantClusteringMembersByRepIdx.beginPass1(clusterCount);
    for(uint64_t i=0; i<positionPairCount; i++) {
        const uint64_t clusterId = variantClusteringDisjointSets->find(i);
        variantClusteringMembersByRepIdx.incrementCount(clusterId);
    }

    // Pass 2: Store the member IDs.
    // This fills the data vector.
    variantClusteringMembersByRepIdx.beginPass2();
    for(uint64_t i=0; i<positionPairCount; i++) {
        const uint64_t clusterId = variantClusteringDisjointSets->find(i);
        variantClusteringMembersByRepIdx.store(clusterId, i);
    }
    variantClusteringMembersByRepIdx.endPass2();
    
    cout << timestamp << "Map built. " << clusterCount << " potential clusters." << endl;

    // 2. Initialize output data structures.
    // It will be used as a lookup table (or boolean array) to mark valid clusters.
    variantClusteringValidClusters.createNew(
        largeDataName("VariantClusteringValidClusters"),
        largeDataPageSize);
    variantClusteringValidClusters.resize(clusterCount);
    // Initialize to 0 (invalid)
    std::fill(variantClusteringValidClusters.begin(), variantClusteringValidClusters.end(), 0);

    // 3. Run threads to perform global clustervalidity check.
    // A cluster is valid if it has at least 2 alleles with coverage >= minAlleleCoverage.
    cout << timestamp << "Running global cluster validity checks." << endl;
    const uint32_t threadCount = std::thread::hardware_concurrency();
    const uint64_t minAlleleCoverage = 5;
    setupLoadBalancing(clusterCount, 1);
    runThreads(&Assembler::computeClusterValidityThreadFunction, threadCount);
    cout << timestamp << "Global cluster validity checks completed." << endl;


    // 4. Prepare for Haplotype Voting
    variantClusteringValidClustersCompatible.createNew(
        largeDataName("VariantClusteringValidClustersCompatible"),
        largeDataPageSize
    );
    variantClusteringValidClustersCompatible.resize(clusterCount);
    std::fill(variantClusteringValidClustersCompatible.begin(), variantClusteringValidClustersCompatible.end(), 0);

    // 4. Initialize Global Haplotype Graph
    // Calculate total oriented reads (2 * readCount)
    const uint64_t orientedReadCount = 2 * reads->readCount();
    
    // Create the graph with orientedReadCount vertices
    cout << timestamp << "Initializing Global Haplotype Graph with " << orientedReadCount << " vertices." << endl;
    globalHaplotypeGraph = std::make_shared<HaplotypeGraph>(orientedReadCount);

    // 5. Run threads to generate votes and add edges directly
    cout << timestamp << "Running compatibility checks and building haplotype graph." << endl;
    setupLoadBalancing(getReads().readCount(), 1); 
    runThreads(&Assembler::createReadGraph5ThreadFunction, threadCount);
    cout << timestamp << "Haplotype graph construction completed." << endl;
    
    cout << timestamp << "Global Haplotype Graph built with " 
         << boost::num_vertices(*globalHaplotypeGraph) << " vertices and "
         << boost::num_edges(*globalHaplotypeGraph) << " edges." << endl;

    // 7. Refine Clusters using the Graph
    cout << timestamp << "Refining clusters using haplotype graph." << endl;
    
    // This vector is used to track the quality or "membership status" of each read within a variant cluster.
    // Here is the breakdown:
    // 1. Initialization: It creates a new memory-mapped vector sized to positionPairCount (the total number of read-variant pairs).
    // 2. Default State: It fills the vector with 0, which represents a "Good" or "Keep" status.
    // 3. Usage: Later in the refineClustersThreadFunction this vector is updated. If a read is identified 
    // as "stray" (meaning it connects more strongly to reads of a different allele in the haplotype graph than to its own), 
    // then its status is changed to 1 (Filter/Stray).
    // Essentially, this prepares the "report card" for every read in every cluster, assuming they are all valid until 
    // proven otherwise by the graph refinement step.
    variantClusteringMemberStatus.createNew(
        largeDataName("VariantClusteringMemberStatus"),
        largeDataPageSize
    );
    variantClusteringMemberStatus.resize(positionPairCount);
    std::fill(variantClusteringMemberStatus.begin(), variantClusteringMemberStatus.end(), 0); // Default 0 = Good

    setupLoadBalancing(clusterCount, 1);
    runThreads(&Assembler::refineClustersThreadFunction, threadCount);
    
    cout << timestamp << "Cluster refinement completed." << endl;

    // Print how many clusters are valid and how many are compatible.
    uint64_t validClusterCount = 0;
    uint64_t compatibleClusterCount = 0;
    for(uint64_t i=0; i<clusterCount; i++) {
        if(variantClusteringValidClusters[i]) validClusterCount++;
        if(variantClusteringValidClustersCompatible[i]) compatibleClusterCount++;
    }
    cout << timestamp << "Valid clusters: " << validClusterCount << endl;
    cout << timestamp << "Compatible clusters: " << compatibleClusterCount << endl;


    //
    // TODO: REMOVE THIS AFTER TESTING.
    // KEEP ALL ALIGNMENTS FOR NOW.
    //
    cout << timestamp << "createReadGraph5 begins (including all alignments)." << endl;

    // Get the total number of stored alignments.
    const uint64_t alignmentCount = alignmentData.size();
    SHASTA_ASSERT(compressedAlignments.size() == alignmentCount);

    // Flag all alignments to be kept (no filtering).
    vector<bool> keepAlignment(alignmentCount, true);

    // Also mark all alignments as included in the read graph.
    for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
        alignmentData[alignmentId].info.isInReadGraph = 1;
    }

    cout << "Keeping all " << alignmentCount << " alignments." << endl;

    // Create the read graph using all alignments.
    createReadGraphUsingSelectedAlignments(keepAlignment);
}


// Global Cluster Validity Check Function to check if a cluster is valid.
// A cluster is valid if it has at least 2 alleles with coverage >= minAlleleCoverage.
void Assembler::computeClusterValidityThreadFunction(uint64_t threadId) {
    uint64_t currentMinAlleleCoverage = this->minAlleleCoverage;
    // Loop over a batch of cluster IDs assigned to this thread
    uint64_t clusterIdBegin, clusterIdEnd;
    while (getNextBatch(clusterIdBegin, clusterIdEnd)) {
        for (uint64_t clusterId = clusterIdBegin; clusterId < clusterIdEnd; clusterId++) {
            // variantClusteringMembersByRepIdx (Cluster -> PositionPairs map)
            const auto& members = variantClusteringMembersByRepIdx[clusterId];
            
            // Optimization: Simple array counters if we assume 1 pos per read per cluster
            // (Or use a thread-local scratch buffer to avoid allocation if uniqueness check is needed)
            std::array<uint32_t, 5> alleleCounts = {0};

            // Iterate members to count alleles
            for (uint64_t memberIdx : members) {
                // variantClusteringPositionPairAlleles (PositionPair -> Allele map)
                const uint8_t allele = variantClusteringPositionPairAlleles[memberIdx];
                alleleCounts[allele]++;
            }

            // Check threshold
            int validAlleles = 0;
            for(int a=0; a<5; a++) 
                if(alleleCounts[a] >= currentMinAlleleCoverage) validAlleles++;
                
            if (validAlleles >= 2) {
                variantClusteringValidClusters[clusterId] = 1;
            }
        }
    }
}


void Assembler::createReadGraph5ThreadFunction(uint64_t threadId)
{
    // Access data
    const auto& positionPairs = variantClusteringPositionPairs;
    auto& disjointSets = *variantClusteringDisjointSets;
    const auto& membersByRepIdx = variantClusteringMembersByRepIdx;

    uint64_t readIdBegin, readIdEnd;
    while (getNextBatch(readIdBegin, readIdEnd)) {
        for (ReadId readId = readIdBegin; readId < readIdEnd; readId++) {
            
            // Only process strand 0
            // The process is symmetric for the two strands.
            OrientedReadId currentOrientedReadId0(readId, 0);
            
            // Find clusters involving this read on strand 0.
            // The position pairs are sorted by read id and strand.
            auto it = std::lower_bound(positionPairs.begin(), positionPairs.end(), make_pair(currentOrientedReadId0, 0u));
            
            struct ClusterOnRead {
                uint64_t clusterId;
                uint32_t position;
                uint64_t indexInPositionPairs;
            };
            vector<ClusterOnRead> clusters;

            // Find clusters that involve this read on strand 0.
            while(it != positionPairs.end() && it->first == currentOrientedReadId0) {
                uint64_t index = it - positionPairs.begin();
                uint64_t clusterId = disjointSets.find(index);
                clusters.push_back({clusterId, it->second, index});
                ++it;
            }

            // If no clusters that involve this read on strand 0, skip it
            if (clusters.empty()) {
                continue;
            }

            vector<ClusterOnRead> filteredClusters;
            for(const auto& cluster : clusters) {
                // O(1) Lookup
                if (variantClusteringValidClusters[cluster.clusterId]) {
                    filteredClusters.push_back(cluster);
                }
            }
            
            // If no valid clusters after filtering, skip this read
            if (filteredClusters.empty()) {
                continue;
            }
            
            // Use filtered clusters for the rest of the analysis
            clusters = std::move(filteredClusters);




            // --------------------------------------------------------------
            // Compatibility Check with Multi-Allelic Splitting
            // --------------------------------------------------------------

            // Check if we have any clusters to process
            if (clusters.empty()) {
                continue;
            }

            // --- Step 1: Generate Virtual Clusters (Nodes for DP) ---
            // A physical cluster with alleles {Target, A, B} becomes 
            // two virtual clusters: (Target vs A) and (Target vs B).
            
            struct VirtualCluster {
                uint64_t clusterId;
                uint8_t targetAllele;
                uint8_t altAllele;
                // ID in the original 'clusters' vector to retrieve members
                uint32_t originalIndex; 
            };
            vector<VirtualCluster> dpNodes;
            dpNodes.reserve(clusters.size() * 2); // Heuristic reserve

            const uint64_t currentMinAlleleCoverage = this->minAlleleCoverage; 

            for(uint32_t i=0; i<clusters.size(); i++) {
                const auto& cluster = clusters[i];
                
                // Identify Target Allele
                uint8_t targetAllele = 255;
                if (cluster.indexInPositionPairs < variantClusteringPositionPairAlleles.size()) {
                    targetAllele = variantClusteringPositionPairAlleles[cluster.indexInPositionPairs];
                }
                SHASTA_ASSERT(targetAllele < 5); // Invalid target allele

                // Count alleles in this cluster to find significant ALTs
                const auto& members = membersByRepIdx[cluster.clusterId];
                std::array<uint32_t, 5> alleleCounts = {0};
                
                for(uint64_t memberIdx : members) {
                    if(memberIdx < variantClusteringPositionPairAlleles.size()) {
                        uint8_t a = variantClusteringPositionPairAlleles[memberIdx];
                        if(a < 5) alleleCounts[a]++;
                    }
                }

                // Create a Virtual Cluster for each significant ALT allele
                // It is guaranteed that there is at least one significant ALT allele.
                // The point is to split multiallelic clusters into binary het clusters.
                for(uint8_t a=0; a<5; a++) {
                    if (a == targetAllele) continue;
                    if (alleleCounts[a] >= currentMinAlleleCoverage) {
                        dpNodes.push_back({cluster.clusterId, targetAllele, a, i});
                    }
                }
            }

            const size_t N = dpNodes.size();
            if (N == 0) continue;

            // --- Step 2: Build Read Phase Vectors for each DP Node ---
            // nodeReads[i] contains list of {ReadId, isPhase1} for dpNodes[i]
            // Phase 0 = Matches Target. Phase 1 = Matches Alt. 

            struct ReadPhase {
                OrientedReadId::Int orientedReadIdValue; // Using the raw integer value
                bool isPhase1; // false=Target, true=Alt
            };
            vector<vector<ReadPhase>> nodeReads(N);

            for(size_t i=0; i<N; i++) {
                const auto& node = dpNodes[i];
                const auto& members = membersByRepIdx[node.clusterId];
                nodeReads[i].reserve(members.size());

                for(uint64_t memberIdx : members) {
                    const auto& pp = positionPairs[memberIdx];
                    
                    OrientedReadId orientedReadId = pp.first;
                    
                    // Skip target oriented read itself (given that the clusters involve this read on strand 0).
                    if(orientedReadId == currentOrientedReadId0) continue;

                    if(memberIdx < variantClusteringPositionPairAlleles.size()) {
                        uint8_t a = variantClusteringPositionPairAlleles[memberIdx];
                        
                        if (a == node.targetAllele) {
                            nodeReads[i].push_back({orientedReadId.getValue(), false}); 
                        } else if (a == node.altAllele) {
                            nodeReads[i].push_back({orientedReadId.getValue(), true});  
                        }
                    }
                }
                
                // Sort by OrientedReadId value for fast intersection
                std::sort(nodeReads[i].begin(), nodeReads[i].end(), 
                    [](const ReadPhase& a, const ReadPhase& b) { 
                        return a.orientedReadIdValue < b.orientedReadIdValue; 
                    });
            }

            // --- Step 3: DP for Longest Compatible Group (LCG) ---
            vector<int> LCG(N, 1);
            vector<int> parent(N, -1);

            for(int i = 0; i < N; i++) {
                for(int j = 0; j < i; j++) {
                    // Don't link two virtual nodes coming from the SAME physical cluster.
                    // This is to avoid linking multiallelic sites coming from the same physical cluster.
                    // We force the DP chain to pick at most one virtual node per physical genomic site, 
                    // ensuring we don't double-count the same site or create invalid paths.
                    // This enforces the constraint: "One genomic location = One node in the path."
                    if (dpNodes[i].clusterId == dpNodes[j].clusterId) continue;

                    // Check compatibility(i, j)
                    // Compatible if NO read conflicts.
                    // Conflict: Read r present in both, but has Phase 0 in one and Phase 1 in other.
                    
                    bool compatible = true;
                    const auto& orientedReadsI = nodeReads[i];
                    const auto& orientedReadsJ = nodeReads[j];
                    
                    auto itI = orientedReadsI.begin();
                    auto itJ = orientedReadsJ.begin();
                    
                    int supportPhase0 = 0; // Reads linking Target -> Target
                    int supportPhase1 = 0; // Reads linking Alt -> Alt

                    // Linear intersection scan to compare sorted sets.
                    while(itI != orientedReadsI.end() && itJ != orientedReadsJ.end()) {
                        if(itI->orientedReadIdValue < itJ->orientedReadIdValue) {
                            ++itI;
                        } else if(itJ->orientedReadIdValue < itI->orientedReadIdValue) {
                            ++itJ;
                        } else {
                            // Same OrientedRead covers both.
                            if(itI->isPhase1 != itJ->isPhase1) {
                                // CONFLICT found! (Phase 0 -> 1 or 1 -> 0).
                                compatible = false; 
                                break; 
                            } else {
                                // Consistent!
                                if (itI->isPhase1 == false) supportPhase0++;
                                else supportPhase1++;
                            }
                            ++itI;
                            ++itJ;
                        }
                    }

                    // COMPATIBILITY CRITERIA:
                    // 1. Must be compatible (no conflicts)
                    // 2. Must have ACTIVE support (shared reads)
                    
                    if (compatible) {
                        // "Strong Linkage" (Support for BOTH alleles)
                        // This ensures we are tracking two real haplotypes, not just linking noise.
                        if (supportPhase0 > 0 && supportPhase1 > 0) {
                            if(LCG[j] + 1 > LCG[i]) {
                                LCG[i] = LCG[j] + 1; // Longest compatible group length
                                parent[i] = j; // Parent node in the DP chain
                            }
                        }
                    }
                }
            }

            // --- Step 4: Traceback and Chain Identification (Multi-Chain) ---
            // This is the core of the DP algorithm.
            // It identifies the longest compatible group of virtual nodes.
            // It also identifies isolated sites that are supported by a sufficiently high number of reads.
            // It then marks the clusters as valid and collects the reads that support the Target and Alt alleles.
 
            vector<bool> isAssigned(N, false);
            vector<pair<int, int>> sortedLCGIndices; // (Length, Index)
            vector<pair<int, int>> isolatedSitesIndices;

            for(int i = 0; i < N; i++) {
                if (LCG[i] > 1) {
                    sortedLCGIndices.push_back({LCG[i], i});
                } else {
                    isolatedSitesIndices.push_back({LCG[i], i});
                }
            }

            // Sort chains by length descending
            std::sort(sortedLCGIndices.rbegin(), sortedLCGIndices.rend());

            vector<vector<int>> validChains;
            
            // 1. Process Multi-Node Chains
            for (const auto& p : sortedLCGIndices) {
                int endNode = p.second;
                if (!isAssigned[endNode]) {
                    vector<int> chain;
                    int curr = endNode;
                    while (curr != -1 && !isAssigned[curr]) {
                        chain.push_back(curr);
                        isAssigned[curr] = true;
                        curr = parent[curr];
                    }
                    // No need to reverse if we just iterate, but for correctness:
                    std::reverse(chain.begin(), chain.end());
                    validChains.push_back(chain);
                }
            }

            // 2. Process Isolated Sites
            // An isolated site is considered informative only if the reference read allele 
            // is supported by a sufficiently high number of reads.
            const uint64_t minIsolatedCoverage = 7; // Stricter threshold for isolated sites
            
            for (const auto& p : isolatedSitesIndices) {
                int nodeIdx = p.second;
                if (!isAssigned[nodeIdx]) {
                    // Check support for this specific Virtual Cluster
                    
                    // Count reads supporting the reference read allele (Target Allele / Phase 0)
                    uint64_t targetAlleleSupport = 0;
                    for (const auto& rp : nodeReads[nodeIdx]) {
                        if (!rp.isPhase1) { // isPhase1=false means Target Allele
                            targetAlleleSupport++;
                        }
                    }

                    if (targetAlleleSupport >= minIsolatedCoverage) {
                        validChains.push_back({nodeIdx});
                        isAssigned[nodeIdx] = true;
                    }
                }
            }

            // --- Step 5: Mark Valid Clusters & Identify Haplotype Reads ---
            
            // Collect all In-Phase reads (supporting Target Allele) 
            // and Out-Of-Phase reads (supporting Alt Allele) from VALID chains.
            vector<uint64_t> inPhaseReads;
            vector<uint64_t> outOfPhaseReads;

            for (const auto& chain : validChains) {
                for (int nodeIdx : chain) {

                    // Mark cluster as valid in global array
                    uint64_t originalClusterId = dpNodes[nodeIdx].clusterId;

                    // Benign race: multiple threads may write 1 simultaneously.
                    // variantClusteringValidClustersCompatible[originalClusterId] = 1;
                    // Mark compatible atomically to avoid data races
                    __sync_bool_compare_and_swap(&variantClusteringValidClustersCompatible[originalClusterId], 0, 1);

                    // Collect reads
                    for (const auto& rp : nodeReads[nodeIdx]) {
                        if (rp.isPhase1 == false) {
                            // Supported Target Allele -> In Phase
                            inPhaseReads.push_back(rp.orientedReadIdValue);
                        } else {
                            // Supported Alt Allele -> Out of Phase
                            outOfPhaseReads.push_back(rp.orientedReadIdValue);
                        }
                    }
                }
            }

            // Sort & Unique for set operations
            std::sort(inPhaseReads.begin(), inPhaseReads.end());
            inPhaseReads.erase(std::unique(inPhaseReads.begin(), inPhaseReads.end()), inPhaseReads.end());

            std::sort(outOfPhaseReads.begin(), outOfPhaseReads.end());
            outOfPhaseReads.erase(std::unique(outOfPhaseReads.begin(), outOfPhaseReads.end()), outOfPhaseReads.end());

            // Final Set: InPhase MINUS OutOfPhase
            // A read is only trusted if it NEVER supported the Alt allele in any valid cluster.
            
            // Calculate total weight (total number of nodes in all valid chains)
            // This is the number of sites that were compatible and were used to cast votes
            uint32_t totalWeight = 0;
            for (const auto& chain : validChains) {
                totalWeight += (uint32_t)chain.size();
            }

            for (const auto& chain : validChains) {
                for (int nodeIdx : chain) {
                    // Mark cluster as valid... (already done above)

                    // Collect reads
                    for (const auto& rp : nodeReads[nodeIdx]) {
                        if (rp.isPhase1 == false) {
                            // Supported Target Allele -> In Phase
                            inPhaseReads.push_back(rp.orientedReadIdValue);
                        } else {
                            // Supported Alt Allele -> Out of Phase
                            outOfPhaseReads.push_back(rp.orientedReadIdValue);
                        }
                    }
                }
            }

            // Sort & Unique for set operations
            std::sort(inPhaseReads.begin(), inPhaseReads.end());
            inPhaseReads.erase(std::unique(inPhaseReads.begin(), inPhaseReads.end()), inPhaseReads.end());

            std::sort(outOfPhaseReads.begin(), outOfPhaseReads.end());
            outOfPhaseReads.erase(std::unique(outOfPhaseReads.begin(), outOfPhaseReads.end()), outOfPhaseReads.end());

            // Final Set: InPhase MINUS OutOfPhase
            vector<OrientedReadId::Int> finalHapOrientedReadIds;
            // Add Target Read itself (always in phase)
            finalHapOrientedReadIds.push_back(currentOrientedReadId0.getValue());
            std::set_difference(
                inPhaseReads.begin(), inPhaseReads.end(),
                outOfPhaseReads.begin(), outOfPhaseReads.end(),
                std::back_inserter(finalHapOrientedReadIds)
            );

            // Sort again to include the target read
            std::sort(finalHapOrientedReadIds.begin(), finalHapOrientedReadIds.end());
            finalHapOrientedReadIds.erase(std::unique(finalHapOrientedReadIds.begin(), finalHapOrientedReadIds.end()), finalHapOrientedReadIds.end());
            

            // Print finalHapOrientedReadIds with one orientedReadId per line if readId is 0
            if(currentOrientedReadId0.getReadId() == 87 and currentOrientedReadId0.getStrand() == 0) {
                cout << "Final Haplotype Oriented Reads: " << endl;
                for(OrientedReadId::Int orientedReadIdValue : finalHapOrientedReadIds) {
                    cout << OrientedReadId::fromValue(orientedReadIdValue).getString() << endl;
                }
            }


            // CAST VOTES DIRECTLY TO GRAPH
            // Add votes: currentOrientedReadId -> otherOrientedReadId
            {
                std::lock_guard<std::mutex> lock(haplotypeGraphMutex);
                for(OrientedReadId::Int otherOrientedReadIdValue : finalHapOrientedReadIds) {
                    if(otherOrientedReadIdValue != currentOrientedReadId0.getValue()) {
                        // Try to add edge with totalWeight
                        auto result = boost::add_edge(currentOrientedReadId0.getValue(), otherOrientedReadIdValue, totalWeight, *globalHaplotypeGraph);
                    }
                }
            }


        }
    }
}


void Assembler::refineClustersThreadFunction(uint64_t threadId) {
    const uint64_t currentMinAlleleCoverage = this->minAlleleCoverage;
    const auto& graph = *globalHaplotypeGraph;
    
    uint64_t clusterIdBegin, clusterIdEnd;
    while (getNextBatch(clusterIdBegin, clusterIdEnd)) {
        for (uint64_t clusterId = clusterIdBegin; clusterId < clusterIdEnd; clusterId++) {
            
            // Only process valid/compatible clusters
            if (!variantClusteringValidClustersCompatible[clusterId]) {
                continue;
            }

            const auto& members = variantClusteringMembersByRepIdx[clusterId];
            
            // 1. Identify Valid Alleles (those with sufficient coverage)
            std::array<uint32_t, 5> alleleCounts = {0};
            for (uint64_t memberIdx : members) {
                if (memberIdx < variantClusteringPositionPairAlleles.size()) {
                    uint8_t a = variantClusteringPositionPairAlleles[memberIdx];
                    if (a < 5) alleleCounts[a]++;
                }
            }
            
            std::vector<uint8_t> validAlleles;
            for(uint8_t a=0; a<5; a++) {
                if (alleleCounts[a] >= currentMinAlleleCoverage) {
                    validAlleles.push_back(a);
                }
            }
            
            // If fewer than 2 valid alleles, no phasing to refine
            if (validAlleles.size() < 2) continue;

            // 2. Partition Reads by Allele using oriented read IDs
            std::array<std::unordered_set<uint32_t>, 5> groups;
            for (uint64_t memberIdx : members) {
                if (memberIdx < variantClusteringPositionPairAlleles.size()) {
                    uint8_t a = variantClusteringPositionPairAlleles[memberIdx];
                    if (a < 5) {
                        const auto& pp = variantClusteringPositionPairs[memberIdx];
                        groups[a].insert(pp.first.getValue());
                    }
                }
            }

            // 3. Process each allele group separately
            for (uint8_t allele = 0; allele < 5; allele++) {
                const auto& group = groups[allele];
                
                // If group is too small, we can't really do intersection of 2 reads + others.
                // If size < 3, we just keep them (default status 0).
                if(group.size() < 3) {
                    continue;
                }

                // Store (ReadVal, Weight) pairs
                std::vector<std::pair<uint32_t, uint32_t>> readWeights;
                readWeights.reserve(group.size());

                // Access weight map
                auto weightMap = boost::get(boost::edge_weight, graph);

                for(uint32_t readVal : group) {
                    uint32_t weight = 0;
                    // Check if vertex exists in graph
                    if (readVal >= boost::num_vertices(graph)) {
                        readWeights.push_back({readVal, 0}); 
                    } else {
                        // Get weight from first outgoing edge
                        auto outEdges = boost::out_edges(readVal, graph);
                        if(outEdges.first != outEdges.second) {
                            auto edge = *outEdges.first;
                            weight = weightMap[edge];
                        }
                    }
                    readWeights.push_back({readVal, weight});
                }

                // Sort by Weight Descending
                std::sort(readWeights.begin(), readWeights.end(), 
                    [](const std::pair<uint32_t, uint32_t>& a, const std::pair<uint32_t, uint32_t>& b) {
                        return a.second > b.second;
                    });

                // If the top 2 have 0 weight (no edges), keep all.
                if (readWeights[0].second == 0 || readWeights.size() < 2) {
                    continue;
                }

                // Select Top 2
                uint32_t r1 = readWeights[0].first;
                uint32_t r2 = readWeights[1].first;
                
                // Filter
                for(const auto& rw : readWeights) {
                    uint32_t readVal = rw.first;
                    
                    // Keep Top 2
                    if (readVal == r1 || readVal == r2) {
                        continue; 
                    }

                    // Check connectivity to R1 and R2
                    bool connectedToR1 = false;
                    bool connectedToR2 = false;

                    if (r1 < boost::num_vertices(graph) && r2 < boost::num_vertices(graph) && readVal < boost::num_vertices(graph)) {
                        
                        // Check R1 neighbors
                        auto neighborsR1 = boost::adjacent_vertices(r1, graph);
                        for(auto it = neighborsR1.first; it != neighborsR1.second; ++it) {
                            if(*it == readVal) { connectedToR1 = true; break; }
                        }
                        if (!connectedToR1) {
                            auto inEdgesR1 = boost::in_edges(r1, graph);
                            for(auto it = inEdgesR1.first; it != inEdgesR1.second; ++it) {
                                if(boost::source(*it, graph) == readVal) { connectedToR1 = true; break; }
                            }
                        }

                        // Check R2 neighbors
                        auto neighborsR2 = boost::adjacent_vertices(r2, graph);
                        for(auto it = neighborsR2.first; it != neighborsR2.second; ++it) {
                            if(*it == readVal) { connectedToR2 = true; break; }
                        }
                        if (!connectedToR2) {
                            auto inEdgesR2 = boost::in_edges(r2, graph);
                            for(auto it = inEdgesR2.first; it != inEdgesR2.second; ++it) {
                                if(boost::source(*it, graph) == readVal) { connectedToR2 = true; break; }
                            }
                        }
                    }

                    if (connectedToR1 && connectedToR2) {
                        // Keep
                    } else {
                        // Mark as Stray
                        // Need to find memberIdx
                        for(uint64_t memberIdx : members) {
                            if (memberIdx < variantClusteringPositionPairAlleles.size()) {
                                if (variantClusteringPositionPairAlleles[memberIdx] == allele && 
                                    variantClusteringPositionPairs[memberIdx].first.getValue() == readVal) {
                                    variantClusteringMemberStatus[memberIdx] = 1; // Mark as Stray
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
