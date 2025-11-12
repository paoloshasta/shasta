// Shasta.
#include "Assembler.hpp"
#include "ProjectedAlignment.hpp"
#include "SHASTA_ASSERT.hpp"
#include "iostream.hpp"
#include "Reads.hpp"
using namespace shasta;

// Helper function called during alignment computation in AssemblerAlign.cpp to collect position pairs
// that have SNP differences
void Assembler::collectVariantClusteringPositionPairs(
    const ProjectedAlignment& projectedAlignment,
    const array<OrientedReadId, 2>& orientedReadIds,
    MemoryMapped::Vector< pair<OrientedReadId, uint32_t> >& positionPairs)
{
    const OrientedReadId orientedReadId0 = orientedReadIds[0];
    const OrientedReadId orientedReadId1 = orientedReadIds[1];
    
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

                    // Get read IDs and lengths for position normalization
                    const ReadId readId0 = orientedReadId0.getReadId();
                    const ReadId readId1 = orientedReadId1.getReadId();
                    const uint64_t readLength0 = getReads().getReadRawSequenceLength(readId0);
                    const uint64_t readLength1 = getReads().getReadRawSequenceLength(readId1);
                    
                    // Normalize positions to strand 0 coordinates
                    uint64_t position0Strand0 = positionInRead0;
                    uint64_t position1Strand0 = positionInRead1;
                    if (orientedReadId0.getStrand() == 1) {
                        position0Strand0 = readLength0 - 1 - positionInRead0;
                    }
                    if (orientedReadId1.getStrand() == 1) {
                        position1Strand0 = readLength1 - 1 - positionInRead1;
                    }
                    
                    // Create pairs for all four strand combinations
                    // OrientedReadId0-0, OrientedReadId0-1, OrientedReadId1-0, OrientedReadId1-1
                    positionPairs.push_back(make_pair(OrientedReadId(readId0, 0), uint32_t(position0Strand0)));
                    positionPairs.push_back(make_pair(OrientedReadId(readId0, 1), uint32_t(readLength0 - 1 - position0Strand0)));
                    positionPairs.push_back(make_pair(OrientedReadId(readId1, 0), uint32_t(position1Strand0)));
                    positionPairs.push_back(make_pair(OrientedReadId(readId1, 1), uint32_t(readLength1 - 1 - position1Strand0)));

                    if ( (orientedReadId0 == OrientedReadId(0, 0) and orientedReadId1 == OrientedReadId(453, 1)) or (orientedReadId0 == OrientedReadId(453, 1) and orientedReadId1 == OrientedReadId(0, 0))) {
                        // print the position pairs
                        cout << "Position pair: " << OrientedReadId(readId0, 0) << " " << position0Strand0 << " " << OrientedReadId(readId1, 0) << " " << position1Strand0 << endl;
                        cout << "Position pair: " << OrientedReadId(readId0, 1) << " " << readLength0 - 1 - position0Strand0 << " " << OrientedReadId(readId1, 1) << " " << readLength1 - 1 - position1Strand0 << endl; 
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