#include "Assembler.hpp"
#include "timestamp.hpp"
using namespace shasta;



void Assembler::createReadGraph5()
{
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

    cout << timestamp << "createReadGraph5 ends." << endl;
}



void Assembler::createReadGraph5ThreadFunction(uint64_t /* threadId */)
{
    SHASTA_ASSERT(0);
}
