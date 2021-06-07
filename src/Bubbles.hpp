#ifndef SHASTA_BUBBLES_HPP
#define SHASTA_BUBBLES_HPP

/*******************************************************************************

Class to describe an analyze bubbles in the assembly graph.

*******************************************************************************/

#include "AssemblyGraph.hpp"

#include "vector.hpp"

namespace shasta {
    class Bubbles;
    class Assembler;
}

class shasta::Bubbles {
public:

    Bubbles(
        const Assembler&
    );

private:

    // For now we only consider diploid bubbles, defined using the
    // following strict criteria:
    // - Source vertex v0 has out-degree 2.
    // - Target vertex v1 has in-degree 2.
    // - There are two parallel edges eA and eB, both v0->v1.
    class Bubble {
    public:
        AssemblyGraph::VertexId v0;
        AssemblyGraph::VertexId v1;
        AssemblyGraph::EdgeId eA;
        AssemblyGraph::EdgeId eB;
    };
    vector<Bubble> bubbles;
    void findBubbles();

    // Figure out if two sequences differ only by copy numbers in
    // a 2- or 3-base repeat.
    static bool isShortRepeatCopyNumberDifference(
        const vector<Base>&,
        const vector<Base>&);

    const Assembler& assembler;

};

#endif
