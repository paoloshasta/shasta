// Shasta.
#include "mode3-AssemblyGraph.hpp"
#include "mode3-LocalAssembly.hpp"
#include "mode3-PrimaryGraph.hpp"
#include "mode3-PhasingTable.hpp"
#include "Assembler.hpp"
#include "copyNumber.hpp"
#include "deduplicate.hpp"
#include "diploidBayesianPhase.hpp"
#include "dominatorTree.hpp"
#include "enumeratePaths.hpp"
#include "findLinearChains.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/graph/strong_components.hpp>

// Standard library.
#include "fstream.hpp"
#include <queue>
#include "tuple.hpp"

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<AssemblyGraph>;


// Create from a connected component of the PrimaryGraph, then call run.
AssemblyGraph::AssemblyGraph(
    const PrimaryGraph& graph,
    uint64_t componentId,
    const Assembler& assembler,
    uint64_t threadCount,
    bool debug) :
    MultithreadedObject<AssemblyGraph>(*this),
    componentId(componentId),
    assembler(assembler)
{
    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    performanceLog << timestamp << "Creating the assembly graph for component " << componentId << endl;
    create(graph, debug);

    // Serialize it so we can restore it to facilitate debugging.
    save("AssemblyGraph-" + to_string(componentId) + ".data");

    performanceLog << timestamp << "Processing the assembly graph for component " << componentId << endl;
    run(threadCount, true, debug);
    performanceLog << timestamp << "Done with the assembly graph for component " << componentId << endl;
}



// Load it from a binary archive, then call run.
AssemblyGraph::AssemblyGraph(
    const string& fileName,
    const Assembler& assembler,
    uint64_t threadCount,
    bool debug) :
    MultithreadedObject<AssemblyGraph>(*this),
    assembler(assembler)
{
    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    load(fileName);
    run(threadCount, true, debug);
}



void AssemblyGraph::run(
    uint64_t threadCount,
    bool assembleSequence,
    bool debug)
{
    // *** EXPOSE WHEN CODE STABILIZES
    const uint64_t detangleToleranceLow = 0;
    const uint64_t detangleToleranceHigh = 2;
    // const uint64_t detangleWithSearchToleranceLow = 1;
    // const uint64_t detangleWithSearchToleranceHigh = 6;
    const bool useBayesianModel = true;
    const double epsilon = 0.1;
    const double minLogP = 20.;
    const uint64_t longBubbleThreshold = 5000;
    const double phaseErrorThreshold = 0.1;
    const double bubbleErrorThreshold = 0.03;
    const uint64_t chainTerminalCommonThreshold = 3;
    // const uint64_t optimizeChainsMinCommon = 3;
    // const uint64_t optimizeChainsK = 100;

    if(debug) write("A");

    // Don't do any detangling before cleanup of bubbles and superbubbles.

    // Cleanup bubbles and superbubbles.
    // Must do compress to make sure all bubbles are in bubble chains.
    compress();
    for(uint64_t iteration=0; ; iteration ++) {
        performanceLog << timestamp << "Iteration " << iteration <<
            " of bubble cleanup begins." << endl;
        const uint64_t cleanedUpBubbleCount = cleanupBubbles(
            false, 1000, chainTerminalCommonThreshold, threadCount);
        if(cleanedUpBubbleCount == 0) {
            break;
        }
        cout << "Cleaned up " << cleanedUpBubbleCount << " bubbles probably caused by errors." << endl;
        compressBubbleChains();
        compress();
    }
    if(debug) write("B");
    cleanupSuperbubbles(false, 30000, chainTerminalCommonThreshold);
    compress();

    // Remove short superbubbles.
    removeShortSuperbubbles(false, 10000, 30000);
    compress();

    // Phase.
    compressBubbleChains();
    if(debug) write("C");
    phaseBubbleChainsUsingPhasingTable(
        debug ? "C" : "",
        phaseErrorThreshold,
        bubbleErrorThreshold,
        longBubbleThreshold);
    compress();

    // For detangling, expand all bubble chains.
    expand();

    // Detangle.
    if(debug) write("D");
    performanceLog << timestamp << "Detangling begins." << endl;
    while(compressSequentialEdges());
    compressBubbleChains();
    detangleEdges(false, detangleToleranceLow, detangleToleranceHigh, useBayesianModel, epsilon, minLogP);
    while(compressSequentialEdges());
    compressBubbleChains();
    detangleVertices(false, detangleToleranceLow, detangleToleranceHigh, useBayesianModel, epsilon, minLogP);
    while(compressSequentialEdges());
    compressBubbleChains();
    detangleEdges(false, detangleToleranceLow, detangleToleranceHigh, useBayesianModel, epsilon, minLogP);
    // detangleShortSuperbubbles(false, 30000, detangleToleranceLow, detangleToleranceHigh);
    performanceLog << timestamp << "Detangling ends." << endl;

    compress();
    compressBubbleChains();
    if(debug) write("E");

#if 0
    // Optimize the chains.
    optimizeChains(
        false,
        optimizeChainsMinCommon,
        optimizeChainsK);
#endif

    // Before final output, renumber the edges contiguously.
    renumberEdges();
    if(debug) write("F");

    if(assembleSequence) {

        // Assemble sequence.
        assembleAllChainsMultithreaded(chainTerminalCommonThreshold, threadCount);
        writeAssemblyDetails();

        if(debug) write("G", true);

    } else {

        // Skip sequence assembly.
        if(debug) write("G");
    }


}



// Initial creation from the PrimaryGraph.
// Each linear chain of edges in the PrimaryGraph after transitive reduction generates
// an AssemblyGraphEdge (BubbleChain) consisting of a single haploid bubble.
void AssemblyGraph::create(const PrimaryGraph& graph, bool debug)
{
    AssemblyGraph& cGraph = *this;

    // Create a filtered version of the PathGraph, containing only the
    // transitive reduction edges.
    class EdgePredicate {
    public:
        bool operator()(const PrimaryGraph::edge_descriptor e) const
        {
            return not (*graph)[e].isNonTransitiveReductionEdge;
        }
        EdgePredicate(const PrimaryGraph& graph) : graph(&graph) {}
        EdgePredicate() : graph(0) {}
    private:
        const PrimaryGraph* graph;
    };
    using FilteredPrimaryGraph = boost::filtered_graph<PrimaryGraph, EdgePredicate>;
    FilteredPrimaryGraph filteredGraph(graph, EdgePredicate(graph));

    // Find linear chains in the PathGraph after transitive reduction.
    vector< vector<PrimaryGraph::edge_descriptor> > inputChains;
    findLinearChains(filteredGraph, 0, inputChains);

    // Each chain generates an edge.
    // Vertices are added as needed.
    std::map<MarkerGraphEdgeId, vertex_descriptor> vertexMap;
    for(const vector<PrimaryGraph::edge_descriptor>& inputChain: inputChains) {
        const PrimaryGraph::vertex_descriptor v0 = source(inputChain.front(), graph);
        const PrimaryGraph::vertex_descriptor v1 = target(inputChain.back(), graph);
        const MarkerGraphEdgeId markerGraphEdgeId0 = graph[v0].edgeId;
        const MarkerGraphEdgeId markerGraphEdgeId1 = graph[v1].edgeId;
        const vertex_descriptor cv0 = getVertex(markerGraphEdgeId0, vertexMap);
        const vertex_descriptor cv1 = getVertex(markerGraphEdgeId1, vertexMap);

        // Create an edge for this input chain.
        edge_descriptor ce;
        tie(ce, ignore) = add_edge(cv0, cv1, cGraph);
        AssemblyGraphEdge& edge = cGraph[ce];
        edge.id = nextEdgeId++;

        // The edge is a degenerate BubbleChain consisting of a single haploid bubble.
        edge.resize(1);                 // BubbleChain has length 1.
        Bubble& bubble = edge.front();
        bubble.resize(1);               // Bubble is haploid.

        // Store the chain.
        Chain& chain = bubble.front();
        for(const PrimaryGraph::edge_descriptor e: inputChain) {
            const PrimaryGraph::vertex_descriptor v = source(e, graph);
            chain.push_back(graph[v].edgeId);
        }
        const PrimaryGraph::edge_descriptor eLast = inputChain.back();
        const PrimaryGraph::vertex_descriptor vLast = target(eLast, graph);
        chain.push_back(graph[vLast].edgeId);
    }
}



// Return the vertex corresponding to a given MarkerGraphEdgeId,
// creating it if it is not in the given vertexMap
AssemblyGraph::vertex_descriptor AssemblyGraph::getVertex(
    MarkerGraphEdgeId markerGraphEdgeId,
    std::map<MarkerGraphEdgeId, vertex_descriptor>& vertexMap)
{
    AssemblyGraph& cGraph = *this;

    auto it = vertexMap.find(markerGraphEdgeId);
    if(it == vertexMap.end()) {
        const vertex_descriptor cv = add_vertex({markerGraphEdgeId}, cGraph);
        vertexMap.insert({markerGraphEdgeId, cv});
        return cv;
    } else {
        return it->second;
    }
}



// Create a new vertex with a given MarkerGraphEdgeId.
AssemblyGraph::vertex_descriptor AssemblyGraph::createVertex(
    MarkerGraphEdgeId markerGraphEdgeId)
{
    return add_vertex({markerGraphEdgeId}, *this);
}



void AssemblyGraph::removeVertex(vertex_descriptor cv)
{
    AssemblyGraph& cGraph = *this;

    SHASTA_ASSERT(in_degree(cv, cGraph) == 0);
    SHASTA_ASSERT(out_degree(cv, cGraph) == 0);

    boost::remove_vertex(cv, cGraph);
}



// Compute vertexIndex for every vertex.
// This numbers vertices consecutively starting at zero.
// This numbering becomes invalid as soon as a vertex is added or removed.
void AssemblyGraph::numberVertices()
{
    AssemblyGraph& cGraph = *this;
    uint64_t index = 0;
    BGL_FORALL_VERTICES(cv, cGraph, AssemblyGraph) {
        cGraph[cv].index = index++;
    }
}



void AssemblyGraph::clearVertexNumbering()
{
    AssemblyGraph& cGraph = *this;
    BGL_FORALL_VERTICES(cv, cGraph, AssemblyGraph) {
        cGraph[cv].index = invalid<uint64_t>;
    }

}


void AssemblyGraph::renumberEdges()
{
    AssemblyGraph& cGraph = *this;
    nextEdgeId = 0;

    BGL_FORALL_EDGES(ce, cGraph, AssemblyGraph) {
        cGraph[ce].id = nextEdgeId++;
    }
}



// Compress parallel edges into bubbles, where possible.
bool AssemblyGraph::compressParallelEdges()
{
    AssemblyGraph& cGraph = *this;
    bool changesWereMade = false;

    // Look for sets of parallel edges v0->v1.
    vector<vertex_descriptor> childrenVertices;
    vector<edge_descriptor> edgesToBeRemoved;
    Bubble newBubble;
    BGL_FORALL_VERTICES(v0, cGraph, AssemblyGraph) {
        if(out_degree(v0, cGraph) < 2) {
            continue;
        }

        // Find distinct children vertices of v0.
        childrenVertices.clear();
        BGL_FORALL_OUTEDGES(v0, e, cGraph, AssemblyGraph) {
            childrenVertices.push_back(target(e, cGraph));
        }
        deduplicate(childrenVertices);

        // Handle the children vertices one at a time.
        for(const vertex_descriptor v1: childrenVertices) {

            // Create the new bubble using parallel edges v0->v1.
            newBubble.clear();
            edgesToBeRemoved.clear();
            BGL_FORALL_OUTEDGES(v0, e, cGraph, AssemblyGraph) {
                if(target(e, cGraph) != v1) {
                    continue;
                }
                AssemblyGraphEdge& edge = cGraph[e];

                // The BubbleChain must have length 1.
                if(edge.size() > 1) {
                    continue;
                }
                const Bubble& oldBubble = edge.front();

                copy(oldBubble.begin(), oldBubble.end(), back_inserter(newBubble));
                edgesToBeRemoved.push_back(e);
            }
            if(edgesToBeRemoved.size() < 2) {
                continue;
            }

            // Create the new edge.
            changesWereMade = true;
            edge_descriptor eNew;
            tie(eNew, ignore) = add_edge(v0, v1, cGraph);
            AssemblyGraphEdge& newEdge = cGraph[eNew];
            newEdge.id = nextEdgeId++;
            newEdge.resize(1);  // Make it a single bubble.
            Bubble& newEdgeBubble = newEdge.front();
            newEdgeBubble = newBubble;
            newEdgeBubble.deduplicate();

            // Remove the old edges.
            for(const edge_descriptor e: edgesToBeRemoved) {
                boost::remove_edge(e, cGraph);
            }

        }
    }
    return changesWereMade;
}



// Remove duplicate chains.
void Bubble::deduplicate()
{
    shasta::deduplicate(*this);
}



// Compress linear sequences of edges (BubbleChains) into longer BubbleChains.
bool AssemblyGraph::compressSequentialEdges()
{
    AssemblyGraph& cGraph = *this;
    bool changesWereMade = false;

    // Find linear chains of edges.
    vector< vector<edge_descriptor> > linearChains;
    findLinearChains(cGraph, 0, linearChains);



    // Each linear chain of more than one edge gets compressed into a single edge (BubbleChain).
    for(const vector<edge_descriptor>& linearChain: linearChains) {
        if(linearChain.size() < 2) {
            continue;
        }

        // Create the new edge.
        changesWereMade = true;
        const vertex_descriptor v0 = source(linearChain.front(), cGraph);
        const vertex_descriptor v1 = target(linearChain.back(), cGraph);
        edge_descriptor ceNew;
        tie(ceNew, ignore) = add_edge(v0, v1, cGraph);
        AssemblyGraphEdge& newEdge = cGraph[ceNew];
        newEdge.id = nextEdgeId++;
        for(const edge_descriptor ce: linearChain) {
            const AssemblyGraphEdge& oldEdge = cGraph[ce];
            copy(oldEdge.begin(), oldEdge.end(), back_inserter(newEdge));
        }

        // Remove the old edges.
        for(const edge_descriptor ce: linearChain) {
            boost::remove_edge(ce, cGraph);
        }

        // Remove the vertices internal to the old edge.
        for(uint64_t i=1; i<linearChain.size(); i++) {
            const vertex_descriptor cv = source(linearChain[i], cGraph);
            cGraph.removeVertex(cv);
        }
    }
    return changesWereMade;
}



// Call compressParallelEdges and compressSequentialEdges iteratively until nothing changes.
bool AssemblyGraph::compress()
{
    bool changesWereMade = false;

    while(true) {
        const bool compressParallelChanges = compressParallelEdges();
        const bool compressSequentialChanges = compressSequentialEdges();

        if(compressParallelChanges or compressSequentialChanges) {
            // Something changed. Continue the iteration loop.
            changesWereMade = true;
            continue;
        } else {
            // Nothing changed at this iteration. Stop iteration loop.
            break;
        }
    }

    return changesWereMade;
}



// Call compress on all BubbleChains to merge adjacent haploid bubbles.
void AssemblyGraph::compressBubbleChains()
{
    AssemblyGraph& cGraph = *this;
    BGL_FORALL_EDGES(e, cGraph, AssemblyGraph) {
        cGraph[e].compress();
    }
}



// This does the opposite of compress. All bubble chains that
// consist of more than one simple haploid bubble are expanded into one
// edge for each edge of each bubble.
// For optimal results it is best to call compressBubbleChains before expand.
void AssemblyGraph::expand()
{
    AssemblyGraph& cGraph = *this;

    // Gather all edges that exist at this point.
    vector<edge_descriptor> initialEdges;
    BGL_FORALL_EDGES(ce, cGraph, AssemblyGraph) {
        initialEdges.push_back(ce);
    }



    // Loop over the initial edges.
    for(const edge_descriptor ce: initialEdges) {
        BubbleChain& bubbleChain = cGraph[ce];

        // If this bubbleChain consists of a single haploid bubble, don't do anything.
        if(bubbleChain.isSimpleChain()) {
            continue;
        }

        // Prepare a vector of the vertices that will be the sources and targets
        // of the edges we will create.
        vector<vertex_descriptor> newVertices;
        newVertices.push_back(source(ce, cGraph));
        for(uint64_t positionInBubbleChain=1; positionInBubbleChain<bubbleChain.size(); positionInBubbleChain++) {
            const vertex_descriptor cv = createVertex(bubbleChain[positionInBubbleChain].front().front());
            newVertices.push_back(cv);
        }
        newVertices.push_back(target(ce, cGraph));

        // Create a new edge for each chain of each bubble in this bubble chain.
        for(uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size(); positionInBubbleChain++) {
            Bubble& bubble = bubbleChain[positionInBubbleChain];
            const vertex_descriptor cv0 = newVertices[positionInBubbleChain];
            const vertex_descriptor cv1 = newVertices[positionInBubbleChain + 1];

            for(uint64_t indexInBubble=0; indexInBubble<bubble.size(); indexInBubble++) {
                Chain& chain = bubble[indexInBubble];

                // Create a new edge for this chain.
                edge_descriptor ceNew;
                tie(ceNew, ignore) = add_edge(cv0, cv1, cGraph);
                AssemblyGraphEdge& edge = cGraph[ceNew];
                edge.id = nextEdgeId++;

                // Store this Chain in the new edge.
                BubbleChain& newBubbleChain = cGraph[ceNew];
                newBubbleChain.resize(1);
                Bubble& newBubble = newBubbleChain.front();
                newBubble.resize(1);
                Chain& newChain = newBubble.front();
                newChain.swap(chain);
            }
        }

        // Now we can remove the BubbleChain.
        boost::remove_edge(ce, cGraph);
    }
}



void AssemblyGraph::write(const string& name, bool writeSequence) const
{
    const string fileNamePrefix = name + "-" + to_string(componentId);

    cout << fileNamePrefix << ": " << num_vertices(*this) <<
        " vertices, " << num_edges(*this) << " edges. Next edge id " << nextEdgeId << endl;

    writeCsv(fileNamePrefix);
    writeGraphviz(fileNamePrefix, true);
    writeGraphviz(fileNamePrefix, false);
    writeGfa(fileNamePrefix);
    writeGfaExpanded(name, writeSequence);
    if(writeSequence) {
        writeFastaExpanded(name);
    }
}



void AssemblyGraph::writeCsv(const string& fileNamePrefix) const
{
    writeChainsDetailsCsv(fileNamePrefix);
    writeChainsCsv(fileNamePrefix);
    writeBubblesCsv(fileNamePrefix);
    writeBubbleChainsCsv(fileNamePrefix);
}



void AssemblyGraph::writeBubbleChainsCsv(const string& fileNamePrefix) const
{
    const AssemblyGraph& cGraph = *this;

    ofstream csv(fileNamePrefix + "-BubbleChains.csv");
    csv << "Id,ComponentId,BubbleChainId,v0,v1,BubbleCount,AverageOffset,MinOffset,MaxOffset,\n";

    BGL_FORALL_EDGES(ce, cGraph, AssemblyGraph) {
        const vertex_descriptor cv0 = source(ce, cGraph);
        const vertex_descriptor cv1 = target(ce, cGraph);
        const BubbleChain& bubbleChain = cGraph[ce];

        uint64_t averageOffset;
        uint64_t minOffset;
        uint64_t maxOffset;
        bubbleChainOffset(bubbleChain, averageOffset, minOffset, maxOffset);

        csv << bubbleChainStringId(ce) << ",";
        csv << componentId << ",";
        csv << cGraph[ce].id << ",";
        csv << cGraph[cv0].edgeId << ",";
        csv << cGraph[cv1].edgeId << ",";
        csv << bubbleChain.size() << ",";
        csv << averageOffset << ",";
        csv << minOffset << ",";
        csv << maxOffset << ",";
        csv << "\n";
    }
}




void AssemblyGraph::writeBubblesCsv(const string& fileNamePrefix) const
{
    const AssemblyGraph& cGraph = *this;

    ofstream csv(fileNamePrefix + "-Bubbles.csv");
    csv << "Id,ComponentId,BubbleChainId,Position in bubble chain,v0,v1,Ploidy,AverageOffset,MinOffset,MaxOffset,\n";

    BGL_FORALL_EDGES(ce, cGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = cGraph[ce];

        for(uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size(); positionInBubbleChain++) {
            const Bubble& bubble = bubbleChain[positionInBubbleChain];
            const Chain& firstChain = bubble.front();

            // Check that all the chains begins/end in the same place.
            for(const Chain& chain: bubble) {
                SHASTA_ASSERT(chain.front() == firstChain.front());
                SHASTA_ASSERT(chain.back() == firstChain.back());
            }

            uint64_t averageOffset;
            uint64_t minOffset;
            uint64_t maxOffset;
            bubbleOffset(bubble, averageOffset, minOffset, maxOffset);

            csv << bubbleStringId(ce, positionInBubbleChain) << ",";
            csv << componentId << ",";
            csv << cGraph[ce].id << ",";
            csv << positionInBubbleChain << ",";
            csv << firstChain.front() << ",";
            csv << firstChain.back() << ",";
            csv << bubble.size() << ",";
            csv << averageOffset << ",";
            csv << minOffset << ",";
            csv << maxOffset << ",";
            csv << "\n";
        }
    }

}


void AssemblyGraph::writeChainsCsv(const string& fileNamePrefix) const
{
    const AssemblyGraph& cGraph = *this;

    ofstream csv(fileNamePrefix + "-Chains.csv");
    csv << "Id,ComponentId,BubbleChainId,Position in bubble chain,Index in bubble,Length,Offset\n";

    BGL_FORALL_EDGES(ce, cGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = cGraph[ce];

        for(uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size(); positionInBubbleChain++) {
            const Bubble& bubble = bubbleChain[positionInBubbleChain];
            const uint64_t ploidy = bubble.size();

            for(uint64_t indexInBubble=0; indexInBubble<ploidy; indexInBubble++) {
                const Chain& chain = bubble[indexInBubble];
                SHASTA_ASSERT(chain.size() >= 2);

                csv << chainStringId(ce, positionInBubbleChain, indexInBubble) << ",";
                csv << componentId << ",";
                csv << cGraph[ce].id << ",";
                csv << positionInBubbleChain << ",";
                csv << indexInBubble << ",";
                csv << chain.size() << ",";
                csv << chainOffset(chain) << ",";
                csv << "\n";
            }
        }
    }

}



void AssemblyGraph::writeChainsDetailsCsv(const string& fileNamePrefix) const
{
    const AssemblyGraph& cGraph = *this;

    ofstream csv(fileNamePrefix + "-ChainsDetails.csv");
    csv << "Id,ComponentId,BubbleChainId,Position in bubble chain,"
        "Index in bubble,Position in chain,MarkerGraphEdgeId,Coverage,Common,Offset\n";

    BGL_FORALL_EDGES(e, cGraph, AssemblyGraph) {
        writeChainDetailsCsv(csv, e, false);
    }
}



void AssemblyGraph::writeChainDetailsCsv(
    ostream& csv,
    edge_descriptor e,
    bool writeHeader) const
{
    const AssemblyGraph& cGraph = *this;
    const BubbleChain& bubbleChain = cGraph[e];

    if(writeHeader) {
        csv << "Id,ComponentId,BubbleChainId,Position in bubble chain,"
            "Index in bubble,Position in chain,MarkerGraphEdgeId,Coverage,Common,Offset\n";
    }

    for(uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size(); positionInBubbleChain++) {
        const Bubble& bubble = bubbleChain[positionInBubbleChain];
        const uint64_t ploidy = bubble.size();

        for(uint64_t indexInBubble=0; indexInBubble<ploidy; indexInBubble++) {
            const Chain& chain = bubble[indexInBubble];
            SHASTA_ASSERT(chain.size() >= 2);

            for(uint64_t positionInChain=0; positionInChain<chain.size(); positionInChain++) {
                const MarkerGraphEdgeId markerGraphEdgeId = chain[positionInChain];
                const uint64_t coverage = assembler.markerGraph.edgeCoverage(markerGraphEdgeId);
                csv << chainStringId(e, positionInBubbleChain, indexInBubble) << ",";
                csv << componentId << ",";
                csv << cGraph[e].id << ",";
                csv << positionInBubbleChain << ",";
                csv << indexInBubble << ",";
                csv << positionInChain << ",";
                csv << markerGraphEdgeId << ",";
                csv << coverage << ",";

                if(positionInChain != 0) {
                    const MarkerGraphEdgeId previousMarkerGraphEdgeId = chain[positionInChain - 1];
                    MarkerGraphEdgePairInfo info;
                    SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(
                        previousMarkerGraphEdgeId, markerGraphEdgeId, info));
                    csv << info.common << ",";
                    if(info.common != 0) {
                        csv << info.offsetInBases << ",";
                    }
                }
                csv << "\n";
            }
        }
    }
}



void AssemblyGraph::writeGraphviz(
    const string& fileNamePrefix,
    bool labels) const
{
    const AssemblyGraph& cGraph = *this;

    ofstream dot;
    if(labels) {
        dot.open(fileNamePrefix + ".dot");
    } else {
        dot.open(fileNamePrefix + "-NoLabels.dot");
    }

    dot << "digraph Component_" << componentId << "{\n";

    BGL_FORALL_VERTICES(cv, cGraph, AssemblyGraph) {
        const MarkerGraphEdgeId edgeId = cGraph[cv].edgeId;
        const uint64_t coverage = assembler.markerGraph.edgeCoverage(edgeId);
        dot << edgeId << "[label=\"" << edgeId << "\\n" << coverage << "\"];\n";
    }



    BGL_FORALL_EDGES(ce, cGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = cGraph[ce];
        const vertex_descriptor cv0 = source(ce, cGraph);
        const vertex_descriptor cv1 = target(ce, cGraph);

        uint64_t averageOffset;
        uint64_t minOffset;
        uint64_t maxOffset;
        bubbleChainOffset(cGraph[ce], averageOffset, minOffset, maxOffset);

        dot << cGraph[cv0].edgeId << "->" << cGraph[cv1].edgeId;

        if(labels) {
            dot << " [label=\"";
            dot << bubbleChainStringId(ce) << "\\noff=" << averageOffset;

            // Additional annotation if this BubbleChain consists of a single
            // haploid bubble.
            const uint64_t bubbleCount = bubbleChain.size();
            if(bubbleCount == 1) {
                const Bubble& bubble = bubbleChain.front();
                const uint64_t ploidy = bubble.size();
                if(ploidy == 1) {
                    const Chain& chain = bubble.front();
                    dot << "\\nlen=" << chain.size();
                    if(chain.size() > 2) {
                        // Compute average coverage for the internal edges.
                        uint64_t coverageSum = 0;
                        for(uint64_t i=1; i<chain.size()-1; i++) {
                            coverageSum += assembler.markerGraph.edgeCoverage(chain[i]);
                        }
                        const double averageCoverage = double(coverageSum) / double(chain.size() - 2);
                        dot << "\\ncov=" << uint64_t(std::round(averageCoverage));

                        dot << "\\n" << chain.second();
                        if(chain.size() > 3) {
                            dot << "\\n" << chain.secondToLast();
                        }
                    }
                }
            }

            dot << "\"]";
        }
        dot << ";\n";
    }

    dot << "}\n";
}



void AssemblyGraph::writeGfa(const string& fileNamePrefix) const
{
    const AssemblyGraph& cGraph = *this;

    ofstream gfa(fileNamePrefix + ".gfa");

    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";

    // Write a segment for each edge.
    BGL_FORALL_EDGES(ce, cGraph, AssemblyGraph) {

        uint64_t averageOffset;
        uint64_t minOffset;
        uint64_t maxOffset;
        bubbleChainOffset(cGraph[ce], averageOffset, minOffset, maxOffset);

        // Record type.
        gfa << "S\t";

        // Name.
        gfa << bubbleChainStringId(ce) << "\t";

        // Sequence.
        gfa << "*\t";

        // Sequence length in bases.
        gfa << "LN:i:" << averageOffset << "\n";
    }

    // For each vertex, write links between each pair of incoming/outgoing edges.
    BGL_FORALL_VERTICES(cv, cGraph, AssemblyGraph) {
        BGL_FORALL_INEDGES(cv, ceIn, cGraph, AssemblyGraph) {
            BGL_FORALL_OUTEDGES(cv, ceOut, cGraph, AssemblyGraph) {
                gfa <<
                    "L\t" <<
                    bubbleChainStringId(ceIn) << "\t+\t" <<
                    bubbleChainStringId(ceOut) << "\t+\t*\n";
            }
        }
    }
}



void AssemblyGraph::writeGfaExpanded(ostream& gfa, bool includeSequence) const
{
    writeGfaHeader(gfa);
    writeGfaSegmentsExpanded(gfa, includeSequence);
    writeGfaLinksExpanded(gfa);
}



void AssemblyGraph::writeGfaSegmentsExpanded(ostream& gfa, bool includeSequence) const
{
    const AssemblyGraph& graph = *this;

    // Loop over BubbleChains. Each Chain of each Bubble generates a GFA segment.
    BGL_FORALL_EDGES(ce, graph, AssemblyGraph) {
        const BubbleChain& bubbleChain = graph[ce];

        // Loop over Bubbles of this chain.
        for(uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size();
            ++positionInBubbleChain) {
            const Bubble& bubble = bubbleChain[positionInBubbleChain];

            // Loop over chains of this bubble.
            for(uint64_t indexInBubble=0; indexInBubble<bubble.size(); indexInBubble++) {
                const Chain& chain = bubble[indexInBubble];
                const uint64_t offset = chainOffset(chain);

                // Record type.
                gfa << "S\t";

                // Name.
                gfa << chainStringId(ce, positionInBubbleChain, indexInBubble) << "\t";

                if(includeSequence) {
                    using shasta::Base;
                    const vector<Base>& sequence = chain.sequence;

                    // Sequence.
                    copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(gfa));
                    gfa << "\t";

                    // Sequence length in bases.
                    gfa << "LN:i:" << sequence.size() << "\n";

                } else {

                    // Sequence.
                    gfa << "*\t";

                    // Sequence length in bases.
                    gfa << "LN:i:" << offset << "\n";
                }
            }
        }
    }
}



void AssemblyGraph::writeGfaLinksExpanded(ostream& gfa) const
{
    const AssemblyGraph& graph = *this;

    // Write links between adjacent Chains of each BubbleChain.
    BGL_FORALL_EDGES(ce, graph, AssemblyGraph) {
        const BubbleChain& bubbleChain = graph[ce];

        // Loop over Bubbles of this chain.
        for(uint64_t positionInBubbleChain=1; positionInBubbleChain<bubbleChain.size();
            ++positionInBubbleChain) {
            const Bubble& bubble0 = bubbleChain[positionInBubbleChain - 1];
            const Bubble& bubble1 = bubbleChain[positionInBubbleChain];
            const uint64_t overlapLength = assembler.markerGraph.edgeSequence[bubble1.front().front()].size();

            for(uint64_t indexInBubble0=0; indexInBubble0<bubble0.size(); indexInBubble0++) {
                const string chain0StringId = chainStringId(ce, positionInBubbleChain-1, indexInBubble0);

                for(uint64_t indexInBubble1=0; indexInBubble1<bubble1.size(); indexInBubble1++) {
                   const string chain1StringId = chainStringId(ce, positionInBubbleChain, indexInBubble1);

                   gfa <<
                       "L\t" <<
                       chain0StringId << "\t+\t" <<
                       chain1StringId << "\t+\t" << overlapLength << "M\n";
                }
            }
        }
    }



    // Write links between Chains in different bubble chains.
    BGL_FORALL_VERTICES(cv, graph, AssemblyGraph) {
        const uint64_t overlapLength = assembler.markerGraph.edgeSequence[graph[cv].edgeId].size();

        BGL_FORALL_INEDGES(cv, ce0, graph, AssemblyGraph) {
            const BubbleChain& bubbleChain0 = graph[ce0];
            const Bubble& bubble0 = bubbleChain0.back();
            BGL_FORALL_OUTEDGES(cv, ce1, graph, AssemblyGraph) {
                const BubbleChain& bubbleChain1 = graph[ce1];
                const Bubble& bubble1 = bubbleChain1.front();

                for(uint64_t indexInBubble0=0; indexInBubble0<bubble0.size(); indexInBubble0++) {
                    const string chain0StringId = chainStringId(ce0, bubbleChain0.size()-1, indexInBubble0);

                    for(uint64_t indexInBubble1=0; indexInBubble1<bubble1.size(); indexInBubble1++) {
                       const string chain1StringId = chainStringId(ce1, 0, indexInBubble1);

                       gfa <<
                           "L\t" <<
                           chain0StringId << "\t+\t" <<
                           chain1StringId << "\t+\t" << overlapLength << "M\n";
                    }
                }
            }
        }
    }


}



void AssemblyGraph::writeGfaHeader(ostream& gfa)
{
    gfa << "H\tVN:Z:1.0\n";
}


// This version writes each chain as a segment, so it shows the
// details of the BubbleChains.
void AssemblyGraph::writeGfaExpanded(
    const string& fileNamePrefix,
    bool includeSequence) const
{
    ofstream gfa(fileNamePrefix + "-" + to_string(componentId) + "-Expanded.gfa");
    writeGfaExpanded(gfa, includeSequence);
}




void AssemblyGraph::writeFastaExpanded(const string& fileNamePrefix) const
{
    ofstream fasta(fileNamePrefix + "-" + to_string(componentId) + "-Expanded.fasta");
    writeFastaExpanded(fasta);
}



void AssemblyGraph::writeFastaExpanded(ostream& fasta) const
{
    const AssemblyGraph& cGraph = *this;


    // Loop over BubbleChains. Each Chain of each Bubble generates a GFA segment.
    BGL_FORALL_EDGES(ce, cGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = cGraph[ce];

        // Loop over Bubbles of this chain.
        for(uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size();
            ++positionInBubbleChain) {
            const Bubble& bubble = bubbleChain[positionInBubbleChain];

            // Loop over chains of this bubble.
            for(uint64_t indexInBubble=0; indexInBubble<bubble.size(); indexInBubble++) {
                const Chain& chain = bubble[indexInBubble];

                using shasta::Base;
                const vector<Base>& sequence = chain.sequence;

                fasta << ">" << chainStringId(ce, positionInBubbleChain, indexInBubble) <<
                    " " << sequence.size() << "\n";
                copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(fasta));
                fasta << "\n";



            }
        }
    }
}



void AssemblyGraph::writeSnapshot(uint64_t& snapshotNumber) const
{
    const string name = to_string(snapshotNumber++);
    write(name);
    writeGfaExpanded(name, false);
}



string AssemblyGraph::bubbleChainStringId(edge_descriptor ce) const
{
    const AssemblyGraph& cGraph = *this;
    const AssemblyGraphEdge& edge = cGraph[ce];
    return to_string(componentId) + "-" + to_string(edge.id);
}



string AssemblyGraph::bubbleStringId(
    edge_descriptor ce,
    uint64_t positionInBubbleChain) const
{
    const AssemblyGraph& cGraph = *this;
    const AssemblyGraphEdge& edge = cGraph[ce];

    return
        to_string(componentId) + "-" +
        to_string(edge.id) + "-" +
        to_string(positionInBubbleChain);
}



string AssemblyGraph::chainStringId(
    edge_descriptor ce,
    uint64_t positionInBubbleChain,
    uint64_t indexInBubble) const
{
    // Locate the bubble chain.
    const AssemblyGraph& cGraph = *this;
    const AssemblyGraphEdge& edge = cGraph[ce];
    const BubbleChain& bubbleChain = edge;

    // Find the ploidy.
    const Bubble& bubble = bubbleChain[positionInBubbleChain];
    const uint64_t ploidy = bubble.size();

    return
        to_string(componentId) + "-" +
        to_string(edge.id) + "-" +
        to_string(positionInBubbleChain) + "-" +
        to_string(indexInBubble) + "-P" +
        to_string(ploidy);
}



uint64_t AssemblyGraph::chainOffset(const Chain& chain) const
{
    const uint64_t length = chain.size();
    SHASTA_ASSERT(length >= 2);

    uint64_t offset = 0;
    for(uint64_t i=1; i<length; i++) {
        const MarkerGraphEdgeId edgeId0 = chain[i-1];
        const MarkerGraphEdgeId edgeId1 = chain[i];

        const uint64_t offsetThisPair = assembler.estimateBaseOffsetUnsafe(edgeId0, edgeId1);

        if(offsetThisPair != invalid<uint64_t>) {
            offset += offsetThisPair;
        }
    }
    return offset;
}



void AssemblyGraph::bubbleOffset(
    const Bubble& bubble,
    uint64_t& averageOffset,
    uint64_t& minOffset,
    uint64_t& maxOffset
    ) const
{
    averageOffset = 0;
    minOffset = std::numeric_limits<uint64_t>::max();
    maxOffset = 0;

    for(const Chain& chain: bubble) {
        const uint64_t offset = chainOffset(chain);

        averageOffset += offset;
        minOffset = min(minOffset, offset);
        maxOffset = max(maxOffset, offset);
    }
    averageOffset /= bubble.size();
}



bool AssemblyGraph::bubbleOffsetNoException(
    const Bubble& bubble,
    uint64_t& averageOffset,
    uint64_t& minOffset,
    uint64_t& maxOffset
    ) const
{
    averageOffset = 0;
    minOffset = std::numeric_limits<uint64_t>::max();
    maxOffset = 0;

    for(const Chain& chain: bubble) {
        const uint64_t offset = chainOffset(chain);
        if(offset == invalid<uint64_t>) {
            return false;
        }

        averageOffset += offset;
        minOffset = min(minOffset, offset);
        maxOffset = max(maxOffset, offset);
    }
    averageOffset /= bubble.size();
    return true;
}



void AssemblyGraph::bubbleChainOffset(
    const BubbleChain& bubbleChain,
    uint64_t& averageOffset,
    uint64_t& minOffset,
    uint64_t& maxOffset
    ) const
{
    averageOffset = 0;
    minOffset = 0;
    maxOffset = 0;

    for(const Bubble& bubble: bubbleChain) {
        uint64_t bubbleAverageOffset;
        uint64_t bubbleMinOffset;
        uint64_t bubbleMaxOffset;
        bubbleOffset(bubble, bubbleAverageOffset, bubbleMinOffset, bubbleMaxOffset);

        averageOffset += bubbleAverageOffset;
        minOffset += bubbleMinOffset;
        maxOffset += bubbleMaxOffset;
    }
}



AssemblyGraph::Superbubbles::Superbubbles(
    AssemblyGraph& cGraph,
    uint64_t maxOffset1     // Used to define superbubbles
    ) :
    cGraph(cGraph)
{
    cGraph.numberVertices();
    const uint64_t vertexCount = num_vertices(cGraph);

    vector<uint64_t> rank(vertexCount);
    vector<uint64_t> parent(vertexCount);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);

    // Compute connected components, using only edges with average offset up to maxOffset1.
    for(uint64_t i=0; i<vertexCount; i++) {
        disjointSets.make_set(i);
    }
    BGL_FORALL_EDGES(ce, cGraph, AssemblyGraph) {
        uint64_t averageOffset;
        uint64_t minOffset;
        uint64_t maxOffset;
        cGraph.bubbleChainOffset(cGraph[ce], averageOffset, minOffset, maxOffset);
        if(averageOffset <= maxOffset1) {
            const vertex_descriptor cv0 = source(ce, cGraph);
            const vertex_descriptor cv1 = target(ce, cGraph);
            disjointSets.union_set(cGraph[cv0].index, cGraph[cv1].index);
        }
    }

    // Gather the vertices in each connected component.
    vector< vector<vertex_descriptor> > components(vertexCount);
    BGL_FORALL_VERTICES(cv, cGraph, AssemblyGraph) {
        const uint64_t componentId = disjointSets.find_set(cGraph[cv].index);
        components[componentId].push_back(cv);
    }

    // The superbubbles are the components with size at least 2.
    for(uint64_t componentId=0; componentId<components.size(); componentId++) {
        const vector<vertex_descriptor> component = components[componentId];
        if(components[componentId].size() > 1) {
            superbubbles.emplace_back(Superbubble(component));
        }
    }

    // Store superbubble ids in the vertices.
    BGL_FORALL_VERTICES(cv, cGraph, AssemblyGraph) {
        cGraph[cv].superbubbleId = invalid<uint64_t>;
    }
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        const vector<vertex_descriptor>& superbubble = getSuperbubble(superbubbleId);
        for(const vertex_descriptor cv: superbubble) {
            cGraph[cv].superbubbleId = superbubbleId;
        }
    }



    // Find entrances and exists of each superbubble.
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        Superbubble& superbubble = getSuperbubble(superbubbleId);

        // Find entrances. These are superbubble vertices with in-edges
        // from outside the superbubble.
        for(const vertex_descriptor cv0: superbubble) {
            BGL_FORALL_INEDGES(cv0, ce, cGraph, AssemblyGraph) {
                const vertex_descriptor cv1 = source(ce, cGraph);
                if(not isInSuperbubble(superbubbleId, cv1)) {
                    superbubble.entrances.push_back(cv0);
                    break;
                }
            }
        }

        // Find exits. These are superbubble vertices with out-edges
        // to outside the superbubble.
        vector<vertex_descriptor> exits;
        for(const vertex_descriptor cv0: superbubble) {
            BGL_FORALL_OUTEDGES(cv0, ce, cGraph, AssemblyGraph) {
                const vertex_descriptor cv1 = target(ce, cGraph);
                if(not isInSuperbubble(superbubbleId, cv1)) {
                    superbubble.exits.push_back(cv0);
                    break;
                }
            }
        }
     }

}



// This uses dominator trees.
// It only finds superbubbles with one entrance and one exit.
AssemblyGraph::Superbubbles::Superbubbles(
    AssemblyGraph& cGraph) :
    cGraph(cGraph)
{
    const bool debug = false;

    // Map vertices to integers.
    std::map<vertex_descriptor, uint64_t> indexMap;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, cGraph, AssemblyGraph) {
        indexMap.insert({v, vertexIndex++});
    }
    auto associativeIndexMap = boost::make_assoc_property_map(indexMap);
    const uint64_t vertexCount = vertexIndex;

    // Vectors used below to compute the dominator tree.
    vector<uint64_t> dfNum(vertexCount);
    vector<vertex_descriptor> parent(vertexCount);
    vector<vertex_descriptor> verticesByDFNum(vertexCount);

    // Tree pairs found on forward and backward dominator tree.
    vector< pair<vertex_descriptor, vertex_descriptor> > forwardPairs;
    vector< pair<vertex_descriptor, vertex_descriptor> > backwardPairs;



    // Compute dominator trees using as entrance each of the
    // vertices with zero in-degree.
    BGL_FORALL_VERTICES(entrance, cGraph, AssemblyGraph) {
        if(in_degree(entrance, cGraph) != 0) {
            continue;
        }

        // Compute the dominator tree.
        fill(dfNum.begin(), dfNum.end(), invalid<uint64_t>);
        fill(parent.begin(), parent.end(), null_vertex());
        fill(verticesByDFNum.begin(), verticesByDFNum.end(), null_vertex());
        std::map<vertex_descriptor, vertex_descriptor> predecessorMap;

        boost::lengauer_tarjan_dominator_tree(
            cGraph,
            entrance,
            boost::make_assoc_property_map(indexMap),
            boost::make_iterator_property_map(dfNum.begin(), associativeIndexMap),
            boost::make_iterator_property_map(parent.begin(), associativeIndexMap),
            verticesByDFNum,
            boost::make_assoc_property_map(predecessorMap));

        if(debug) {
            cout << "Forward dominator tree with entrance at " << cGraph[entrance].edgeId << endl;
        }
        for(const auto& p: predecessorMap) {
            const vertex_descriptor cv0 = p.second;
            const vertex_descriptor cv1 = p.first;
            forwardPairs.push_back({cv0, cv1});
            if(debug) {
                cout << "F " << cGraph[cv0].edgeId << "->" << cGraph[cv1].edgeId << endl;
            }
        }
    }



    // Compute dominator trees on the reverse graph using as entrance each of the
    // vertices with zero in-degree on the reverse graph
    // (that is, zero out-degree on the AssemblyGraph).
    using ReverseAssemblyGraph = boost::reverse_graph<AssemblyGraph>;
    ReverseAssemblyGraph reverseGraph(cGraph);
    BGL_FORALL_VERTICES(entrance, reverseGraph, ReverseAssemblyGraph) {
        if(in_degree(entrance, reverseGraph) != 0) {
            continue;
        }

        // Compute the dominator tree.
        fill(dfNum.begin(), dfNum.end(), invalid<uint64_t>);
        fill(parent.begin(), parent.end(), null_vertex());
        fill(verticesByDFNum.begin(), verticesByDFNum.end(), null_vertex());
        std::map<vertex_descriptor, vertex_descriptor> predecessorMap;

        boost::lengauer_tarjan_dominator_tree(
            reverseGraph,
            entrance,
            boost::make_assoc_property_map(indexMap),
            boost::make_iterator_property_map(dfNum.begin(), associativeIndexMap),
            boost::make_iterator_property_map(parent.begin(), associativeIndexMap),
            verticesByDFNum,
            boost::make_assoc_property_map(predecessorMap));

        if(debug) {
            cout << "Backward dominator tree with exit at " << cGraph[entrance].edgeId << endl;
        }
        for(const auto& p: predecessorMap) {
            const vertex_descriptor cv0 = p.first;
            const vertex_descriptor cv1 = p.second;
            backwardPairs.push_back({cv0, cv1});
            if(debug) {
                cout << "B " << cGraph[cv0].edgeId << "->" << cGraph[cv1].edgeId << endl;
            }
        }
    }

    // Compute strongly connected components.
    std::map<vertex_descriptor, uint64_t> componentMap;
    boost::strong_components(
        cGraph,
        boost::make_assoc_property_map(componentMap),
        boost::vertex_index_map(boost::make_assoc_property_map(indexMap)));

    // Gather the vertices in each strong component.
    vector< vector<vertex_descriptor> > strongComponents(vertexCount);
    for(const auto& p: componentMap) {
        const vertex_descriptor v = p.first;
        const uint64_t componentId = p.second;
        SHASTA_ASSERT(componentId < vertexCount);
        strongComponents[componentId].push_back(v);
    }



    // The pairs that appear both in forwardPairs and backwardPairs define our superbubbles
    deduplicate(forwardPairs);
    deduplicate(backwardPairs);
    vector< pair<vertex_descriptor, vertex_descriptor> > bidirectionalPairs;
    std::set_intersection(
        forwardPairs.begin(), forwardPairs.end(),
        backwardPairs.begin(), backwardPairs.end(),
        back_inserter(bidirectionalPairs)
        );

    if(debug) {
        cout << "Bidirectional pairs:" << endl;
        for(const auto& p: bidirectionalPairs) {
            const vertex_descriptor cv0 = p.first;
            const vertex_descriptor cv1 = p.second;
            cout << cGraph[cv0].edgeId << "->" << cGraph[cv1].edgeId << endl;
        }
    }

    // Each bidirectional pair generates a superbubble if
    // the out-degree of the entrance and
    // the in-degree of the exit are greater than 1,
    // unless the entrance or exit or any of the
    // superbubble vertices are in a non-trivial strong component..
    for(const auto& p: bidirectionalPairs) {
        const vertex_descriptor cv0 = p.first;
        const vertex_descriptor cv1 = p.second;
        if(out_degree(cv0, cGraph) <= 1) {
            continue;
        }
        if(in_degree(cv1, cGraph) <= 1) {
            continue;
        }
        if(strongComponents[componentMap[cv0]].size() > 1) {
            // The entrance is in a non-trivial strong component.
            continue;
        }
        if(strongComponents[componentMap[cv1]].size() > 1) {
            // The exit is in a non-trivial strong component.
            continue;
        }
        superbubbles.resize(superbubbles.size() + 1);
        Superbubble& superbubble = superbubbles.back();
        superbubble.entrances.push_back(cv0);
        superbubble.exits.push_back(cv1);
        superbubble.fillInFromEntranceAndExit(cGraph);

        if(debug) {
            cout << "Tentative superbubble with entrance " << cGraph[cv0].edgeId <<
                " exit " << cGraph[cv1].edgeId << " and " << superbubble.size() <<
                " vertices total." << endl;
        }

        // If any vertices in the superbubble are in a non-trivial
        // strong component, remove it.
        for(const vertex_descriptor cv: superbubble) {
            if(strongComponents[componentMap[cv]].size() > 1) {
                superbubbles.pop_back();
                if(debug) {
                    cout << "This superbubble will not be stored because some vertices are in a non-trivial strong component." << endl;
                }
                break;
            }
        }
    }

    if(debug) {
        cout << "Superbubble entrance/exit pairs:" << endl;
        for(const Superbubble& superbubble: superbubbles) {
            const vertex_descriptor cv0 = superbubble.entrances.front();
            const vertex_descriptor cv1 = superbubble.exits.front();;
            cout << cGraph[cv0].edgeId << "->" << cGraph[cv1].edgeId << endl;
        }
    }
}



// Fill in the superbubble given a single entrance and exit.
void AssemblyGraph::Superbubble::fillInFromEntranceAndExit(const AssemblyGraph& cGraph)
{
    SHASTA_ASSERT(empty());
    SHASTA_ASSERT(entrances.size() == 1);
    SHASTA_ASSERT(exits.size() == 1);

    const vertex_descriptor entrance = entrances.front();
    const vertex_descriptor exit = exits.front();

    // Do a BFS starting at the entrance and stopping at the exit.
    std::set<vertex_descriptor> internalVertices;
    std::queue<vertex_descriptor> q;
    q.push(entrance);
    while(not q.empty()) {
        const vertex_descriptor cv0 = q.front();
        q.pop();
        BGL_FORALL_OUTEDGES(cv0, e, cGraph, AssemblyGraph) {
            const vertex_descriptor cv1 = target(e, cGraph);
            if(cv1 != exit) {
                if(not internalVertices.contains(cv1)) {
                    internalVertices.insert(cv1);
                    q.push(cv1);
                }
            }
        }
    }

    push_back(entrance);
    copy(internalVertices.begin(), internalVertices.end(), back_inserter(*this));
    push_back(exit);

}



AssemblyGraph::Superbubbles::~Superbubbles()
{
    cGraph.clearVertexNumbering();
}



// Remove short superbubbles with one entry and one exit.
bool AssemblyGraph::removeShortSuperbubbles(
    bool debug,
    uint64_t maxOffset1,    // Used to define superbubbles
    uint64_t maxOffset2)    // Compared against the offset between entry and exit
{
    AssemblyGraph& cGraph = *this;
    bool changesWereMade = false;

    // Find the superbubbles.
    Superbubbles superbubbles(cGraph, maxOffset1);

    // Loop over the superbubbles.
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        Superbubble& superbubble = superbubbles.getSuperbubble(superbubbleId);
        SHASTA_ASSERT(superbubble.size() > 1);

        if(debug) {
            cout << "Found a superbubble with " << superbubble.size() << " vertices:";
            for(const vertex_descriptor v: superbubble) {
                cout << " " << cGraph[v].edgeId;
            }
            cout << endl;
        }

        // Skip it if it has more than one entrance or exit.
        if(not(superbubble.entrances.size()==1 and superbubble.exits.size()==1)) {
            if(debug) {
                cout << "This superbubble will not be removed because it has " <<
                    superbubble.entrances.size() << " entrances and " <<
                    superbubble.exits.size() << " exits." << endl;
            }
            continue;
        }

        const vertex_descriptor entrance = superbubble.entrances.front();
        const vertex_descriptor exit = superbubble.exits.front();
        if(entrance == exit) {
            if(debug) {
                cout << "This superbubble will not be removed because it the entrance vertex"
                    " is the same as the exit vertex." << endl;
            }
            continue;
        }

        // Check the base offset between the entrance and the exit.
        MarkerGraphEdgePairInfo info;
        SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(cGraph[entrance].edgeId, cGraph[exit].edgeId, info));
        if(info.common == 0) {
            if(debug) {
                cout << "This superbubble will not be removed because "
                    "there are no common oriented reads between the entrance and the exit." << endl;
            }
            continue;
        }
        if(info.offsetInBases > int64_t(maxOffset2)) {
            if(debug) {
                cout << "This superbubble will not be removed because offsetInBases is " <<
                    info.offsetInBases << endl;
            }
            continue;
        }

#if 1
        // If a trivial superbubble, skip it.
        // Trivial means:
        // - Has two vertices of which one is the entrance and one is the exit.
        // - There is only one edge between the two.
        if(superbubble.size() == 2) {
            uint64_t edgeCount = 0;
            BGL_FORALL_OUTEDGES(entrance, e, cGraph, AssemblyGraph) {
                if(target(e, cGraph) == exit) {
                    ++edgeCount;
                }
            }
            if(edgeCount == 1) {
                if(debug) {
                    cout << "This superbubble will not be removed because it is trivial." << endl;
                }
                continue;
            }
        }
#endif
        if(debug) {
            cout << "This superbubble will be removed." << endl;
        }

        // Remove all vertices and edges internal to the superbubble.
        for(const vertex_descriptor cv: superbubble) {
            if(cv!=entrance and cv!=exit) {
                boost::clear_vertex(cv, cGraph);
                cGraph.removeVertex(cv);
            }
        }
        // We must also remove edges between the entrance and the exit.
        vector<edge_descriptor> entranceToExitEdges;
        BGL_FORALL_OUTEDGES(entrance, ce, cGraph, AssemblyGraph) {
            if(target(ce, cGraph) == exit) {
                entranceToExitEdges.push_back(ce);
            }
        }
        for(const edge_descriptor ce: entranceToExitEdges) {
            boost::remove_edge(ce, cGraph);
        }
        vector<edge_descriptor> exitToEntranceEdges;
        BGL_FORALL_OUTEDGES(exit, ce, cGraph, AssemblyGraph) {
            if(target(ce, cGraph) == entrance) {
                exitToEntranceEdges.push_back(ce);
            }
        }
        for(const edge_descriptor ce: exitToEntranceEdges) {
            boost::remove_edge(ce, cGraph);
        }

        // Generate an edge between the entrance and the exit.
        // This will be a BubbleChain consisting of a single haploid Bubble.
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(entrance, exit, cGraph);
        AssemblyGraphEdge& newEdge = cGraph[eNew];
        newEdge.id = nextEdgeId++;
        BubbleChain& bubbleChain = newEdge;
        bubbleChain.resize(1);
        Bubble& bubble = bubbleChain.front();
        bubble.resize(1);
        Chain& chain = bubble.front();
        chain.push_back(cGraph[entrance].edgeId);
        chain.push_back(cGraph[exit].edgeId);

        changesWereMade = true;
    }

    return changesWereMade;
}



// Cleanup/simplify superbubbles that are likely to be caused by errors,
// completely or in part.
void AssemblyGraph::cleanupSuperbubbles(
    bool debug,
    uint64_t maxOffset1,    // Used to define superbubbles
    uint64_t maxOffset2,    // Compared against the offset between entry and exit
    uint64_t chainTerminalCommonThreshold

)
{
    AssemblyGraph& cGraph = *this;

    if(debug) {
        cout << "cleanupSuperbubbles begins." << endl;
    }

    // Find the superbubbles.
    Superbubbles superbubbles(cGraph, maxOffset1);

    // The bubbles constructed in this way are guaranteed to not overlap,
    // so we don't have to worry about overlapping bubbles.
    std::set<vertex_descriptor> previousSuperbubblesVertices;

    // Loop over the superbubbles.
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        cleanupSuperbubble(debug, superbubbles, superbubbleId,
            maxOffset2, chainTerminalCommonThreshold, previousSuperbubblesVertices);
    }
    if(debug) {
        cout << "cleanupSuperbubbles ends." << endl;
    }
}



// This version of superbubble cleanup uses dominator trees to define superbubbles,
// instead of computing connected components using edges of length uo tp maxOffset1.
void AssemblyGraph::cleanupSuperbubbles(
    bool debug,
    uint64_t maxOffset2,     // Compared against the offset between entry and exit
    uint64_t chainTerminalCommonThreshold
    )
{
    performanceLog << timestamp << "AssemblyGraph::cleanupSuperbubbles begins." << endl;
    AssemblyGraph& cGraph = *this;

    if(debug) {
        cout << "cleanupSuperbubbles begins." << endl;
    }

    // Find the superbubbles using dominator trees.
    Superbubbles superbubbles(cGraph);

    // The superbubbles found in this way can have overlaps.
    // To deal with this, we process superbubbles in order of increasing size
    // and keep track of the vertices.
    // If a bubble contains a previously encountered vertex, don't process it.
    // Note cleanupSuperbubble does not create any new vertices,
    // so keeping track of the vertex descriptors that were removed is save.
    std::set<vertex_descriptor> previousSuperbubblesVertices;

    // Sort the superbubbles in order of increasing size.
    vector< pair<uint64_t, uint64_t> > superbubbleTable;    // (superbubbleId, size)
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        const Superbubble& superbubble = superbubbles.getSuperbubble(superbubbleId);
        superbubbleTable.push_back({superbubbleId, superbubble.size()});
    }
    sort(superbubbleTable.begin(), superbubbleTable.end(),
        OrderPairsBySecondOnly<uint64_t, uint64_t>());

    // Loop over the superbubbles in order of increasing size.
    for(const auto& p: superbubbleTable) {
        const uint64_t superbubbleId = p.first;
        cleanupSuperbubble(debug, superbubbles, superbubbleId, maxOffset2,
            chainTerminalCommonThreshold, previousSuperbubblesVertices);
    }
    if(debug) {
        cout << "cleanupSuperbubbles ends." << endl;
    }
    performanceLog << timestamp << "AssemblyGraph::cleanupSuperbubbles ends." << endl;

}



// Cleanup/simplify a superbubble that is likely to be caused by errors,
// completely or in part.
// This handles superbubbles caused by two marker graph bubbles with
// no primary edges in between.
void AssemblyGraph::cleanupSuperbubble(
    bool debug,
    const Superbubbles& superbubbles,
    uint64_t superbubbleId,
    uint64_t maxOffset2,        // Compared against the offset between entry and exit
    uint64_t chainTerminalCommonThreshold,
    std::set<vertex_descriptor>& previousSuperbubblesVertices)
{
    AssemblyGraph& cGraph = *this;
    const Superbubble& superbubble = superbubbles.getSuperbubble(superbubbleId);

#if 0
    debug = (superbubble.entrances.size() == 1 and
        (cGraph[superbubble.entrances.front()].edgeId == 16093908 or
        cGraph[superbubble.entrances.front()].edgeId == 9555933));
#endif

    if(debug) {
        cout << "Working on a superbubble with " << superbubble.size() << " vertices:";
        for(const vertex_descriptor v: superbubble) {
            cout << " " << cGraph[v].edgeId;
        }
        cout << endl;
    }

    // See if it overlaps any vertices of previous superbubbles.
    bool overlaps = false;
    for(const vertex_descriptor v: superbubble) {
        if(previousSuperbubblesVertices.contains(v)) {
            if(debug) {
                cout << "This superbubble ignored because it contains vertex " << cGraph[v].edgeId <<
                    " which is in a previously processed superbubble." << endl;
            }
            overlaps = true;
            break;
        }
    }
    for(const vertex_descriptor v: superbubble) {
        previousSuperbubblesVertices.insert(v);
    }
    if(overlaps) {
        return;
    }

    // Skip it if it has more than one entrance or exit.
    if(not(superbubble.entrances.size()==1 and superbubble.exits.size()==1)) {
        if(debug) {
            cout << "This superbubble will be skipped because it has " <<
                superbubble.entrances.size() << " entrances and " <<
                superbubble.exits.size() << " exits." << endl;
        }
        return;
    }

    const vertex_descriptor entrance = superbubble.entrances.front();
    const vertex_descriptor exit = superbubble.exits.front();
    if(debug) {
        cout << "Entrance " << cGraph[entrance].edgeId << endl;
        cout << "Exit " << cGraph[exit].edgeId << endl;
    }

    if(entrance == exit) {
        if(debug) {
            cout << "This superbubble will be skipped because the entrance vertex"
                " is the same as the exit vertex." << endl;
        }
        return;
    }



    // Check the base offset between the entrance and the exit.
    MarkerGraphEdgePairInfo info;
    SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(cGraph[entrance].edgeId, cGraph[exit].edgeId, info));
    if(info.common == 0) {
        if(debug) {
            cout << "This superbubble will be skipped because "
                "there are no common oriented reads between the entrance and the exit." << endl;
        }
        return;
    }
    if(info.offsetInBases > int64_t(maxOffset2)) {
        if(debug) {
            cout << "This superbubble will be skipped because offsetInBases is " <<
                info.offsetInBases << endl;
        }
        return;
    }

    // If a trivial superbubble, skip it.
    // Trivial means:
    // - Has two vertices of which one is the entrance and one is the exit.
    // - There is only one edge between the two.
    if(superbubble.size() == 2) {
        uint64_t edgeCount = 0;
        BGL_FORALL_OUTEDGES(entrance, e, cGraph, AssemblyGraph) {
            if(target(e, cGraph) == exit) {
                ++edgeCount;
            }
        }
        if(edgeCount == 1) {
            if(debug) {
                cout << "This superbubble be skipped because it is trivial." << endl;
            }
            return;
        }
    }

    // Find the out-edges of the entrance that go inside the superbubble.
    vector<edge_descriptor> entranceOutEdges;
    BGL_FORALL_OUTEDGES(entrance, ce, cGraph, AssemblyGraph) {
        const vertex_descriptor cv = target(ce, cGraph);
        if(superbubbles.isInSuperbubble(superbubbleId, cv)) {
            entranceOutEdges.push_back(ce);
        }
    }
    sort(entranceOutEdges.begin(), entranceOutEdges.end());

    // Find the in-edges of the exit that come from inside the superbubble.
    vector<edge_descriptor> exitInEdges;
    BGL_FORALL_INEDGES(exit, ce, cGraph, AssemblyGraph) {
        const vertex_descriptor cv = source(ce, cGraph);
        if(superbubbles.isInSuperbubble(superbubbleId, cv)) {
            exitInEdges.push_back(ce);
        }
    }
    sort(exitInEdges.begin(), exitInEdges.end());

    if(debug) {
        cout << "Entrance out-edges to inside the superbubble:";
        for(const edge_descriptor ce: entranceOutEdges) {
            cout << " " << bubbleChainStringId(ce);
        }
        cout << endl;
        cout << "Exit in-edges from inside the superbubble:";
        for(const edge_descriptor ce: exitInEdges) {
            cout << " " << bubbleChainStringId(ce);
        }
        cout << endl;
    }

    // If there are common edges between the entranceOutEdges and exitInEdges,
    // skip this superbubble.
    {
        vector<edge_descriptor> commonEdges;
        std::set_intersection(
            entranceOutEdges.begin(), entranceOutEdges.end(),
            exitInEdges.begin(), exitInEdges.end(),
            back_inserter(commonEdges));

        if(not commonEdges.empty()) {
            if(debug) {
                cout << "This superbubble will be skipped because there are " <<
                    commonEdges.size() << " common edges between the out-edges of the entrance "
                    "and the in-edges of the exit." << endl;
            }
            return;
        }
    }


    // We will consider replacing this superbubble with either its "entrance bubble"
    // or its "exit bubble":
    // - The "entrance bubble" is obtained by removing all edges
    //   except for the out-edges of the entrance, and joining them directly with the exit.
    // - The "exit bubble" is obtained by removing all edges
    //   except for the in-edges of the exit, and joining the entry directly with them.



    // If there are exactly two entranceOutEdges, construct the entrance bubble.
    // This can only be done if the two entranceOutEdges consist of simple chains.
    Bubble entranceBubble;
    if(entranceOutEdges.size() == 2) {

        // See if the two entranceOutEdges consist of simple chains.
        bool canDo = true;
        for(const edge_descriptor ce: entranceOutEdges) {
            if(not cGraph[ce].isSimpleChain()) {
                canDo = false;
                break;
            }
        }

        // Only continue creating the entranceBubble if both entranceOutEdges
        // consist of single chains.
        if(canDo) {

            // Construct the two chains of the entranceBubble and assemble their sequence.
            entranceBubble.resize(2);
            ofstream noCsv;
            for(uint64_t i=0; i<2; i++) {
                const edge_descriptor entranceOutEdge = entranceOutEdges[i];
                Chain& chain = entranceBubble[i];
                chain = cGraph[entranceOutEdge].getOnlyChain();
                chain.push_back(cGraph[exit].edgeId);
                assembleChain(chain, chainTerminalCommonThreshold);
            }

            if(debug) {
                cout << "Entrance bubble:" << endl;
                for(uint64_t i=0; i<2; i++) {
                    const Chain& chain = entranceBubble[i];
                    cout << "Entrance bubble chain " << i << ":";
                    for (const MarkerGraphEdgeId edgeId: chain) {
                        cout << " " << edgeId;
                    }
                    cout << endl;
                }
                for(uint64_t i=0; i<2; i++) {
                    const Chain& chain = entranceBubble[i];
                    cout << ">Entrance-" << i << " " << chain.sequence.size() << "\n";
                    copy(chain.sequence.begin(), chain.sequence.end(), ostream_iterator<shasta::Base>(cout));
                    cout << "\n";
                }
            }

            // If the sequences differ just by a copy number of short periodicity,
            // the entrance bubble is probably causes by errors and so we don't wat to use it.
            const uint64_t period = isCopyNumberDifference(entranceBubble[0].sequence, entranceBubble[1].sequence, 4);
            if(debug) {
                cout << "Period " << period << "\n";
            }
            if(period != 0) {
                entranceBubble.clear();
            }
        }
    }



    // If there are exactly two exitEdges, construct the exit bubble.
    // This can only be done if the two exitInEdges consist of simple chains.
    Bubble exitBubble;
    if(exitInEdges.size() == 2) {

        // See if the two exitInEdges consist of simple chains.
        bool canDo = true;
        for(const edge_descriptor ce: exitInEdges) {
            if(not cGraph[ce].isSimpleChain()) {
                canDo = false;
                break;
            }
        }

        // Only continue creating the exitBubble if both exitInEdges
        // consist of single chains.
        if(canDo) {

            // Construct the two chains of the exitBubble and assemble their sequence.
            exitBubble.resize(2);
            ofstream noCsv;
            for(uint64_t i=0; i<2; i++) {
                const edge_descriptor exitInEdge = exitInEdges[i];
                Chain& chain = exitBubble[i];
                chain.push_back(cGraph[entrance].edgeId);
                const Chain& exitChain = cGraph[exitInEdge].getOnlyChain();
                copy(exitChain.begin(), exitChain.end(), back_inserter(chain));
                assembleChain(chain, chainTerminalCommonThreshold);
            }

            if(debug) {
                cout << "Exit bubble:" << endl;
                for(uint64_t i=0; i<2; i++) {
                    const Chain& chain = exitBubble[i];
                    cout << "Exit bubble chain " << i << ":";
                    for (const MarkerGraphEdgeId edgeId: chain) {
                        cout << " " << edgeId;
                    }
                    cout << endl;
                }
                for(uint64_t i=0; i<2; i++) {
                    const Chain& chain = exitBubble[i];
                    cout << ">Exit-" << i << " " << chain.sequence.size() << "\n";
                    copy(chain.sequence.begin(), chain.sequence.end(), ostream_iterator<shasta::Base>(cout));
                    cout << "\n";
                }
            }

            // If the sequences differ just by a copy number of short periodicity,
            // the exit bubble is probably causes by errors and so we don't wat to use it.
            const uint64_t period = isCopyNumberDifference(exitBubble[0].sequence, exitBubble[1].sequence, 4);
            if(debug) {
                cout << "Period " << period << "\n";
            }
            if(period != 0) {
                exitBubble.clear();
            }
        }
    }


    // Handle the case where both the entrance and the exit bubble look usable.
    if(entranceBubble.size() == 2 and exitBubble.size() == 2) {

        // If the entrance and exit bubbles have the same assembled sequences, we can just keep one of them.
        const auto& entrance0 = entranceBubble[0].sequence;
        const auto& entrance1 = entranceBubble[1].sequence;
        const auto& exit0 = exitBubble[0].sequence;
        const auto& exit1 = exitBubble[1].sequence;
        if(
            (entrance0 == exit0 and entrance1 == exit1)
            or
            (entrance0 == exit1 and entrance1 == exit0)) {
            if(debug) {
                cout << "The entrance and exit bubbles are equivalent." << endl;
                cout << "Keeping only the entrance bubble." << endl;
            }
            exitBubble.clear();
        } else {

            // In other cases it is difficult to pick which bubble is best to keep,
            // so we remove both of them.
            // This is no worse than letting removeShortBubbles remove it.
            // The sequence assembly process will still pick the best sequence
            // for each haplotype, but these bubbles are excluded from the
            // phasing/detangling process.
            entranceBubble.clear();
            exitBubble.clear();

            if(debug) {
                cout << "Both the entrance and the exit bubble are usable but both will be removed." << endl;
            }

        }
    }



    // Figure out which ones of the entrance/exit bubbles is usable.
    SHASTA_ASSERT(entranceBubble.size() == 0 or entranceBubble.size() == 2);
    SHASTA_ASSERT(exitBubble.size() == 0 or exitBubble.size() == 2);
    const bool entranceBubbleIsGood = (entranceBubble.size() == 2);
    const bool exitBubbleIsGood = (exitBubble.size() == 2);


    if(entranceBubbleIsGood) {
        if(exitBubbleIsGood) {
            if(debug) {
                cout << "Both the entrance bubble and the exit bubble are good." << endl;
            }
            SHASTA_ASSERT(0);
        } else {
            if(debug) {
                cout << "Only the entrance bubble is good." << endl;
            }
        }
    } else {
        if(exitBubbleIsGood) {
            if(debug) {
                cout << "Only the exit bubble is good." << endl;
            }
        } else {
            if(debug) {
                cout << "Neither the entrance bubble nor the exit bubble are good." << endl;
            }
        }
    }


    // Remove all vertices and edges internal to the superbubble.
    for(const vertex_descriptor cv: superbubble) {
        if(cv != entrance and cv != exit) {
            clear_vertex(cv, cGraph);
            remove_vertex(cv, cGraph);
        }
    }

    // Create the new edge and bubble chain between the entrance and the exit that will replace
    // the superbubble.
    edge_descriptor ce;
    tie(ce, ignore) = add_edge(entrance, exit, cGraph);
    AssemblyGraphEdge& edge = cGraph[ce];
    edge.id = nextEdgeId++;
    BubbleChain& bubbleChain = edge;
    SHASTA_ASSERT(not (entranceBubbleIsGood and exitBubbleIsGood));
    if(entranceBubbleIsGood or exitBubbleIsGood) {
        const Bubble& newBubble = entranceBubbleIsGood ? entranceBubble : exitBubble;
        SHASTA_ASSERT(newBubble.size() == 2);
        bubbleChain.push_back(newBubble);
    } else {
        Chain newChain;
        newChain.push_back(cGraph[entrance].edgeId);
        newChain.push_back(cGraph[exit].edgeId);
        Bubble newBubble;
        newBubble.push_back(newChain);
        bubbleChain.push_back(newBubble);
    }

}



#if 0
bool AssemblyGraph::detangleVerticesStrict(bool debug)
{
    if(debug) {
        cout << "Detangling vertices." << endl;
    }
    AssemblyGraph& cGraph = *this;

    vector<vertex_descriptor> allVertices;
    BGL_FORALL_VERTICES(cv, cGraph, AssemblyGraph) {
        allVertices.push_back(cv);
    }

    uint64_t detangledCount = 0;
    for(const vertex_descriptor cv: allVertices) {
        if(detangleVertexStrict(cv, debug)) {
            ++detangledCount;
        }
    }

    if(debug) {
        cout << "Detangled " << detangledCount << " vertices." << endl;

    }

    return detangledCount > 0;
}
#endif



bool AssemblyGraph::detangleVertices(
    bool debug,
    uint64_t detangleToleranceLow,
    uint64_t detangleToleranceHigh,
    bool useBayesianModel,
    double epsilon,
    double minLogP)
{
    if(debug) {
        cout << "Detangling vertices." << endl;
    }
    AssemblyGraph& cGraph = *this;

    vector<vertex_descriptor> allVertices;
    BGL_FORALL_VERTICES(cv, cGraph, AssemblyGraph) {
        allVertices.push_back(cv);
    }

    uint64_t detangledCount = 0;
    for(const vertex_descriptor cv: allVertices) {
        if(detangleVertex(cv, debug, detangleToleranceLow, detangleToleranceHigh,
            useBayesianModel, epsilon, minLogP)) {
            ++detangledCount;
        }
    }

    if(debug) {
        cout << "Detangled " << detangledCount << " vertices." << endl;
    }

    return detangledCount > 0;
}



bool AssemblyGraph::detangleVerticesGeneral(
    bool debug,
    uint64_t detangleToleranceLow,
    uint64_t detangleToleranceHigh,
    bool useBayesianModel,
    double epsilon,
    double minLogP)
{
    if(debug) {
        cout << "Detangling vertices (general detangling)." << endl;
    }
    AssemblyGraph& cGraph = *this;

    vector<vertex_descriptor> allVertices;
    BGL_FORALL_VERTICES(cv, cGraph, AssemblyGraph) {
        allVertices.push_back(cv);
    }

    uint64_t detangledCount = 0;
    for(const vertex_descriptor cv: allVertices) {
        if(detangleVertexGeneral(cv, debug, detangleToleranceLow, detangleToleranceHigh,
            useBayesianModel, epsilon, minLogP)) {
            ++detangledCount;
        }
    }

    if(debug) {
        cout << "Detangled " << detangledCount << " vertices." << endl;

    }

    return detangledCount > 0;
}


// Compute the tangle matrix given in-edges and out-edges.
// The last bubble of each in-edge and the first bubble
// of each out-edge must be haploid.
void AssemblyGraph::computeTangleMatrix(
    const vector<edge_descriptor>& inEdges,
    const vector<edge_descriptor>& outEdges,
    vector< vector<uint64_t> >& tangleMatrix,
    bool setToZeroForComplementaryPairs
    ) const
{
    const AssemblyGraph& cGraph = *this;

    tangleMatrix.clear();
    tangleMatrix.resize(inEdges.size(), vector<uint64_t>(outEdges.size()));

    for(uint64_t i0=0; i0<inEdges.size(); i0++) {
        const edge_descriptor ce0 = inEdges[i0];
        const BubbleChain& bubbleChain0 = cGraph[ce0];
        const Bubble& bubble0 = bubbleChain0.lastBubble();
        SHASTA_ASSERT(bubble0.isHaploid());
        const Chain& chain0 = bubble0.front();
        SHASTA_ASSERT(chain0.size() >= 2);
        const MarkerGraphEdgeId markerGraphEdgeId0 = chain0[chain0.size() - 2];  // Exclude last

        for(uint64_t i1=0; i1<outEdges.size(); i1++) {
            const edge_descriptor ce1 = outEdges[i1];
            const BubbleChain& bubbleChain1 = cGraph[ce1];
            const Bubble& bubble1 = bubbleChain1.firstBubble();
            SHASTA_ASSERT(bubble1.isHaploid());
            const Chain& chain1 = bubble1.front();
            SHASTA_ASSERT(chain1.size() >= 2);
            const MarkerGraphEdgeId markerGraphEdgeId1 = chain1[1];  // Exclude first

            if(setToZeroForComplementaryPairs and
                assembler.markerGraph.reverseComplementEdge[markerGraphEdgeId0] == markerGraphEdgeId1) {
                tangleMatrix[i0][i1] = 0;
            } else {
                MarkerGraphEdgePairInfo info;
                SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(markerGraphEdgeId0, markerGraphEdgeId1, info));
                tangleMatrix[i0][i1] = info.common;
            }
        }
    }
}



#if 0
// This works if the following is true:
// - For all incoming edges (bubble chains) of cv, the last bubble is haploid.
// - For all outgoing edges (bubble chains) of cv, the first bubble is haploid.
bool AssemblyGraph::detangleVertexStrict(
    vertex_descriptor cv, bool debug)
{
    AssemblyGraph& cGraph = *this;

    // Gather the in-edges and check that the last bubble is haploid.
    vector<edge_descriptor> inEdges;
    BGL_FORALL_INEDGES(cv, ce, cGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.lastBubble().isHaploid()) {
            return false;
        }
        inEdges.push_back(ce);
    }

    // Gather the out-edges and check that the first bubble is haploid.
    vector<edge_descriptor> outEdges;
    BGL_FORALL_OUTEDGES(cv, ce, cGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.firstBubble().isHaploid()) {
            return false;
        }
        outEdges.push_back(ce);
    }

    if(inEdges.size() == 1 and outEdges.size() == 1) {
        return false;
    }

    // Compute the tangle matrix.
    vector< vector<uint64_t> > tangleMatrix;
    computeTangleMatrix(inEdges, outEdges, tangleMatrix, false);

    if(debug) {
        cout << "Tangle matrix for vertex " << cGraph[cv].edgeId << endl;
        for(uint64_t i0=0; i0<inEdges.size(); i0++) {
            for(uint64_t i1=0; i1<outEdges.size(); i1++) {
                cout << bubbleChainStringId(inEdges[i0]) << " " <<
                    bubbleChainStringId(outEdges[i1]) << " " <<
                    tangleMatrix[i0][i1] << endl;
            }
        }
    }

    // If the tangle matrix contains no zeros, there is nothing to do.
    bool foundZero = false;
    for(uint64_t i0=0; i0<inEdges.size(); i0++) {
        for(uint64_t i1=0; i1<outEdges.size(); i1++) {
            if(tangleMatrix[i0][i1] == 0) {
                foundZero = true;
                break;
            }
        }
        if(foundZero) {
            break;
        }
    }
    if(not foundZero) {
        return false;
    }

    // To avoid breaking contiguity, we require each column and each row of the
    // tangle matrix to have at least one non-zero element.
    // This means that each in-edge will be "merged" with at least one out-edge,
    // and each out-edge will be "merged" with at least one in-edge.
    for(uint64_t i0=0; i0<inEdges.size(); i0++) {
        bool foundNonZero = false;
        for(uint64_t i1=0; i1<outEdges.size(); i1++) {
            if(tangleMatrix[i0][i1] != 0) {
                foundNonZero = true;
                break;
            }
        }
        if(not foundNonZero) {
            return false;
        }
    }
    for(uint64_t i1=0; i1<outEdges.size(); i1++) {
        bool foundNonZero = false;
        for(uint64_t i0=0; i0<inEdges.size(); i0++) {
            if(tangleMatrix[i0][i1] != 0) {
                foundNonZero = true;
                break;
            }
        }
        if(not foundNonZero) {
            return false;
        }
    }

    if(debug) {
        cout << "This vertex will be detangled " << inEdges.size() << " by " << outEdges.size() << endl;
    }



    // Each non-zero element of the tangle matrix generates a new edge,
    // obtained by "merging" an in-edge with an out-edge.
    for(uint64_t i0=0; i0<inEdges.size(); i0++) {
        const edge_descriptor ce0 = inEdges[i0];
        const BubbleChain& bubbleChain0 = cGraph[ce0];
        const Bubble& bubble0 = bubbleChain0.lastBubble();
        SHASTA_ASSERT(bubble0.isHaploid());
        const Chain& chain0 = bubble0.front();
        SHASTA_ASSERT(chain0.size() >= 2);
        for(uint64_t i1=0; i1<outEdges.size(); i1++) {
            if(tangleMatrix[i0][i1] == 0) {
                continue;
            }
            const edge_descriptor ce1 = outEdges[i1];
            const BubbleChain& bubbleChain1 = cGraph[ce1];
            const Bubble& bubble1 = bubbleChain1.firstBubble();
            SHASTA_ASSERT(bubble1.isHaploid());
            const Chain& chain1 = bubble1.front();
            SHASTA_ASSERT(chain1.size() >= 2);

            edge_descriptor eNew;
            tie(eNew, ignore) = add_edge(source(ce0, cGraph), target(ce1, graph), cGraph);
            AssemblyGraphEdge& newEdge = cGraph[eNew];
            newEdge.id = nextEdgeId++;
            BubbleChain& newBubbleChain = newEdge;

            if(debug) {
                cout << "Merging " <<
                    bubbleChainStringId(ce0) << " " <<
                    bubbleChainStringId(ce1) << " into " <<
                    bubbleChainStringId(eNew) << endl;
            }

            // Create the new BubbleChain. It is obtained by joining
            // bubbleChain0 and bubbleChain1, with vertex cv
            // removed from the end of bubbleChain0
            // and from the beginning of bubbleChain1.
            // Here we use the above assumption that
            // the last bubble of bubbleChain0 and the first bubble of bubbleChain1
            // are haploid.
            newBubbleChain = bubbleChain0;

            // Remove cv from the end.
            Bubble& newBubbleLast = newBubbleChain.back();
            SHASTA_ASSERT(newBubbleLast.size() == 1);
            Chain& newChainLast = newBubbleLast.front();
            SHASTA_ASSERT(newChainLast.back() == cGraph[cv].edgeId);
            newChainLast.resize(newChainLast.size() - 1);

            // Append chain1, except for cv.
            SHASTA_ASSERT(chain1.front() == cGraph[cv].edgeId);
            copy(chain1.begin() + 1, chain1.end(), back_inserter(newChainLast));

            // Append the rest of bubbleChain1.
            copy(bubbleChain1.begin() + 1, bubbleChain1.end(), back_inserter(newBubbleChain));
        }

    }

    // Now we can remove cv and all of its in-edges and out-edges.
    clear_vertex(cv, cGraph);
    cGraph.removeVertex(cv);

    return true;
}
#endif



bool AssemblyGraph::detangleVertex(
    vertex_descriptor cv,
    bool debug,
    uint64_t detangleToleranceLow,
    uint64_t detangleToleranceHigh,
    bool useBayesianModel,
    double epsilon,
    double minLogP)
{
    AssemblyGraph& cGraph = *this;

    if(debug) {
        cout << "Attempting to detangle vertex " << cGraph[cv].edgeId << endl;
    }


    // Gather the in-edges and check that the last bubble is haploid.
    vector<edge_descriptor> inEdges;
    BGL_FORALL_INEDGES(cv, ce, cGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.lastBubble().isHaploid()) {
            if(debug) {
                cout << "Not detangled because the last bubble of in-edge " <<
                    bubbleChainStringId(ce) << " is not haploid." << endl;
            }
            return false;
        }
        inEdges.push_back(ce);
    }

    // Gather the out-edges and check that the first bubble is haploid.
    vector<edge_descriptor> outEdges;
    BGL_FORALL_OUTEDGES(cv, ce, cGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.firstBubble().isHaploid()) {
            if(debug) {
                cout << "Not detangled because the first bubble of out-edge " <<
                    bubbleChainStringId(ce) << " is not haploid." << endl;
            }
            return false;
        }
        outEdges.push_back(ce);
    }

    if(inEdges.size() == 0 or outEdges.size() == 0) {
        if(debug) {
            cout << "Not detangling due to degree (case 1)." << endl;
        }
        return false;
    }
    if(inEdges.size() < 2 and outEdges.size() < 2) {
        if(debug) {
            cout << "Not detangling due to degree (case 2)." << endl;
        }
        return false;
    }

    // Compute the tangle matrix.
    vector< vector<uint64_t> > tangleMatrix;
    computeTangleMatrix(inEdges, outEdges, tangleMatrix, false);

    if(debug) {
        cout << "Tangle matrix for vertex " << cGraph[cv].edgeId << endl;

        cout << "In-edges: ";
        for(const edge_descriptor ce: inEdges) {
            cout << " " << bubbleChainStringId(ce);
        }
        cout << endl;

        cout << "Out-edges: ";
        for(const edge_descriptor ce: outEdges) {
            cout << " " << bubbleChainStringId(ce);
        }
        cout << endl;

        for(uint64_t i0=0; i0<inEdges.size(); i0++) {
            for(uint64_t i1=0; i1<outEdges.size(); i1++) {
                cout << bubbleChainStringId(inEdges[i0]) << " " <<
                    bubbleChainStringId(outEdges[i1]) << " " <<
                    tangleMatrix[i0][i1] << endl;
            }
        }
    }



    // Do the detangling based on the tangle matrix.
    if(useBayesianModel and inEdges.size() == 2 and outEdges.size() == 2) {

        // Use the 2 by 2 Bayesian model for detangling.
        array< array<uint64_t, 2>, 2> tangleMatrix22;
        for(uint64_t i=0; i<2; i++) {
            for(uint64_t j=0; j<2; j++) {
                tangleMatrix22[i][j] = tangleMatrix[i][j];
            }
        }

        // Compute logarithmic probability ratio of in-phase and out-of-phase
        // against random.
        double logPin;
        double logPout;
        tie(logPin, logPout) = diploidBayesianPhase(tangleMatrix22, epsilon);
        if(debug) {
            cout << "logPin = " << logPin << ", logPout = " << logPout << endl;
        }

        // const bool isInPhase    = (logPin  >= minLogP) and ((logPin - logPout) >= minLogP);
        // const bool isOutOfPhase = (logPout >= minLogP) and ((logPout - logPin) >= minLogP);
        // Ignore the random hypothesis.
        const bool isInPhase    = (logPin - logPout) >= minLogP;
        const bool isOutOfPhase = (logPout - logPin) >= minLogP;

        if(isInPhase or isOutOfPhase) {

            // We can detangle.
            if(debug) {
                cout << "This vertex will be detangled." << endl;
            }

            // Create truncated versions of the inEdges and outEdges.
            vector<vertex_descriptor> inVertices;
            for(const edge_descriptor ce: inEdges) {
                inVertices.push_back(cloneAndTruncateAtEnd(ce));
            }
            vector<vertex_descriptor> outVertices;
            for(const edge_descriptor ce: outEdges) {
                outVertices.push_back(cloneAndTruncateAtBeginning(ce));
            }

            if(isInPhase) {
                connect(inVertices[0], outVertices[0]);
                connect(inVertices[1], outVertices[1]);
            } else {
                connect(inVertices[0], outVertices[1]);
                connect(inVertices[1], outVertices[0]);
            }

            // Now we can remove cv and all of its in-edges and out-edges.
            clear_vertex(cv, cGraph);
            cGraph.removeVertex(cv);
            return true;

        } else {

            // Ambiguous. Don't detangle.
            if(debug) {
                cout << "This vertex will not be detangled." << endl;
            }
            return false;
        }

    } else {

        // Don't use the Bayesian model.
        // Instead, do simple counting of tangle matrix elements.

        // Count the number of significant, ambiguous, and negligible elements
        // in the tangle matrix.
        uint64_t significantCount = 0;
        uint64_t ambiguousCount = 0;
        uint64_t negligibleCount = 0;
        for(uint64_t i0=0; i0<inEdges.size(); i0++) {
            for(uint64_t i1=0; i1<outEdges.size(); i1++) {
                const uint64_t t = tangleMatrix[i0][i1];
                if(t <= detangleToleranceLow) {
                    ++negligibleCount;
                } else if(t >= detangleToleranceHigh) {
                    ++significantCount;
                } else {
                    ++ambiguousCount;
                }
            }
        }

        // If the tangle matrix contains any ambiguous elements, do nothing.
        if(ambiguousCount > 0) {
            return false;
        }

        // There are no ambiguous elements.
        // If there are no negligible element, that is all elements of the tangle matrix are significant,
        // there is nothing to do.
        if(negligibleCount == 0) {
            return false;
        }

        // To avoid breaking contiguity, we require each column and each row of the
        // tangle matrix to have at least one significant element.
        // This means that each in-edge will be "merged" with at least one out-edge,
        // and each out-edge will be "merged" with at least one in-edge.
        for(uint64_t i0=0; i0<inEdges.size(); i0++) {
            bool foundSignificant = false;
            for(uint64_t i1=0; i1<outEdges.size(); i1++) {
                if(tangleMatrix[i0][i1] >= detangleToleranceHigh) {
                    foundSignificant = true;
                    break;
                }
            }
            if(not foundSignificant) {
                return false;
            }
        }
        for(uint64_t i1=0; i1<outEdges.size(); i1++) {
            bool foundSignificant = false;
            for(uint64_t i0=0; i0<inEdges.size(); i0++) {
                if(tangleMatrix[i0][i1] >= detangleToleranceHigh) {
                    foundSignificant = true;
                    break;
                }
            }
            if(not foundSignificant) {
                return false;
            }
        }

        if(debug) {
            cout << "This vertex will be detangled " << inEdges.size() << " by " << outEdges.size() << endl;
        }

        // Create truncated versions of the inEdges and outEdges.
        vector<vertex_descriptor> inVertices;
        for(const edge_descriptor ce: inEdges) {
            inVertices.push_back(cloneAndTruncateAtEnd(ce));
        }
        vector<vertex_descriptor> outVertices;
        for(const edge_descriptor ce: outEdges) {
            outVertices.push_back(cloneAndTruncateAtBeginning(ce));
        }

        // Each significant element of the tangle matrix generates a new edge.
        for(uint64_t i0=0; i0<inEdges.size(); i0++) {
            for(uint64_t i1=0; i1<outEdges.size(); i1++) {
                if(tangleMatrix[i0][i1] >= detangleToleranceHigh) {
                    connect(inVertices[i0], outVertices[i1]);
                }
            }
        }

        // Now we can remove cv and all of its in-edges and out-edges.
        clear_vertex(cv, cGraph);
        cGraph.removeVertex(cv);
        return true;
    }


#if 0
    // Each significant element of the tangle matrix generates a new edge,
    // obtained by "merging" an in-edge with an out-edge.
    for(uint64_t i0=0; i0<inEdges.size(); i0++) {
        const edge_descriptor ce0 = inEdges[i0];
        const BubbleChain& bubbleChain0 = cGraph[ce0];
        const Bubble& bubble0 = bubbleChain0.lastBubble();
        SHASTA_ASSERT(bubble0.isHaploid());
        const Chain& chain0 = bubble0.front();
        SHASTA_ASSERT(chain0.size() >= 2);
        for(uint64_t i1=0; i1<outEdges.size(); i1++) {
            if(tangleMatrix[i0][i1] < detangleToleranceHigh) {
                continue;
            }
            const edge_descriptor ce1 = outEdges[i1];
            const BubbleChain& bubbleChain1 = cGraph[ce1];
            const Bubble& bubble1 = bubbleChain1.firstBubble();
            SHASTA_ASSERT(bubble1.isHaploid());
            const Chain& chain1 = bubble1.front();
            SHASTA_ASSERT(chain1.size() >= 2);

            edge_descriptor eNew;
            tie(eNew, ignore) = add_edge(source(ce0, cGraph), target(ce1, graph), cGraph);
            AssemblyGraphEdge& newEdge = cGraph[eNew];
            newEdge.id = nextEdgeId++;
            BubbleChain& newBubbleChain = newEdge;

            if(debug) {
                cout << "Merging " <<
                    bubbleChainStringId(ce0) << " " <<
                    bubbleChainStringId(ce1) << " into " <<
                    bubbleChainStringId(eNew) << endl;
            }

            // Create the new BubbleChain. It is obtained by joining
            // bubbleChain0 and bubbleChain1, with vertex cv
            // removed from the end of bubbleChain0
            // and from the beginning of bubbleChain1.
            // Here we use the above assumption that
            // the last bubble of bubbleChain0 and the first bubble of bubbleChain1
            // are haploid.
            newBubbleChain = bubbleChain0;

            // Remove cv from the end.
            Bubble& newBubbleLast = newBubbleChain.back();
            SHASTA_ASSERT(newBubbleLast.size() == 1);
            Chain& newChainLast = newBubbleLast.front();
            SHASTA_ASSERT(newChainLast.back() == cGraph[cv].edgeId);
            newChainLast.resize(newChainLast.size() - 1);

            // Append chain1, except for cv.
            SHASTA_ASSERT(chain1.front() == cGraph[cv].edgeId);
            copy(chain1.begin() + 1, chain1.end(), back_inserter(newChainLast));

            // Append the rest of bubbleChain1.
            copy(bubbleChain1.begin() + 1, bubbleChain1.end(), back_inserter(newBubbleChain));
        }

    }
#endif


    SHASTA_ASSERT(0);
}



// Ths version can handle the case where the last bubble of an in-edge
// or the first bubble of an out-edge is not haploid.
// It works like this:
// - Compute a generalized tangle matrix taking using the next to last
//   MarkerGraphEdgeId of each incoming chain
//   and the second MarkerGraphEdgeId of each outgoing chain.
// - If detangling is possible based on this generalized tangle matrix,
//   split the last bubble of each incoming edge and the first
//   bubble of each outgoing edge. After this operation,
//   the last bubble of each in-edge is haploid and the first bubble
//   of each out-edge is haploid.
// - Call detangleVertex to do the detangling.
bool AssemblyGraph::detangleVertexGeneral(
    vertex_descriptor cv,
    bool debug,
    uint64_t detangleToleranceLow,
    uint64_t detangleToleranceHigh,
    bool useBayesianModel,
    double epsilon,
    double minLogP)
{
    AssemblyGraph& cGraph = *this;

#if 0
    // Use detangleVertex, if possible.
    bool involvesNonHaploidBubbles = false;
    BGL_FORALL_INEDGES(cv, ce, cGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.lastBubble().isHaploid()) {
            involvesNonHaploidBubbles = true;
        }
    }
    BGL_FORALL_OUTEDGES(cv, ce, cGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.firstBubble().isHaploid()) {
            involvesNonHaploidBubbles = true;
        }
    }
    if(not involvesNonHaploidBubbles) {
        if(debug) {
            cout << "No non-haploid bubbles involved, using detangleVertex." << endl;
        }
        return detangleVertex(cv, debug, detangleToleranceLow, detangleToleranceHigh);
    }
#endif

    if(in_degree(cv, cGraph) < 2 or out_degree(cv, cGraph) < 2) {
        return false;
    }

    if(debug) {
        cout << "Attempting general detangling for vertex " << cGraph[cv].edgeId << endl;
    }

    class ChainInfo {
    public:
        edge_descriptor ce;
        uint64_t indexInBubble;
        MarkerGraphEdgeId edgeId;
    };
    vector<ChainInfo> inChains;
    BGL_FORALL_INEDGES(cv, ce, cGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = cGraph[ce];
        const Bubble& bubble = bubbleChain.lastBubble();
        for(uint64_t indexInBubble=0; indexInBubble<bubble.size(); indexInBubble++) {
            const Chain& chain = bubble[indexInBubble];
            inChains.push_back({ce, indexInBubble, chain[chain.size() - 2]});
        }
    }
    vector<ChainInfo> outChains;
    BGL_FORALL_OUTEDGES(cv, ce, cGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = cGraph[ce];
        const Bubble& bubble = bubbleChain.firstBubble();
        for(uint64_t indexInBubble=0; indexInBubble<bubble.size(); indexInBubble++) {
            const Chain& chain = bubble[indexInBubble];
            outChains.push_back({ce, indexInBubble, chain[1]});
        }
    }

    if(debug) {

        cout << "In:" << endl;
        for(const ChainInfo& chainInfo: inChains) {
            cout << bubbleChainStringId(chainInfo.ce) << " " <<
                chainInfo.indexInBubble << " " <<
                chainInfo.edgeId << endl;
        }

        cout << "Out:" << endl;
        for(const ChainInfo& chainInfo: outChains) {
            cout << bubbleChainStringId(chainInfo.ce) << " " <<
                chainInfo.indexInBubble << " " <<
                chainInfo.edgeId << endl;
        }
    }

    // Compute a generalized tangle matrix.
    vector<vector<uint64_t> > tangleMatrix(inChains.size(), vector<uint64_t>(outChains.size()));
    for(uint64_t i0=0; i0<inChains.size(); i0++) {
        const MarkerGraphEdgeId markerGraphEdgeId0 = inChains[i0].edgeId;

        for(uint64_t i1=0; i1<outChains.size(); i1++) {
            const MarkerGraphEdgeId markerGraphEdgeId1 = outChains[i1].edgeId;

            MarkerGraphEdgePairInfo info;
            SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(markerGraphEdgeId0, markerGraphEdgeId1, info));
            tangleMatrix[i0][i1] = info.common;
        }
    }

    if(debug) {
        cout << "Tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inChains.size(); i0++) {
            const ChainInfo& chainInfo0 = inChains[i0];
            for(uint64_t i1=0; i1<outChains.size(); i1++) {
                const ChainInfo& chainInfo1 = outChains[i1];

                cout <<
                    bubbleChainStringId(chainInfo0.ce) << " " <<
                    chainInfo0.indexInBubble << " " <<
                    chainInfo0.edgeId << " " <<
                    bubbleChainStringId(chainInfo1.ce) << " " <<
                    chainInfo1.indexInBubble << " " <<
                    chainInfo1.edgeId << " " <<
                    tangleMatrix[i0][i1] << endl;
            }
        }

    }


    // Figure out if we can detangle.
    if(useBayesianModel and
        (inChains.size() == 2) and
        (outChains.size() == 2)) {

        // Use the 2 by 2 Bayesian model for detangling.
        array< array<uint64_t, 2>, 2> tangleMatrix22;
        for(uint64_t i=0; i<2; i++) {
            for(uint64_t j=0; j<2; j++) {
                tangleMatrix22[i][j] = tangleMatrix[i][j];
            }
        }

        // Compute logarithmic probability ratio of in-phase and out-of-phase
        // against random.
        double logPin;
        double logPout;
        tie(logPin, logPout) = diploidBayesianPhase(tangleMatrix22, epsilon);
        if(debug) {
            cout << "logPin = " << logPin << ", logPout = " << logPout << endl;
        }

        const bool isInPhase    = (logPin  >= minLogP) and ((logPin - logPout) >= minLogP);
        const bool isOutOfPhase = (logPout >= minLogP) and ((logPout - logPin) >= minLogP);
        if(not (isInPhase or isOutOfPhase)) {
            if(debug) {
                cout << "Ambiguous, don't detangle." << endl;
            }
            return false;
        }

    } else {

        // Not using the Bayesian model.
        // Count the number of significant, ambiguous, and negligible elements
        // in the tangle matrix.
        uint64_t significantCount = 0;
        uint64_t ambiguousCount = 0;
        uint64_t negligibleCount = 0;
        for(uint64_t i0=0; i0<inChains.size(); i0++) {
            for(uint64_t i1=0; i1<outChains.size(); i1++) {
                const uint64_t t = tangleMatrix[i0][i1];
                if(t <= detangleToleranceLow) {
                    ++negligibleCount;
                } else if(t >= detangleToleranceHigh) {
                    ++significantCount;
                } else {
                    ++ambiguousCount;
                }
            }
        }

        // If the tangle matrix contains any ambiguous elements, do nothing.
        if(ambiguousCount > 0) {
            if(debug) {
                cout << "Tangle matrix is ambiguous." << endl;
            }
            return false;
        }
        // There are no ambiguous elements.
        // If there are no negligible element, that is all elements of the tangle matrix are significant,
        // there is nothing to do.
        if(negligibleCount == 0) {
            return false;
        }

        // To avoid breaking contiguity, we require each column and each row of the
        // tangle matrix to have at least one significant element.
        // This means that each in-edge will be "merged" with at least one out-edge,
        // and each out-edge will be "merged" with at least one in-edge.
        for(uint64_t i0=0; i0<inChains.size(); i0++) {
            bool foundSignificant = false;
            for(uint64_t i1=0; i1<outChains.size(); i1++) {
                if(tangleMatrix[i0][i1] >= detangleToleranceHigh) {
                    foundSignificant = true;
                    break;
                }
            }
            if(not foundSignificant) {
                return false;
            }
        }
        for(uint64_t i1=0; i1<outChains.size(); i1++) {
            bool foundSignificant = false;
            for(uint64_t i0=0; i0<inChains.size(); i0++) {
                if(tangleMatrix[i0][i1] >= detangleToleranceHigh) {
                    foundSignificant = true;
                    break;
                }
            }
            if(not foundSignificant) {
                return false;
            }
        }
    }

    if(debug) {
        cout << "This vertex can be detangled after some splitting of bubble chains." << endl;
    }


    // Make sure the last bubble of all in-edges is haploid.
    in_edge_iterator itIn, itInEnd;
    tie(itIn, itInEnd) = in_edges(cv, cGraph);
    while(itIn != itInEnd) {
        const edge_descriptor ce = *itIn;
        ++itIn; // Increment before possibly removing this edge!
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.lastBubble().isHaploid()) {
            if(debug) {
                cout << "In-edge " << bubbleChainStringId(ce) <<
                    " needs to be split at the end." << endl;
            }
            splitBubbleChainAtEnd(ce);
        }
    }

    // Make sure the first bubble of all out-edges is haploid.
    out_edge_iterator itOut, itOutEnd;
    tie(itOut, itOutEnd) = out_edges(cv, cGraph);
    while(itOut != itOutEnd) {
        const edge_descriptor ce = *itOut;
        ++itOut; // Increment before possibly removing this edge!
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.firstBubble().isHaploid()) {
            if(debug) {
                cout << "Out-edge " << bubbleChainStringId(ce) <<
                    " needs to be split at the beginning." << endl;
            }
            splitBubbleChainAtBeginning(ce);
        }
    }

    // Now we can detangle using detangleVertex.
    if(debug) {
        cout << "Calling detangleVertex." << endl;
    }
    return detangleVertex(cv, debug, detangleToleranceLow, detangleToleranceHigh,
        useBayesianModel, epsilon, minLogP);
}



// Split the first bubble of a bubble chain.
// Used by detangleVertexGeneral to eliminate
// non-haploid bubbles adjacent to a vertex to be detangled.
void AssemblyGraph::splitBubbleChainAtBeginning(edge_descriptor ce)
{
    AssemblyGraph& cGraph = *this;

    const BubbleChain& bubbleChain = cGraph[ce];
    const Bubble& firstBubble = bubbleChain.firstBubble();
    SHASTA_ASSERT(not firstBubble.isHaploid());

    const vertex_descriptor cv0 = source(ce, cGraph);
    const vertex_descriptor cv1 = target(ce, cGraph);



    // General case where the bubble chain has more than one bubble.
    // Generate a new edge containing the bubble chain except for
    // the first bubble, plus one new edge for each chain in the firstBubble.
    if(bubbleChain.size() > 1) {

        // Generate one new edge containing the bubble chain except for
        // the first bubble.
        AssemblyGraphEdge newEdge;
        newEdge.id = nextEdgeId++;
        copy(bubbleChain.begin() + 1, bubbleChain.end(), back_inserter(newEdge));
        const vertex_descriptor cv2 = createVertex(newEdge.front().front().front());
        boost::add_edge(cv2, cv1, newEdge, cGraph);

        // Generate a  new edge for each chain in the firstBubble.
        for(const Chain& chain: firstBubble) {
            AssemblyGraphEdge newEdge;
            newEdge.resize(1);  // The new edge has only one bubble.
            Bubble& newBubble = newEdge.front();
            newEdge.id = nextEdgeId++;
            newBubble.push_back(chain);
            boost::add_edge(cv0, cv2, newEdge, cGraph);
        }
    }


    // Special case where the bubble chain has one bubble.
    // We generate one new edge for each chain in the firstBubble.
    else {

        // Generate a new edge for each chain in the firstBubble.
        for(const Chain& chain: firstBubble) {
            AssemblyGraphEdge newEdge;
            newEdge.resize(1);  // The new edge has only one bubble.
            Bubble& newBubble = newEdge.front();
            newEdge.id = nextEdgeId++;
            newBubble.push_back(chain);
            boost::add_edge(cv0, cv1, newEdge, cGraph);
        }
    }

    // Now we can remove the original bubble chain.
    boost::remove_edge(ce, cGraph);
}



// Split the last bubble of a bubble chain.
// Used by detangleVertexGeneral to eliminate
// non-haploid bubbles adjacent to a vertex to be detangled.
void AssemblyGraph::splitBubbleChainAtEnd(edge_descriptor ce)
{
    AssemblyGraph& cGraph = *this;

    const BubbleChain& bubbleChain = cGraph[ce];
    const Bubble& lastBubble = bubbleChain.lastBubble();
    SHASTA_ASSERT(not lastBubble.isHaploid());

    const vertex_descriptor cv0 = source(ce, cGraph);
    const vertex_descriptor cv1 = target(ce, cGraph);



    // General case where the bubble chain has more than one bubble.
    // Generate a new edge containing the bubble chain except for
    // the last bubble, plus one new edge for each chain in the lastBubble.
    if(bubbleChain.size() > 1) {

        // Generate one new edge containing the bubble chain except for
        // the last bubble.
        AssemblyGraphEdge newEdge;
        newEdge.id = nextEdgeId++;
        copy(bubbleChain.begin(), bubbleChain.end()-1, back_inserter(newEdge));
        const vertex_descriptor cv2 = createVertex(newEdge.back().front().back());
        boost::add_edge(cv0, cv2, newEdge, cGraph);

        // Generate a  new edge for each chain in the lastBubble.
        for(const Chain& chain: lastBubble) {
            AssemblyGraphEdge newEdge;
            newEdge.resize(1);  // The new edge has only one bubble.
            Bubble& newBubble = newEdge.front();
            newEdge.id = nextEdgeId++;
            newBubble.push_back(chain);
            boost::add_edge(cv2, cv1, newEdge, cGraph);
        }
    }


    // Special case where the bubble chain has one bubble.
    // We generate one new edge for each chain in the lastBubble.
    else {

        // Generate a new edge for each chain in the lastBubble.
        for(const Chain& chain: lastBubble) {
            AssemblyGraphEdge newEdge;
            newEdge.resize(1);  // The new edge has only one bubble.
            Bubble& newBubble = newEdge.front();
            newEdge.id = nextEdgeId++;
            newBubble.push_back(chain);
            boost::add_edge(cv0, cv1, newEdge, cGraph);
        }
    }

    // Now we can remove the original bubble chain.
    boost::remove_edge(ce, cGraph);
}



bool AssemblyGraph::detangleEdges(
    bool debug,
    uint64_t detangleToleranceLow,
    uint64_t detangleToleranceHigh,
    bool useBayesianModel,
    double epsilon,
    double minLogP)
{
    if(debug) {
        cout << "Detangling edges." << endl;
    }

    AssemblyGraph& cGraph = *this;

    // To safely iterate over edges while removing edges we must use edge ids
    // as unique identifiers, because edge descriptors can be reused as edges are
    // deleted ndw new edges are created.
    std::map<uint64_t, edge_descriptor> edgeMap;
    BGL_FORALL_EDGES(ce, cGraph, AssemblyGraph) {
        edgeMap.insert({cGraph[ce].id, ce});
    }

    uint64_t detangleCount = 0;;
    for(auto it=edgeMap.begin(); it!=edgeMap.end(); /* Incremented safely by detangleEdgeStrict */) {
        if(detangleEdge(debug, edgeMap, it, detangleToleranceLow, detangleToleranceHigh,
            useBayesianModel, epsilon, minLogP)) {
            ++detangleCount;
        }
    }

    if(debug) {
        cout << "Detangled " << detangleCount << " edges." << endl;
    }

    return detangleCount > 0;
}



bool AssemblyGraph::detangleEdgesGeneral(
    bool debug,
    uint64_t detangleToleranceLow,
    uint64_t detangleToleranceHigh,
    bool useBayesianModel,
    double epsilon,
    double minLogP)
{
    if(debug) {
        cout << "Detangling edges." << endl;
    }

    AssemblyGraph& cGraph = *this;

    // To safely iterate over edges while removing edges we must use edge ids
    // as unique identifiers, because edge descriptors can be reused as edges are
    // deleted ndw new edges are created.
    std::map<uint64_t, edge_descriptor> edgeMap;
    BGL_FORALL_EDGES(ce, cGraph, AssemblyGraph) {
        edgeMap.insert({cGraph[ce].id, ce});
    }

    uint64_t detangleCount = 0;;
    for(auto it=edgeMap.begin(); it!=edgeMap.end(); /* Incremented safely by detangleEdgeStrict */) {
        if(detangleEdgeGeneral(debug, edgeMap, it, detangleToleranceLow, detangleToleranceHigh,
            useBayesianModel, epsilon, minLogP)) {
            ++detangleCount;
        }
    }

    if(debug) {
        cout << "Detangled " << detangleCount << " edges." << endl;
    }

    return detangleCount > 0;
}



bool AssemblyGraph::detangleEdge(
    bool debug,
    std::map<uint64_t, edge_descriptor>& edgeMap,
    std::map<uint64_t, edge_descriptor>::iterator& it,
    uint64_t detangleToleranceLow,
    uint64_t detangleToleranceHigh,
    bool useBayesianModel,
    double epsilon,
    double minLogP)
{
    AssemblyGraph& cGraph = *this;
    const edge_descriptor ce = it->second;
    ++it;
    // edgeMap.erase(cGraph[ce].id);

    // Only try detangling if the edge consists of a single haploid bubble.
    // Otherwise detangling would lose information.
    BubbleChain& bubbleChain = cGraph[ce];
    if(bubbleChain.size() > 1) {
        return false;
    }
    if(bubbleChain.front().size() > 1) {
        return false;
    }

    // Tangle matrix elements <= detangleToleranceLow are treated as negigible.
    // Tangle matrix elements >= detangleToleranceHigh are treated as significant.
    // Tangle matrix elements in between are considered ambiguous.
    SHASTA_ASSERT(detangleToleranceHigh > detangleToleranceLow);

    const vertex_descriptor cv0 = source(ce, cGraph);
    const vertex_descriptor cv1 = target(ce, cGraph);

    if(out_degree(cv0, cGraph) != 1) {
        return false;
    }
    if(in_degree(cv1, cGraph) != 1) {
        return false;
    }

    if(debug) {
        cout << "Attempting to detangle edge " << bubbleChainStringId(ce) << endl;
    }

    // Gather the in-edges and check that the last bubble is haploid.
    // Ignore in-edges coming from cv1 (back-edges).
    vector<edge_descriptor> inEdges;
    vector<edge_descriptor> backEdges;
    BGL_FORALL_INEDGES(cv0, ce, cGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.lastBubble().isHaploid()) {
            if(debug) {
                cout << "Not detangling because the last bubble of in-edge " <<
                    bubbleChainStringId(ce) << " is not haploid." << endl;
            }
            return false;
        }
        if(source(ce, cGraph) != cv1) {
            inEdges.push_back(ce);
        } else {
            backEdges.push_back(ce);
        }
    }

    // Gather the out-edges and check that the first bubble is haploid.
    // Ignore out-edges going to cv0 (back-edges).
    vector<edge_descriptor> outEdges;
    BGL_FORALL_OUTEDGES(cv1, ce, cGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.firstBubble().isHaploid()) {
            if(debug) {
                cout << "Not detangling because the first bubble of out-edge " <<
                    bubbleChainStringId(ce) << " is not haploid." << endl;
            }
            return false;
        }
        if(target(ce, cGraph) != cv0) {
            outEdges.push_back(ce);
        }
    }

    if(inEdges.size() == 0 or outEdges.size() == 0) {
        if(debug) {
            cout << "Not detangling due to degree (case 1)." << endl;
        }
        return false;
    }
    if(inEdges.size() < 2 and outEdges.size() < 2) {
        if(debug) {
            cout << "Not detangling due to degree (case 2)." << endl;
        }
        return false;
    }
    if(inEdges.size() != outEdges.size()) {
        if(debug) {
            cout << "Not detangling due to degree (case 3)." << endl;
        }
        return false;
    }

    // Compute the tangle matrix.
    vector< vector<uint64_t> > tangleMatrix;
    computeTangleMatrix(inEdges, outEdges, tangleMatrix, false);

    if(debug) {
        cout << "Computing tangle matrix for edge " << bubbleChainStringId(ce) << endl;

        cout << "In-edges: ";
        for(const edge_descriptor ce: inEdges) {
            cout << " " << bubbleChainStringId(ce);
        }
        cout << endl;

        cout << "Out-edges: ";
        for(const edge_descriptor ce: outEdges) {
            cout << " " << bubbleChainStringId(ce);
        }
        cout << endl;

        cout << "Tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inEdges.size(); i0++) {
            const edge_descriptor ce0 = inEdges[i0];
            for(uint64_t i1=0; i1<outEdges.size(); i1++) {
                const edge_descriptor ce1 = outEdges[i1];
                cout <<
                    bubbleChainStringId(ce0) << " " <<
                    bubbleChainStringId(ce1) << " " <<
                    tangleMatrix[i0][i1];
                if(tangleMatrix[i0][i1] == 0) {
                    cout << " zero tangle matrix element";
                }
                cout << endl;
            }
        }
    }



    // Detangle based on the contents of the tangle matrix.
    if(useBayesianModel and inEdges.size() == 2 and outEdges.size() == 2) {

        // Use the 2 by 2 Bayesian model for detangling.
        array< array<uint64_t, 2>, 2> tangleMatrix22;
        for(uint64_t i=0; i<2; i++) {
            for(uint64_t j=0; j<2; j++) {
                tangleMatrix22[i][j] = tangleMatrix[i][j];
            }
        }

        // Compute logarithmic probability ratio of in-phase and out-of-phase
        // against random.
        double logPin;
        double logPout;
        tie(logPin, logPout) = diploidBayesianPhase(tangleMatrix22, epsilon);
        if(debug) {
            cout << "logPin = " << logPin << ", logPout = " << logPout << endl;
        }

        // const bool isInPhase    = (logPin  >= minLogP) and ((logPin - logPout) >= minLogP);
        // const bool isOutOfPhase = (logPout >= minLogP) and ((logPout - logPin) >= minLogP);
        // Ignore the random hypothesis.
        const bool isInPhase    = (logPin - logPout) >= minLogP;
        const bool isOutOfPhase = (logPout - logPin) >= minLogP;

        if(isInPhase or isOutOfPhase) {

            // We can detangle.

            // Create truncated versions of the inEdges and outEdges.
            vector<vertex_descriptor> inVertices;
            for(const edge_descriptor ce: inEdges) {
                inVertices.push_back(cloneAndTruncateAtEnd(ce));
            }
            vector<vertex_descriptor> outVertices;
            for(const edge_descriptor ce: outEdges) {
                outVertices.push_back(cloneAndTruncateAtBeginning(ce));
            }

            if(isInPhase) {
                const edge_descriptor e0 = connect(inVertices[0], outVertices[0]);
                const edge_descriptor e1 = connect(inVertices[1], outVertices[1]);
                if(debug) {
                    cout << "In phase: created " << bubbleChainStringId(e0) << " and " <<
                        bubbleChainStringId(e1) << endl;
                }
            } else {
                const edge_descriptor e0 = connect(inVertices[0], outVertices[1]);
                const edge_descriptor e1 = connect(inVertices[1], outVertices[0]);
                if(debug) {
                    cout << "Out of phase phase: created " << bubbleChainStringId(e0) << " and " <<
                        bubbleChainStringId(e1) << endl;
                }
            }

        } else {

            // Ambiguous. Don't detangle.
            if(debug) {
                cout << "Ambiguous. NOt detangling." << endl;
            }
            return false;
        }

    } else {



        // We are not using the Bayesian model.

        // Count the number of significant, ambiguous, and negligible elements
        // in the tangle matrix.
        uint64_t significantCount = 0;
        uint64_t ambiguousCount = 0;
        uint64_t negligibleCount = 0;
        for(uint64_t i0=0; i0<inEdges.size(); i0++) {
            for(uint64_t i1=0; i1<outEdges.size(); i1++) {
                const uint64_t t = tangleMatrix[i0][i1];
                if(t <= detangleToleranceLow) {
                    ++negligibleCount;
                } else if(t >= detangleToleranceHigh) {
                    ++significantCount;
                } else {
                    ++ambiguousCount;
                }
            }
        }

        // If the tangle matrix contains any ambiguous elements, do nothing.
        if(ambiguousCount > 0) {
            return false;
        }

        // There are no ambiguous elements.
        // If there are no negligible element, that is all elements of the tangle matrix are significant,
        // there is nothing to do.
        if(negligibleCount == 0) {
            return false;
        }

        // To avoid breaking contiguity, we require each column and each row of the
        // tangle matrix to have at least one significant element.
        // This means that each in-edge will be "merged" with at least one out-edge,
        // and each out-edge will be "merged" with at least one in-edge.
        // ACTUALY, FOR MORE ROBUSTNESS REQUIRE EXACTLY OEN SIGNIFICANT ELEMENT PER ROW AND COLUMN.
        for(uint64_t i0=0; i0<inEdges.size(); i0++) {
            uint64_t significantCount = 0;
            for(uint64_t i1=0; i1<outEdges.size(); i1++) {
                if(tangleMatrix[i0][i1] >= detangleToleranceHigh) {
                    ++significantCount;
                }
            }
            if(significantCount != 1) {
                return false;
            }
        }
        for(uint64_t i1=0; i1<outEdges.size(); i1++) {
            uint64_t significantCount = 0;
            for(uint64_t i0=0; i0<inEdges.size(); i0++) {
                if(tangleMatrix[i0][i1] >= detangleToleranceHigh) {
                    ++significantCount;
                }
            }
            if(significantCount != 1) {
                return false;
            }
        }

    #if 0
        // In an in-edge is also an out-edge, don't detangle.
        for(const edge_descriptor ce: inEdges) {
            if(find(outEdges.begin(), outEdges.end(), ce) != outEdges.end()) {
                if(debug) {
                    cout << "Not degangled because an in-edge is also an out-edge." << endl;
                }
                return false;
            }
        }
    #endif

        if(debug) {
            cout << "This edge will be detangled " << inEdges.size() << " by " << outEdges.size() << endl;
        }

        // Create truncated versions of the inEdges and outEdges.
        vector<vertex_descriptor> inVertices;
        for(const edge_descriptor ce: inEdges) {
            inVertices.push_back(cloneAndTruncateAtEnd(ce));
        }
        vector<vertex_descriptor> outVertices;
        for(const edge_descriptor ce: outEdges) {
            outVertices.push_back(cloneAndTruncateAtBeginning(ce));
        }


        // Each significant element of the tangle matrix generates a new edge.
        for(uint64_t i0=0; i0<inEdges.size(); i0++) {
            for(uint64_t i1=0; i1<outEdges.size(); i1++) {
                if(tangleMatrix[i0][i1] >= detangleToleranceHigh) {
                    const edge_descriptor ceNew = connect(inVertices[i0], outVertices[i1]);
                    if(debug) {
                        cout << "Created " << bubbleChainStringId(ceNew) << endl;
                    }
                }
            }
        }
    }


    // Now we can remove cv0, cv1, ce, and all of the in-edges and out-edges.
    // We have to do this while safely incrementing the edge iterator to point to the
    // next edge that was not removed.
    // We already incremented the iterator to point past ce.
    boost::remove_edge(ce, cGraph);
    for(const edge_descriptor ce: inEdges) {
        if(it != edgeMap.end() and cGraph[ce].id == it->first) {
            ++it;
        }
        edgeMap.erase(cGraph[ce].id);
        boost::remove_edge(ce, cGraph);
    }
    for(const edge_descriptor ce: outEdges) {
        if(it != edgeMap.end() and cGraph[ce].id == it->first) {
            ++it;
        }
        edgeMap.erase(cGraph[ce].id);
        boost::remove_edge(ce, cGraph);
    }
    for(const edge_descriptor ce: backEdges) {
        if(it != edgeMap.end() and cGraph[ce].id == it->first) {
            ++it;
        }
        edgeMap.erase(cGraph[ce].id);
        boost::remove_edge(ce, cGraph);
    }
    cGraph.removeVertex(cv0);
    cGraph.removeVertex(cv1);

    return true;
}



bool AssemblyGraph::detangleEdgeGeneral(
    bool debug,
    std::map<uint64_t, edge_descriptor>& edgeMap,
    std::map<uint64_t, edge_descriptor>::iterator& it,
    uint64_t detangleToleranceLow,
    uint64_t detangleToleranceHigh,
    bool useBayesianModel,
    double epsilon,
    double minLogP)
{
    // detangleEdgeGeneral does not have code to use the Bayesian model
    // for the 2 by 2 case. See detangleEdge.
    SHASTA_ASSERT(not useBayesianModel);

    AssemblyGraph& cGraph = *this;
    const edge_descriptor ce = it->second;
    ++it;

    // Tangle matrix elements <= detangleToleranceLow are treated as negigible.
    // Tangle matrix elements >= detangleToleranceHigh are treated as significant.
    // Tangle matrix elements in between are considered ambiguous.
    SHASTA_ASSERT(detangleToleranceHigh > detangleToleranceLow);

    const vertex_descriptor cv0 = source(ce, cGraph);
    const vertex_descriptor cv1 = target(ce, cGraph);

    if(out_degree(cv0, cGraph) != 1) {
        return false;
    }
    if(in_degree(cv1, cGraph) != 1) {
        return false;
    }

    if(debug) {
        cout << "Attempting general detangling of edge " << bubbleChainStringId(ce) << endl;
    }

    class ChainInfo {
    public:
        edge_descriptor ce;
        uint64_t indexInBubble;
        MarkerGraphEdgeId edgeId;
    };
    vector<ChainInfo> inChains;
    BGL_FORALL_INEDGES(cv0, ce, cGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = cGraph[ce];
        const Bubble& bubble = bubbleChain.lastBubble();
        for(uint64_t indexInBubble=0; indexInBubble<bubble.size(); indexInBubble++) {
            const Chain& chain = bubble[indexInBubble];
            inChains.push_back({ce, indexInBubble, chain[chain.size() - 2]});
        }
    }
    vector<ChainInfo> outChains;
    BGL_FORALL_OUTEDGES(cv1, ce, cGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = cGraph[ce];
        const Bubble& bubble = bubbleChain.firstBubble();
        for(uint64_t indexInBubble=0; indexInBubble<bubble.size(); indexInBubble++) {
            const Chain& chain = bubble[indexInBubble];
            outChains.push_back({ce, indexInBubble, chain[1]});
        }
    }

    if(debug) {

        cout << "In:" << endl;
        for(const ChainInfo& chainInfo: inChains) {
            cout << bubbleChainStringId(chainInfo.ce) << " " <<
                chainInfo.indexInBubble << " " <<
                chainInfo.edgeId << endl;
        }

        cout << "Out:" << endl;
        for(const ChainInfo& chainInfo: outChains) {
            cout << bubbleChainStringId(chainInfo.ce) << " " <<
                chainInfo.indexInBubble << " " <<
                chainInfo.edgeId << endl;
        }
    }

    if(inChains.size() != outChains.size()) {
        if(debug) {
            cout << "Not detangling due to degree." << endl;
        }
        return false;
    }


    // Compute a generalized tangle matrix.
    vector<vector<uint64_t> > tangleMatrix(inChains.size(), vector<uint64_t>(outChains.size()));
    for(uint64_t i0=0; i0<inChains.size(); i0++) {
        const MarkerGraphEdgeId markerGraphEdgeId0 = inChains[i0].edgeId;

        for(uint64_t i1=0; i1<outChains.size(); i1++) {
            const MarkerGraphEdgeId markerGraphEdgeId1 = outChains[i1].edgeId;

            MarkerGraphEdgePairInfo info;
            SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(markerGraphEdgeId0, markerGraphEdgeId1, info));
            tangleMatrix[i0][i1] = info.common;
        }
    }

    if(debug) {
        cout << "Tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inChains.size(); i0++) {
            const ChainInfo& chainInfo0 = inChains[i0];
            for(uint64_t i1=0; i1<outChains.size(); i1++) {
                const ChainInfo& chainInfo1 = outChains[i1];

                cout <<
                    bubbleChainStringId(chainInfo0.ce) << " " <<
                    chainInfo0.indexInBubble << " " <<
                    chainInfo0.edgeId << " " <<
                    bubbleChainStringId(chainInfo1.ce) << " " <<
                    chainInfo1.indexInBubble << " " <<
                    chainInfo1.edgeId << " " <<
                    tangleMatrix[i0][i1] << endl;
            }
        }

    }

    // Count the number of significant, ambiguous, and negligible elements
    // in the tangle matrix.
    uint64_t significantCount = 0;
    uint64_t ambiguousCount = 0;
    uint64_t negligibleCount = 0;
    for(uint64_t i0=0; i0<inChains.size(); i0++) {
        for(uint64_t i1=0; i1<outChains.size(); i1++) {
            const uint64_t t = tangleMatrix[i0][i1];
            if(t <= detangleToleranceLow) {
                ++negligibleCount;
            } else if(t >= detangleToleranceHigh) {
                ++significantCount;
            } else {
                ++ambiguousCount;
            }
        }
    }

    // If the tangle matrix contains any ambiguous elements, do nothing.
    if(ambiguousCount > 0) {
        return false;
    }

    // There are no ambiguous elements.
    // If there are no negligible element, that is all elements of the tangle matrix are significant,
    // there is nothing to do.
    if(negligibleCount == 0) {
        return false;
    }

    // To avoid breaking contiguity, we require each column and each row of the
    // tangle matrix to have at least one significant element.
    // This means that each in-edge will be "merged" with at least one out-edge,
    // and each out-edge will be "merged" with at least one in-edge.
    for(uint64_t i0=0; i0<inChains.size(); i0++) {
        bool foundSignificant = false;
        for(uint64_t i1=0; i1<outChains.size(); i1++) {
            if(tangleMatrix[i0][i1] >= detangleToleranceHigh) {
                foundSignificant = true;
                break;
            }
        }
        if(not foundSignificant) {
            return false;
        }
    }
    for(uint64_t i1=0; i1<outChains.size(); i1++) {
        bool foundSignificant = false;
        for(uint64_t i0=0; i0<inChains.size(); i0++) {
            if(tangleMatrix[i0][i1] >= detangleToleranceHigh) {
                foundSignificant = true;
                break;
            }
        }
        if(not foundSignificant) {
            return false;
        }
    }

    if(debug) {
        cout << "This edge can be detangled after some splitting of bubble chains." << endl;
    }

    // Make sure the last bubble of all in-edges is haploid.
    in_edge_iterator itIn, itInEnd;
    tie(itIn, itInEnd) = in_edges(cv0, cGraph);
    while(itIn != itInEnd) {
        const edge_descriptor ce = *itIn;
        ++itIn; // Increment before possibly removing this edge!
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.lastBubble().isHaploid()) {
            if(debug) {
                cout << "In-edge " << bubbleChainStringId(ce) <<
                    " needs to be split at the end." << endl;
            }
            splitBubbleChainAtEnd(ce);
        }
    }

    // Make sure the first bubble of all out-edges is haploid.
    out_edge_iterator itOut, itOutEnd;
    tie(itOut, itOutEnd) = out_edges(cv1, cGraph);
    while(itOut != itOutEnd) {
        const edge_descriptor ce = *itOut;
        ++itOut; // Increment before possibly removing this edge!
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.firstBubble().isHaploid()) {
            if(debug) {
                cout << "Out-edge " << bubbleChainStringId(ce) <<
                    " needs to be split at the beginning." << endl;
            }
            splitBubbleChainAtBeginning(ce);
        }
    }

    // Now we can detangle using detangleEdge.
    if(debug) {
        cout << "Calling detangleEdge." << endl;
    }
    --it;   // Because detangleEdge increments it again.
    return detangleEdge(debug, edgeMap, it, detangleToleranceLow, detangleToleranceHigh,
        useBayesianModel, epsilon, minLogP);
}


#if 0
bool AssemblyGraph::detangleEdgesWithSearch(
    bool debug,
    uint64_t detangleToleranceLow,
    uint64_t detangleToleranceHigh)
{
    if(debug) {
        cout << "Detangling edges with search." << endl;
    }

    AssemblyGraph& cGraph = *this;

    // To safely iterate over edges while removing edges we must use edge ids
    // as unique identifiers, because edge descriptors can be reused as edges are
    // deleted ndw new edges are created.
    std::map<uint64_t, edge_descriptor> edgeMap;
    BGL_FORALL_EDGES(ce, cGraph, AssemblyGraph) {
        edgeMap.insert({cGraph[ce].id, ce});
    }

    uint64_t detangleCount = 0;;
    for(auto it=edgeMap.begin(); it!=edgeMap.end(); /* Incremented safely by detangleEdgeStrict */) {
        if(detangleEdgeWithSearch(debug, edgeMap, it, detangleToleranceLow, detangleToleranceHigh)) {
            ++detangleCount;
        }
    }

    if(debug) {
        cout << "Detangled " << detangleCount << " edges." << endl;
    }

    return detangleCount > 0;
}



bool AssemblyGraph::detangleEdgeWithSearch(
    bool debug,
    std::map<uint64_t, edge_descriptor>& edgeMap,
    std::map<uint64_t, edge_descriptor>::iterator& it,
    uint64_t detangleToleranceLow,
    uint64_t detangleToleranceHigh)
{
    AssemblyGraph& cGraph = *this;
    const edge_descriptor ce = it->second;
    ++it;

    // Only try detangling if the edge consists of a single haploid bubble.
    // Otherwise detangling would lose information.
    BubbleChain& bubbleChain = cGraph[ce];
    if(bubbleChain.size() > 1) {
        return false;
    }
    if(bubbleChain.front().size() > 1) {
        return false;
    }


    const vertex_descriptor cv0 = source(ce, cGraph);
    const vertex_descriptor cv1 = target(ce, cGraph);

    if(out_degree(cv0, cGraph) != 1) {
        return false;
    }
    if(in_degree(cv1, cGraph) != 1) {
        return false;
    }

    if(debug) {
        cout << "Attempting to detangle edge " << bubbleChainStringId(ce) << " with search." << endl;
    }

    // Gather the in-edges and check that the last bubble is haploid.
    vector<edge_descriptor> inEdges;
    vector<edge_descriptor> backEdges;
    BGL_FORALL_INEDGES(cv0, ce, cGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.lastBubble().isHaploid()) {
            if(debug) {
                cout << "Not detangling because the last bubble of in-edge " <<
                    bubbleChainStringId(ce) << " is not haploid." << endl;
            }
            return false;
        }
        if(source(ce, cGraph) != cv1) {
            inEdges.push_back(ce);
        } else {
            backEdges.push_back(ce);
        }
    }

    // Gather the out-edges and check that the first bubble is haploid.
    // Ignore out-edges going to cv0 (back-edges).
    vector<edge_descriptor> outEdges;
    BGL_FORALL_OUTEDGES(cv1, ce, cGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.firstBubble().isHaploid()) {
            if(debug) {
                cout << "Not detangling because the first bubble of out-edge " <<
                    bubbleChainStringId(ce) << " is not haploid." << endl;
            }
            return false;
        }
        if(target(ce, cGraph) != cv0) {
            outEdges.push_back(ce);
        }
    }

    if(inEdges.size() == 0 or outEdges.size() == 0) {
        if(debug) {
            cout << "Not detangling due to degree (case 1)." << endl;
        }
        return false;
    }
    if(inEdges.size() != 2 and outEdges.size() != 2) {
        if(debug) {
            cout << "Not detangling due to degree (case 2)." << endl;
        }
        return false;
    }
    if(inEdges.size() != outEdges.size()) {
        if(debug) {
            cout << "Not detangling due to degree (case 3)." << endl;
        }
        return false;
    }


    // Get the second to last MarkerGraphEdgeIds of the incoming chains.
    array<MarkerGraphEdgeId, 2> in;
    for(uint64_t i=0; i<2; i++) {
        const Chain& chain = cGraph[inEdges[i]].back().front();
        in[i] = chain.secondToLast();
    }

    // Get the second MarkerGraphEdgeIds of the outgoing chains.
    array<MarkerGraphEdgeId, 2> out;
    for(uint64_t i=0; i<2; i++) {
        const Chain& chain = cGraph[outEdges[i]].front().front();
        out[i] = chain.second();
    }
    if(debug) {
        cout << "in " << bubbleChainStringId(inEdges[0]) << " " << bubbleChainStringId(inEdges[1]) << endl;
        cout << "out " << bubbleChainStringId(outEdges[0]) << " " << bubbleChainStringId(outEdges[1]) << endl;
        cout << "in " << in[0] << " " << in[1] << endl;
        cout << "out " << out[0] << " " << out[1] << endl;
    }

    array<array<vector<MarkerGraphEdgeId>, 2>, 2> detanglingCandidates;
    GlobalPathGraph::searchForDetangling(
        in, out,
        detangleToleranceHigh, detangleToleranceLow,
        assembler, detanglingCandidates);
    for(uint64_t i0=0; i0<2; i0++) {
        for(uint64_t i1=0; i1<2; i1++) {
            const auto& hits = detanglingCandidates[i0][i1];
            cout << "Found " << hits.size() << " hits for " << i0 << " " << i1 << ":" << endl;
            if(not hits.empty()) {
                copy(hits.begin(), hits.end(), ostream_iterator<MarkerGraphEdgeId>(cout, " "));
                cout << endl;
            }
        }
    }

    return false;

#if 0
    // Compute the tangle matrix.
    vector< vector<uint64_t> > tangleMatrix;
    computeTangleMatrix(inEdges, outEdges, tangleMatrix, false);

    if(debug) {
        cout << "Computing tangle matrix for edge " << bubbleChainStringId(ce) << endl;

        cout << "In-edges: ";
        for(const edge_descriptor ce: inEdges) {
            cout << " " << bubbleChainStringId(ce);
        }
        cout << endl;

        cout << "Out-edges: ";
        for(const edge_descriptor ce: outEdges) {
            cout << " " << bubbleChainStringId(ce);
        }
        cout << endl;

        cout << "Tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inEdges.size(); i0++) {
            const edge_descriptor ce0 = inEdges[i0];
            for(uint64_t i1=0; i1<outEdges.size(); i1++) {
                const edge_descriptor ce1 = outEdges[i1];
                cout <<
                    bubbleChainStringId(ce0) << " " <<
                    bubbleChainStringId(ce1) << " " <<
                    tangleMatrix[i0][i1];
                if(tangleMatrix[i0][i1] == 0) {
                    cout << " zero tangle matrix element";
                }
                cout << endl;
            }
        }
    }



    // Detangle based on the contents of the tangle matrix.
    if(useBayesianModel and inEdges.size() == 2 and outEdges.size() == 2) {

        // Use the 2 by 2 Bayesian model for detangling.
        array< array<uint64_t, 2>, 2> tangleMatrix22;
        for(uint64_t i=0; i<2; i++) {
            for(uint64_t j=0; j<2; j++) {
                tangleMatrix22[i][j] = tangleMatrix[i][j];
            }
        }

        // Compute logarithmic probability ratio of in-phase and out-of-phase
        // against random.
        double logPin;
        double logPout;
        tie(logPin, logPout) = diploidBayesianPhase(tangleMatrix22, epsilon);
        if(debug) {
            cout << "logPin = " << logPin << ", logPout = " << logPout << endl;
        }

        // const bool isInPhase    = (logPin  >= minLogP) and ((logPin - logPout) >= minLogP);
        // const bool isOutOfPhase = (logPout >= minLogP) and ((logPout - logPin) >= minLogP);
        // Ignore the random hypothesis.
        const bool isInPhase    = (logPin - logPout) >= minLogP;
        const bool isOutOfPhase = (logPout - logPin) >= minLogP;

        if(isInPhase or isOutOfPhase) {

            // We can detangle.

            // Create truncated versions of the inEdges and outEdges.
            vector<vertex_descriptor> inVertices;
            for(const edge_descriptor ce: inEdges) {
                inVertices.push_back(cloneAndTruncateAtEnd(ce));
            }
            vector<vertex_descriptor> outVertices;
            for(const edge_descriptor ce: outEdges) {
                outVertices.push_back(cloneAndTruncateAtBeginning(ce));
            }

            if(isInPhase) {
                const edge_descriptor e0 = connect(inVertices[0], outVertices[0]);
                const edge_descriptor e1 = connect(inVertices[1], outVertices[1]);
                if(debug) {
                    cout << "In phase: created " << bubbleChainStringId(e0) << " and " <<
                        bubbleChainStringId(e1) << endl;
                }
            } else {
                const edge_descriptor e0 = connect(inVertices[0], outVertices[1]);
                const edge_descriptor e1 = connect(inVertices[1], outVertices[0]);
                if(debug) {
                    cout << "Out of phase phase: created " << bubbleChainStringId(e0) << " and " <<
                        bubbleChainStringId(e1) << endl;
                }
            }

        } else {

            // Ambiguous. Don't detangle.
            if(debug) {
                cout << "Ambiguous. NOt detangling." << endl;
            }
            return false;
        }

    } else {



        // We are not using the Bayesian model.

        // Count the number of significant, ambiguous, and negligible elements
        // in the tangle matrix.
        uint64_t significantCount = 0;
        uint64_t ambiguousCount = 0;
        uint64_t negligibleCount = 0;
        for(uint64_t i0=0; i0<inEdges.size(); i0++) {
            for(uint64_t i1=0; i1<outEdges.size(); i1++) {
                const uint64_t t = tangleMatrix[i0][i1];
                if(t <= detangleToleranceLow) {
                    ++negligibleCount;
                } else if(t >= detangleToleranceHigh) {
                    ++significantCount;
                } else {
                    ++ambiguousCount;
                }
            }
        }

        // If the tangle matrix contains any ambiguous elements, do nothing.
        if(ambiguousCount > 0) {
            return false;
        }

        // There are no ambiguous elements.
        // If there are no negligible element, that is all elements of the tangle matrix are significant,
        // there is nothing to do.
        if(negligibleCount == 0) {
            return false;
        }

        // To avoid breaking contiguity, we require each column and each row of the
        // tangle matrix to have at least one significant element.
        // This means that each in-edge will be "merged" with at least one out-edge,
        // and each out-edge will be "merged" with at least one in-edge.
        // ACTUALY, FOR MORE ROBUSTNESS REQUIRE EXACTLY OEN SIGNIFICANT ELEMENT PER ROW AND COLUMN.
        for(uint64_t i0=0; i0<inEdges.size(); i0++) {
            uint64_t significantCount = 0;
            for(uint64_t i1=0; i1<outEdges.size(); i1++) {
                if(tangleMatrix[i0][i1] >= detangleToleranceHigh) {
                    ++significantCount;
                }
            }
            if(significantCount != 1) {
                return false;
            }
        }
        for(uint64_t i1=0; i1<outEdges.size(); i1++) {
            uint64_t significantCount = 0;
            for(uint64_t i0=0; i0<inEdges.size(); i0++) {
                if(tangleMatrix[i0][i1] >= detangleToleranceHigh) {
                    ++significantCount;
                }
            }
            if(significantCount != 1) {
                return false;
            }
        }

    #if 0
        // In an in-edge is also an out-edge, don't detangle.
        for(const edge_descriptor ce: inEdges) {
            if(find(outEdges.begin(), outEdges.end(), ce) != outEdges.end()) {
                if(debug) {
                    cout << "Not degangled because an in-edge is also an out-edge." << endl;
                }
                return false;
            }
        }
    #endif

        if(debug) {
            cout << "This edge will be detangled " << inEdges.size() << " by " << outEdges.size() << endl;
        }

        // Create truncated versions of the inEdges and outEdges.
        vector<vertex_descriptor> inVertices;
        for(const edge_descriptor ce: inEdges) {
            inVertices.push_back(cloneAndTruncateAtEnd(ce));
        }
        vector<vertex_descriptor> outVertices;
        for(const edge_descriptor ce: outEdges) {
            outVertices.push_back(cloneAndTruncateAtBeginning(ce));
        }


        // Each significant element of the tangle matrix generates a new edge.
        for(uint64_t i0=0; i0<inEdges.size(); i0++) {
            for(uint64_t i1=0; i1<outEdges.size(); i1++) {
                if(tangleMatrix[i0][i1] >= detangleToleranceHigh) {
                    const edge_descriptor ceNew = connect(inVertices[i0], outVertices[i1]);
                    if(debug) {
                        cout << "Created " << bubbleChainStringId(ceNew) << endl;
                    }
                }
            }
        }
    }


    // Now we can remove cv0, cv1, ce, and all of the in-edges and out-edges.
    // We have to do this while safely incrementing the edge iterator to point to the
    // next edge that was not removed.
    // We already incremented the iterator to point past ce.
    boost::remove_edge(ce, cGraph);
    for(const edge_descriptor ce: inEdges) {
        if(it != edgeMap.end() and cGraph[ce].id == it->first) {
            ++it;
        }
        edgeMap.erase(cGraph[ce].id);
        boost::remove_edge(ce, cGraph);
    }
    for(const edge_descriptor ce: outEdges) {
        if(it != edgeMap.end() and cGraph[ce].id == it->first) {
            ++it;
        }
        edgeMap.erase(cGraph[ce].id);
        boost::remove_edge(ce, cGraph);
    }
    for(const edge_descriptor ce: backEdges) {
        if(it != edgeMap.end() and cGraph[ce].id == it->first) {
            ++it;
        }
        edgeMap.erase(cGraph[ce].id);
        boost::remove_edge(ce, cGraph);
    }
    cGraph.removeVertex(cv0);
    cGraph.removeVertex(cv1);

    return true;
#endif
}
#endif


// Detangle short superbubbles with any number of entrances and exits.
bool AssemblyGraph::detangleShortSuperbubbles(
    bool debug,
    uint64_t maxOffset1,    // Used to define superbubbles
    uint64_t detangleToleranceLow,
    uint64_t detangleToleranceHigh)
{
    AssemblyGraph& cGraph = *this;

    // Find the superbubbles.
    Superbubbles superbubbles(cGraph, maxOffset1);

    // Loop over the superbubbles.
    bool changesWereMade = false;
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        if(detangleShortSuperbubble(debug,
            superbubbles, superbubbleId, detangleToleranceLow, detangleToleranceHigh)) {
            changesWereMade = true;
        }
    }

    return changesWereMade;
}



bool AssemblyGraph::detangleShortSuperbubble(
    bool debug,
    const Superbubbles& superbubbles,
    uint64_t superbubbleId,
    uint64_t detangleToleranceLow,
    uint64_t detangleToleranceHigh)
{
    AssemblyGraph& cGraph = *this;
    const Superbubble& superbubble = superbubbles.getSuperbubble(superbubbleId);

    if(debug) {
        cout << "Found a superbubble with " << superbubble.size() <<
            " vertices:";
        for(const vertex_descriptor cv: superbubble) {
            cout << " " << cGraph[cv].edgeId;
        }
        cout << endl;
    }

    // Fill in the in-edges and out-edges.
    // These cannot be computed while constructing the superbubbles
    // as they can change when other superbubbles are detangled.
    vector<edge_descriptor> inEdges;
    vector<edge_descriptor> outEdges;
    for(const vertex_descriptor cv0: superbubble) {
        BGL_FORALL_INEDGES(cv0, ce, cGraph, AssemblyGraph) {
            const vertex_descriptor cv1 = source(ce, cGraph);
            if(not superbubbles.isInSuperbubble(superbubbleId, cv1)) {
                inEdges.push_back(ce);
            }
        }
        BGL_FORALL_OUTEDGES(cv0, ce, cGraph, AssemblyGraph) {
            const vertex_descriptor cv1 = target(ce, cGraph);
            if(not superbubbles.isInSuperbubble(superbubbleId, cv1)) {
                outEdges.push_back(ce);
            }
        }
    }
    const uint64_t inDegree = inEdges.size();
    const uint64_t outDegree = outEdges.size();

    if(debug) {
        cout << inDegree << " in-edges:";
        for(const edge_descriptor ce: inEdges) {
            cout << " " << bubbleChainStringId(ce);
        }
        cout << endl;
        cout << outDegree << " out-edges:";
        for(const edge_descriptor ce: outEdges) {
            cout << " " << bubbleChainStringId(ce);
        }
        cout << endl;
    }

    if(inDegree == 0 or outDegree == 0) {
        if(debug) {
            cout << "Not detangling due to degree (case 1)." << endl;
        }
        return false;
    }

#if 0
    // Skip this check. We still want to remove the superbubble if possible.
    if(inDegree < 2 and outDegree < 2) {
        if(debug) {
            cout << "Not detangling due to degree (case 2)." << endl;
        }
        return false;
    }
#endif

    // This requires the last bubble of each in-edge
    // and the first bubble of each out-edge to be haploid.
    bool canDo = true;
    for(const edge_descriptor ce: inEdges) {
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.lastBubble().isHaploid()) {
            if(debug) {
                cout << "Not detangling because the last bubble of in-edge " <<
                    bubbleChainStringId(ce) << " is not haploid." << endl;
            }
            canDo = false;
            break;
        }
    }
    for(const edge_descriptor ce: outEdges) {
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.firstBubble().isHaploid()) {
            if(debug) {
                cout << "Not detangling because the first bubble of out-edge " <<
                    bubbleChainStringId(ce) << " is not haploid." << endl;
            }
            canDo = false;
            break;
        }
    }
    if(not canDo) {
        return false;
    }


    // Compute the tangle matrix.
    vector< vector<uint64_t> > tangleMatrix;
    computeTangleMatrix(inEdges, outEdges, tangleMatrix, true);

    if(debug) {
        cout << "Tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            const edge_descriptor inEdge = inEdges[i0];
            for(uint64_t i1=0; i1<outDegree; i1++) {
                const edge_descriptor outEdge = outEdges[i1];

                cout << bubbleChainStringId(inEdge) << " " <<
                    bubbleChainStringId(outEdge) << " " << tangleMatrix[i0][i1];

                cout << endl;
            }
        }
    }

    // Count the number of significant, ambiguous, and negligible elements
    // in the tangle matrix.
    uint64_t significantCount = 0;
    uint64_t ambiguousCount = 0;
    uint64_t negligibleCount = 0;
    for(uint64_t i0=0; i0<inDegree; i0++) {
        for(uint64_t i1=0; i1<outDegree; i1++) {
            const uint64_t t = tangleMatrix[i0][i1];
            if(t <= detangleToleranceLow) {
                ++negligibleCount;
            } else if(t >= detangleToleranceHigh) {
                ++significantCount;
            } else {
                ++ambiguousCount;
            }
        }
    }

    // If the tangle matrix contains any ambiguous elements, do nothing.
    if(ambiguousCount > 0) {
        if(debug) {
            cout << "Not detangled because the tangle matrix contains ambiguous elements." << endl;
        }
        return false;
    }

#if 0
    // (Skip this check - we still want to get rid of the superbubble in that case too!)
    // There are no ambiguous elements.
    // If there are no negligible element, that is all elements of the tangle matrix are significant,
    // there is nothing to do.
    if(negligibleCount == 0) {
        if(debug) {
            cout << "Not detangled because the tangle matrix contains no negligible elements." << endl;
        }
        return false;
    }
#endif

    // To avoid breaking contiguity, we require each column and each row of the
    // tangle matrix to have at least one significant element.
    // This means that each in-edge will be "merged" with at least one out-edge,
    // and each out-edge will be "merged" with at least one in-edge.
    bool ok = true;
    for(uint64_t i0=0; i0<inDegree; i0++) {
        bool foundSignificant = false;
        for(uint64_t i1=0; i1<outDegree; i1++) {
            if(tangleMatrix[i0][i1] >= detangleToleranceHigh) {
                foundSignificant = true;
                break;
            }
        }
        if(not foundSignificant) {
            ok = false;
            break;
        }
    }
    for(uint64_t i1=0; i1<outDegree; i1++) {
        bool foundSignificant = false;
        for(uint64_t i0=0; i0<inDegree; i0++) {
            if(tangleMatrix[i0][i1] >= detangleToleranceHigh) {
                foundSignificant = true;
                break;
            }
        }
        if(not foundSignificant) {
            ok = false;
            break;
        }
    }
    if(not ok) {
        if(debug) {
            cout << "Not detangled to avoid breaking contiguity." << endl;
        }
        return false;
    }

    if(debug) {
        cout << "This superbubble will be detangled." << endl;
    }

    // Create truncated versions of the inEdges and outEdges.
    vector<vertex_descriptor> inVertices;
    for(const edge_descriptor ce: inEdges) {
        inVertices.push_back(cloneAndTruncateAtEnd(ce));
    }
    vector<vertex_descriptor> outVertices;
    for(const edge_descriptor ce: outEdges) {
        outVertices.push_back(cloneAndTruncateAtBeginning(ce));
    }



    // Each significant element of the tangle matrix generates a new edge.
    for(uint64_t i0=0; i0<inEdges.size(); i0++) {
        for(uint64_t i1=0; i1<outEdges.size(); i1++) {
            if(tangleMatrix[i0][i1] >= detangleToleranceHigh) {
                connect(inVertices[i0], outVertices[i1]);
            }
        }
    }
    if(debug) {
        cout << "After creating new edges, nextEdgeId is " << nextEdgeId << endl;
    }


#if 0
    // Each significant element of the tangle matrix generates a new edge,
    // obtained by "merging" an in-edge with an out-edge.
    for(uint64_t i0=0; i0<inEdges.size(); i0++) {
        const edge_descriptor ce0 = inEdges[i0];
        const BubbleChain& bubbleChain0 = cGraph[ce0];
        const Bubble& bubble0 = bubbleChain0.lastBubble();
        SHASTA_ASSERT(bubble0.isHaploid());
        const Chain& chain0 = bubble0.front();
        SHASTA_ASSERT(chain0.size() >= 2);
        for(uint64_t i1=0; i1<outEdges.size(); i1++) {
            if(tangleMatrix[i0][i1] < detangleToleranceHigh) {
                continue;
            }
            const edge_descriptor ce1 = outEdges[i1];
            const BubbleChain& bubbleChain1 = cGraph[ce1];
            const Bubble& bubble1 = bubbleChain1.firstBubble();
            SHASTA_ASSERT(bubble1.isHaploid());
            const Chain& chain1 = bubble1.front();
            SHASTA_ASSERT(chain1.size() >= 2);

            edge_descriptor eNew;
            tie(eNew, ignore) = add_edge(source(ce0, cGraph), target(ce1, cGraph), cGraph);
            AssemblyGraphEdge& newEdge = cGraph[eNew];
            newEdge.id = nextEdgeId++;
            BubbleChain& newBubbleChain = newEdge;

            if(debug) {
                cout << "Merging " <<
                    bubbleChainStringId(ce0) << " " <<
                    bubbleChainStringId(ce1) << " into " <<
                    bubbleChainStringId(eNew) << endl;
            }

            // Create the new BubbleChain. It is obtained by joining
            // bubbleChain0 and bubbleChain1, with vertex cv
            // removed from the end of bubbleChain0
            // and from the beginning of bubbleChain1.
            // Here we use the above assumption that
            // the last bubble of bubbleChain0 and the first bubble of bubbleChain1
            // are haploid.
            newBubbleChain = bubbleChain0;

            // Remove the last marker graph edge, which is in the superbubble.
            Bubble& newBubbleLast = newBubbleChain.back();
            SHASTA_ASSERT(newBubbleLast.size() == 1);
            Chain& newChainLast = newBubbleLast.front();
            newChainLast.resize(newChainLast.size() - 1);

            // Append chain1, except for the first marker graph edge, which is in the superbubble.
            copy(chain1.begin() + 1, chain1.end(), back_inserter(newChainLast));

            // Append the rest of bubbleChain1.
            copy(bubbleChain1.begin() + 1, bubbleChain1.end(), back_inserter(newBubbleChain));
        }

    }
#endif

    // Now we can remove all the vertices in the superbubble.
    for(const vertex_descriptor cv: superbubble) {
        clear_vertex(cv, cGraph);
        remove_vertex(cv, cGraph);
    }

    return true;
}



// The above versions require the last bubble of superbubble in-edges
// and the first bubble of superbubble out-edges to be haploid.
// This version does not.
bool AssemblyGraph::detangleShortSuperbubblesGeneral(
    bool debug,
    uint64_t maxOffset1,    // Used to define superbubbles
    uint64_t detangleToleranceLow,
    uint64_t detangleToleranceHigh)
{
    AssemblyGraph& cGraph = *this;

    // Find the superbubbles.
    Superbubbles superbubbles(cGraph, maxOffset1);

    // Loop over the superbubbles.
    bool changesWereMade = false;
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        if(detangleShortSuperbubbleGeneral(debug,
            superbubbles, superbubbleId, detangleToleranceLow, detangleToleranceHigh)) {
            changesWereMade = true;
        }
    }

    return changesWereMade;

}



bool AssemblyGraph::detangleShortSuperbubbleGeneral(
    bool debug,
    const Superbubbles& superbubbles,
    uint64_t superbubbleId,
    uint64_t detangleToleranceLow,
    uint64_t detangleToleranceHigh)
{
    AssemblyGraph& cGraph = *this;
    const Superbubble& superbubble = superbubbles.getSuperbubble(superbubbleId);

    if(debug) {
        cout << "General detangling of a superbubble with " << superbubble.size() <<
            " vertices:";
        for(const vertex_descriptor cv: superbubble) {
            cout << " " << cGraph[cv].edgeId;
        }
        cout << endl;
        cout << "nextEdgeId is " << nextEdgeId << endl;
    }

    // Gather the in-edges and out-edges.
    // These cannot be computed while constructing the superbubbles
    // as they can change when other superbubbles are detangled.
    vector<edge_descriptor> inEdges;
    vector<edge_descriptor> outEdges;
    for(const vertex_descriptor cv0: superbubble) {
        BGL_FORALL_INEDGES(cv0, ce, cGraph, AssemblyGraph) {
            const vertex_descriptor cv1 = source(ce, cGraph);
            if(not superbubbles.isInSuperbubble(superbubbleId, cv1)) {
                inEdges.push_back(ce);
            }
        }
        BGL_FORALL_OUTEDGES(cv0, ce, cGraph, AssemblyGraph) {
            const vertex_descriptor cv1 = target(ce, cGraph);
            if(not superbubbles.isInSuperbubble(superbubbleId, cv1)) {
                outEdges.push_back(ce);
            }
        }
    }

    if(debug) {
        cout << inEdges.size() << " in-edges:";
        for(const edge_descriptor ce: inEdges) {
            cout << " " << bubbleChainStringId(ce);
        }
        cout << endl;
        cout << outEdges.size() << " out-edges:";
        for(const edge_descriptor ce: outEdges) {
            cout << " " << bubbleChainStringId(ce);
        }
        cout << endl;
    }



    // See if we can use detangleShortSuperbubble.
    // This is the case if:
    // - The last bubble of every in-edge is haploid.
    // - The first bubble of every out-edge is haploid.
    bool isSimpleCase = true;
    for(const edge_descriptor ce: inEdges) {
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.lastBubble().isHaploid()) {
            isSimpleCase = false;
        }
    }
    for(const edge_descriptor ce: outEdges) {
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.firstBubble().isHaploid()) {
            isSimpleCase = false;
        }
    }
    if(isSimpleCase) {
        if(debug) {
            cout << "Simple case: calling detangleShortSuperbubble." << endl;
        }
        return detangleShortSuperbubble(debug, superbubbles, superbubbleId,
            detangleToleranceLow, detangleToleranceHigh);
    }
    if(debug) {
        cout << "General case." << endl;
    }


    // We will compute a generalized tangle matrix that takes into
    // account all incoming and outgoing chains.

    // Gather the second to last marker graph edge of each chain of
    // the last bubble of each in-edge.
    vector<MarkerGraphEdgeId> in;
    for(const edge_descriptor ce: inEdges) {
        const BubbleChain& bubbleChain = cGraph[ce];
        const Bubble& lastBubble = bubbleChain.lastBubble();
        for(const Chain& chain: lastBubble) {
            in.push_back(chain.secondToLast());
        }
    }

    // Gather the second marker graph edge of each chain of
    // the first bubble of each out-edge.
    vector<MarkerGraphEdgeId> out;
    for(const edge_descriptor ce: outEdges) {
        const BubbleChain& bubbleChain = cGraph[ce];
        const Bubble& firstBubble = bubbleChain.firstBubble();
        for(const Chain& chain: firstBubble) {
            out.push_back(chain.second());
        }
    }

    if(debug) {
        cout << in.size() << " incoming marker graph edges:";
        for(const MarkerGraphEdgeId edgeId: in) {
            cout << " " << edgeId;
        }
        cout << endl;
        cout << out.size() << " outgoing marker graph edges:";
        for(const MarkerGraphEdgeId edgeId: out) {
            cout << " " << edgeId;
        }
        cout << endl;
    }

    // Compute the tangle matrix.
    vector< vector<uint64_t> > tangleMatrix(in.size(), vector<uint64_t>(out.size()));
    for(uint64_t i0=0; i0<in.size(); i0++) {
        const MarkerGraphEdgeId edgeId0 = in[i0];
        for(uint64_t i1=0; i1<out.size(); i1++) {
            const MarkerGraphEdgeId edgeId1 = out[i1];

            if(edgeId1 == assembler.markerGraph.reverseComplementEdge[edgeId0]) {
                tangleMatrix[i0][i1] = 0;
            } else {
                MarkerGraphEdgePairInfo info;
                SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, info));
                tangleMatrix[i0][i1] = info.common;
            }
        }
    }

    if(debug) {
        cout << "Tangle matrix:" << endl;
        for(uint64_t i0=0; i0<in.size(); i0++) {
            const MarkerGraphEdgeId edgeId0 = in[i0];
            for(uint64_t i1=0; i1<out.size(); i1++) {
                const MarkerGraphEdgeId edgeId1 = out[i1];
                cout << edgeId0 << " " << edgeId1 << " " << tangleMatrix[i0][i1] << endl;
            }
        }
    }

    // Count the number of significant, ambiguous, and negligible elements
    // in the tangle matrix.
    uint64_t significantCount = 0;
    uint64_t ambiguousCount = 0;
    uint64_t negligibleCount = 0;
    for(uint64_t i0=0; i0<in.size(); i0++) {
        for(uint64_t i1=0; i1<out.size(); i1++) {
            const uint64_t t = tangleMatrix[i0][i1];
            if(t <= detangleToleranceLow) {
                ++negligibleCount;
            } else if(t >= detangleToleranceHigh) {
                ++significantCount;
            } else {
                ++ambiguousCount;
            }
        }
    }

    // If the tangle matrix contains any ambiguous elements, do nothing.
    if(ambiguousCount > 0) {
        if(debug) {
            cout << "Not detangled because the tangle matrix contains ambiguous elements." << endl;
        }
        return false;
    }

#if 0
    // (Skip this check - we still want to get rid of the superbubble in that case too!)
    // There are no ambiguous elements.
    // If there are no negligible element, that is all elements of the tangle matrix are significant,
    // there is nothing to do.
    if(negligibleCount == 0) {
        if(debug) {
            cout << "Not detangled because the tangle matrix contains no negligible elements." << endl;
        }
        return false;
    }
#endif

    // To avoid breaking contiguity, we require each column and each row of the
    // tangle matrix to have at least one significant element.
    // This means that each in-edge will be "merged" with at least one out-edge,
    // and each out-edge will be "merged" with at least one in-edge.
    bool ok = true;
    for(uint64_t i0=0; i0<in.size(); i0++) {
        bool foundSignificant = false;
        for(uint64_t i1=0; i1<out.size(); i1++) {
            if(tangleMatrix[i0][i1] >= detangleToleranceHigh) {
                foundSignificant = true;
                break;
            }
        }
        if(not foundSignificant) {
            ok = false;
            break;
        }
    }
    for(uint64_t i1=0; i1<out.size(); i1++) {
        bool foundSignificant = false;
        for(uint64_t i0=0; i0<in.size(); i0++) {
            if(tangleMatrix[i0][i1] >= detangleToleranceHigh) {
                foundSignificant = true;
                break;
            }
        }
        if(not foundSignificant) {
            ok = false;
            break;
        }
    }
    if(not ok) {
        if(debug) {
            cout << "Not detangled to avoid breaking contiguity." << endl;
        }
        return false;
    }

    if(debug) {
        cout << "This superbubble will be detangled but splitting of some edges is required." << endl;
    }

    // Make sure the last bubble of every in-edge is haploid.
    for(const edge_descriptor ce: inEdges) {
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.lastBubble().isHaploid()) {
            if(debug) {
                cout << "In-edge " << bubbleChainStringId(ce) <<
                    " needs to be split at the end." << endl;
            }
            splitBubbleChainAtEnd(ce);
        }
    }

    // Make sure the first bubble of every out-edge is haploid.
    for(const edge_descriptor ce: outEdges) {
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.firstBubble().isHaploid()) {
            if(debug) {
                cout << "Out-edge " << bubbleChainStringId(ce) <<
                    " needs to be split at the beginning." << endl;
            }
            splitBubbleChainAtBeginning(ce);
        }
    }
    if(debug) {
        cout << "After splitting of edges, nextEdgeId is " << nextEdgeId << endl;
    }

    // Now we can call detangleShortSuperbubble.
    if(debug) {
        cout << "Calling detangleShortSuperbubble after splitting some edges." << endl;
    }
    return detangleShortSuperbubble(debug, superbubbles, superbubbleId,
        detangleToleranceLow, detangleToleranceHigh);
}



// Special treatment to detangle back edges that were too long
// to be handled by detangleEdges.
bool AssemblyGraph::detangleBackEdges(
    uint64_t detangleToleranceLow,
    uint64_t detangleToleranceHigh)
{
    cout << "Detangling back edges." << endl;
    AssemblyGraph& cGraph = *this;

    // To safely iterate over edges while removing edges we must use edge ids
    // as unique identifiers, because edge descriptors can be reused as edges are
    // deleted a new edges are created.
    std::map<uint64_t, edge_descriptor> edgeMap;
    BGL_FORALL_EDGES(ce, cGraph, AssemblyGraph) {
        edgeMap.insert({cGraph[ce].id, ce});
    }

    uint64_t detangleCount = 0;;
    for(auto it=edgeMap.begin(); it!=edgeMap.end(); /* Incremented safely by detangleEdgeStrict */) {
        if(detangleBackEdge(edgeMap, it, detangleToleranceLow, detangleToleranceHigh)) {
            ++detangleCount;
        }
    }
    cout << "Detangled " << detangleCount << " back edges." << endl;

    return detangleCount > 0;

}



// Special treatment to detangle back edges that were too long
// to be handled by detangleEdge.
bool AssemblyGraph::detangleBackEdge(
    std::map<uint64_t, edge_descriptor>& edgeMap,
    std::map<uint64_t, edge_descriptor>::iterator& it,
    uint64_t detangleToleranceLow,
    uint64_t detangleToleranceHigh)
{
    AssemblyGraph& cGraph = *this;
    const edge_descriptor ce = it->second;
    ++it;
    // edgeMap.erase(cGraph[ce].id);

    const bool debug = false;

    // Tangle matrix elements <= detangleToleranceLow are treated as negligible.
    // Tangle matrix elements >= detangleToleranceHigh are treated as significant.
    // Tangle matrix elements in between are considered ambiguous.
    SHASTA_ASSERT(detangleToleranceHigh > detangleToleranceLow);

    const vertex_descriptor cv0 = source(ce, cGraph);
    const vertex_descriptor cv1 = target(ce, cGraph);

    // Check the degrees.
    if(out_degree(cv0, cGraph) != 1) {
        return false;
    }
    if(in_degree(cv1, cGraph) != 1) {
        return false;
    }

    // Look for a back edge.
    vector<edge_descriptor> backEdges;
    BGL_FORALL_OUTEDGES(cv1, ce, cGraph, AssemblyGraph) {
        if(target(ce, cGraph) == cv0) {
            backEdges.push_back(ce);
        }
    }
    if(backEdges.empty()) {
        return false;
    }

    // Only attempt to handle the case with a single back-edge.
    if(backEdges.size() != 1) {
        return false;
    }
    const edge_descriptor ceBack = backEdges.front();

    if(debug) {
        cout << "Attempting to detangle edge " << bubbleChainStringId(ce) <<
            " with back-edge " << bubbleChainStringId(ceBack) << endl;
    }

    // The back-edge is both an in-edge and an out-edge.
    // Store it at the first position of both inEdges and outEdges.

    // Gather the in-edges.
    vector<edge_descriptor> inEdges(1, ceBack);
    BGL_FORALL_INEDGES(cv0, ce, cGraph, AssemblyGraph) {
        if(ce == ceBack) {
            continue;
        }
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.lastBubble().isHaploid()) {
            if(debug) {
                cout << "Not detangling because the last bubble of in-edge " <<
                    bubbleChainStringId(ce) << " is not haploid." << endl;
            }
            return false;
        }
        inEdges.push_back(ce);
    }

    // Gather the out-edges.
    vector<edge_descriptor> outEdges(1, ceBack);
    BGL_FORALL_OUTEDGES(cv1, ce, cGraph, AssemblyGraph) {
        if(ce == ceBack) {
            continue;
        }
        const BubbleChain& bubbleChain = cGraph[ce];
        if(not bubbleChain.firstBubble().isHaploid()) {
            if(debug) {
                cout << "Not detangling because the first bubble of out-edge " <<
                    bubbleChainStringId(ce) << " is not haploid." << endl;
            }
            return false;
        }
        outEdges.push_back(ce);
    }


    if(debug) {

        // Position 0 of the inEdges and outEdges stores the back-edge.

        cout << "In-edges: ";
        for(uint64_t i=1; i<inEdges.size(); i++) {
            const edge_descriptor ce = inEdges[i];
            cout << " " << bubbleChainStringId(ce);
        }
        cout << endl;

        cout << "Out-edges: ";
        for(uint64_t i=1; i<outEdges.size(); i++) {
            const edge_descriptor ce = outEdges[i];
            cout << " " << bubbleChainStringId(ce);
        }
        cout << endl;
    }
    // Compute the tangle matrix.
    vector< vector<uint64_t> > tangleMatrix;
    computeTangleMatrix(inEdges, outEdges, tangleMatrix, false);

    if(debug) {
        cout << "Tangle matrix:" << endl;
        for(uint64_t i0=0; i0<inEdges.size(); i0++) {
            const edge_descriptor ce0 = inEdges[i0];
            for(uint64_t i1=0; i1<outEdges.size(); i1++) {
                const edge_descriptor ce1 = outEdges[i1];
                cout <<
                    bubbleChainStringId(ce0) << " " <<
                    bubbleChainStringId(ce1) << " " <<
                    tangleMatrix[i0][i1];
                cout << endl;
            }
        }
    }

    return false;
}



void AssemblyGraph::phaseBubbleChainsUsingPhasingGraph(
    bool debug,
    uint64_t n, // Maximum number of Chain MarkerGraphEdgeIds to use when computing tangle matrices.
    uint64_t lowThreshold,
    uint64_t highThreshold,
    bool useBayesianModel,
    double epsilon,
    double minLogP,
    uint64_t longBubbleThreshold)
{
    AssemblyGraph& cGraph = *this;

    if(debug) {
        cout << "phaseBubbleChainsUsingPhasingGraph begins." << endl;
    }

    vector<edge_descriptor> allEdges;
    BGL_FORALL_EDGES(ce, cGraph, AssemblyGraph) {
        allEdges.push_back(ce);
    }

    for(const edge_descriptor ce: allEdges) {
        phaseBubbleChainUsingPhasingGraph(ce, n, lowThreshold, highThreshold, useBayesianModel, epsilon, minLogP, longBubbleThreshold, debug);
    }

    if(debug) {
        cout << "phaseBubbleChainsUsingPhasingGraph ends." << endl;
    }
}



void AssemblyGraph::phaseBubbleChainsUsingPhasingTable(
    const string& debugOutputFileNamePrefix,
    double phaseErrorThreshold,
    double bubbleErrorThreshold,
    uint64_t longBubbleThreshold)
{
    AssemblyGraph& cGraph = *this;

    const bool debug = not debugOutputFileNamePrefix.empty();
    if(debug) {
        cout << "phaseBubbleChainsUsingPhasingTable begins." << endl;
    }
    performanceLog << timestamp << "AssemblyGraph::phaseBubbleChainsUsingPhasingTable begins." << endl;

    // If debug output was requested, make sure we have a directory
    // where the debug output files will go.
    string directoryName;
    if(debug) {
        directoryName = debugOutputFileNamePrefix + "-PhasingTables";
        std::filesystem::create_directory(directoryName);
    }

    vector<edge_descriptor> allEdges;
    BGL_FORALL_EDGES(ce, cGraph, AssemblyGraph) {
        allEdges.push_back(ce);
    }

    for(const edge_descriptor ce: allEdges) {
        phaseBubbleChainUsingPhasingTable(
            debug ? (directoryName + "/" + bubbleChainStringId(ce)) : "",
            ce, phaseErrorThreshold, bubbleErrorThreshold, longBubbleThreshold);
    }

    if(debug) {
        cout << "phaseBubbleChainsUsingPhasingTable ends." << endl;
    }
    performanceLog << timestamp << "AssemblyGraph::phaseBubbleChainsUsingPhasingTable ends." << endl;

}



void AssemblyGraph::phaseBubbleChainUsingPhasingGraph(
    edge_descriptor ce,
    uint64_t n, // Maximum number of Chain MarkerGraphEdgeIds to use when computing tangle matrices.
    uint64_t lowThreshold,
    uint64_t highThreshold,
    bool useBayesianModel,
    double epsilon,
    double minLogP,
    uint64_t longBubbleThreshold,
    bool debug)
{
    AssemblyGraph& cGraph = *this;
    BubbleChain& bubbleChain = cGraph[ce];

    // debug = debug and (cGraph[ce].id == 500048);

    if(debug) {
        cout << "Phasing " << bubbleChainStringId(ce) << endl;
    }

    const bool detailedDebug = debug; // (cGraph[ce].id == 49557);

    // If this bubble chain has a single bubble, there is nothing to do.
    if(bubbleChain.size() == 1) {
        if(debug) {
            cout << "Not phased because it has only one bubble." << endl;
        }
        return;
    }

    // Table to contain the Phasing graph vertex corresponding to each diploid bubble.
    // Indexed by the bubble position in the bubble chains, and contains
    // PhasingGraph::null_vertex() for non-diploid bubbles.
    vector<PhasingGraph::vertex_descriptor> vertexTable(bubbleChain.size(), PhasingGraph::null_vertex());

    // Create the PhasingGraph and its vertices, one for
    // each diploid bubble in the bubble chain.
    PhasingGraph phasingGraph;
    for(uint64_t i=0; i<bubbleChain.size(); i++) {
        if(bubbleChain[i].isDiploid()) {
            vertexTable[i] = add_vertex({i, 0}, phasingGraph);
        }
    }

    // Write a histogram of the bubbles in this bubble chain by ploidy.
    if(debug) {
        cout << "Phasing a bubble chain with " << bubbleChain.size() << " bubbles." << endl;
        vector<uint64_t> histogram;
        for(const Bubble& bubble: bubbleChain) {
            const uint64_t ploidy = bubble.size();
            if(histogram.size() <= ploidy) {
                histogram.resize(ploidy + 1);
            }
            ++histogram[ploidy];
        }
        for(uint64_t ploidy=1; ploidy<histogram.size(); ploidy++) {
            const uint64_t frequency = histogram[ploidy];
            if(frequency) {
                cout << frequency << " bubbles of ploidy " << ploidy << endl;
            }
        }
    }

#if 0
    // If this bubble chain has less than two diploid bubbles, there is nothing to do.
    uint64_t diploidBubblesCount = 0;
    for(const Bubble& bubble: bubbleChain) {
        if(bubble.size() == 2) {
            ++diploidBubblesCount;
        }
    }
    if(diploidBubblesCount < 2) {
        if(debug) {
            cout << "Not phased because it has less than 2 diploid bubbles." << endl;
        }
        return;
    }
#endif

    // Add edges of the phasing graph.
    for(uint64_t i0=0; i0<bubbleChain.size()-1; i0++) {
        const PhasingGraph::vertex_descriptor pv0 = vertexTable[i0];
        if(pv0 == PhasingGraph::null_vertex()) {
            continue;
        }

        // Gather the next-to-last two marker graph edges for the two chains
        // of this bubble.
        const Bubble& bubble0 = bubbleChain[i0];
        SHASTA_ASSERT(bubble0.size() == 2);
        const Chain& chain00 = bubble0[0];
        const Chain& chain01 = bubble0[1];
        const array<MarkerGraphEdgeId, 2> edges0 =
            {chain00[chain00.size()-2], chain01[chain01.size()-2]};

        for(uint64_t i1=i0+1; i1<bubbleChain.size(); i1++) {
            const PhasingGraph::vertex_descriptor pv1 = vertexTable[i1];
            if(pv1 == PhasingGraph::null_vertex()) {
                continue;
            }

            // Gather the next-to-last two marker graph edges for the two chains
            // of this bubble.
            const Bubble& bubble1 = bubbleChain[i1];
            SHASTA_ASSERT(bubble1.size() == 2);
            const Chain& chain10 = bubble1[0];
            const Chain& chain11 = bubble1[1];
            const array<MarkerGraphEdgeId, 2> edges1 =
                {chain10[1], chain11[1]};

            // Compute the tangle matrix.
            TangleMatrix tangleMatrix;
            if(n == 1) {
                for(uint64_t j0=0; j0<2; j0++) {
                    for(uint64_t j1=0; j1<2; j1++) {
                        MarkerGraphEdgePairInfo info;
                        SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(
                            edges0[j0], edges1[j1], info));
                        tangleMatrix[j0][j1] = info.common;
                    }
                }
            } else {
                computeTangleMatrix(
                    {&chain00, &chain01},
                    {&chain10, &chain11},
                    n, tangleMatrix);
            }

            // Analyze the tangle matrix.
            int64_t phase;
            uint64_t minConcordant;
            uint64_t maxDiscordant;
            uint64_t total;
            double logPInPhase;
            double logPOutOfPhase;
            tangleMatrix.analyze(
                lowThreshold,
                highThreshold,
                phase,
                minConcordant,
                maxDiscordant,
                total,
                epsilon,
                logPInPhase,
                logPOutOfPhase);

            // If no common reads, stop the loop on i1.
            if(total == 0) {
                break;
            }

            if(detailedDebug) {
                cout << "Tangle matrix " << i0 << " " << i1 << ": " <<
                    tangleMatrix[0][0] << " " <<
                    tangleMatrix[0][1] << " " <<
                    tangleMatrix[1][0] << " " <<
                    tangleMatrix[1][1] << endl;
                cout << "minConcordant " << minConcordant << endl;
                cout << "maxDiscordant " << maxDiscordant << endl;
                cout << "log[p(in-phase)/p(random)] = " << logPInPhase <<
                    " dB, log[p(out-of-phase)/p(random)] = " << logPOutOfPhase << " dB." << endl;
            }

            // If using the Bayesian model, redefine the phase based on logPInPhase and logPOutOfPhase.
            if(useBayesianModel) {
                if((logPInPhase > minLogP) and (logPInPhase - logPOutOfPhase) > minLogP) {
                    phase = +1;
                } else  if((logPOutOfPhase > minLogP) and (logPOutOfPhase - logPInPhase) > minLogP) {
                    phase = -1;
                } else {
                    phase = 0;
                }
            }

            // If not ambiguous, add an edge to the PhasingGraph.
            if(phase != 0) {
                boost::add_edge(pv0, pv1, {phase, minConcordant, maxDiscordant, logPInPhase, logPOutOfPhase}, phasingGraph);

                if(detailedDebug) {
                    cout << " Added phasing graph edge " <<
                        phasingGraph[pv0].positionInBubbleChain << " " <<
                        phasingGraph[pv1].positionInBubbleChain << " with minConcordant " <<
                        minConcordant << ", maxDiscordant " << maxDiscordant << endl;
                }
            } else {
                if(detailedDebug) {
                    cout << " No phasing graph edge for " <<
                        phasingGraph[pv0].positionInBubbleChain << " " <<
                        phasingGraph[pv1].positionInBubbleChain << endl;
                }
            }

        }
    }

    if(debug) {
        const uint64_t vertexCount = num_vertices(phasingGraph);
        const uint64_t edgeCount = num_edges(phasingGraph);
        const double connectivity = 2. * double(edgeCount) / double(vertexCount);
        cout << "The phasing graph has " << vertexCount <<
            " vertices and " << edgeCount << " edges."
            " Average connectivity " << connectivity << endl;
    }

    phasingGraph.phase1(false, useBayesianModel);



    // Use the PhasedComponents in the PhasingGraph to create
    // a new BubbleChain that will replace the existing one.
    phaseBubbleChainUsingPhasedComponents(debug, ce, phasingGraph.phasedComponents, longBubbleThreshold);
}



// Use PhasedComponents to create a new BubbleChain that will replace the existing one.
void AssemblyGraph::phaseBubbleChainUsingPhasedComponents(
    bool debug,
    edge_descriptor e,
    const vector<shared_ptr<PhasedComponent> >& phasedComponents,
    uint64_t longBubbleThreshold)
{
    AssemblyGraph& cGraph = *this;
    BubbleChain& bubbleChain = cGraph[e];

    BubbleChain newBubbleChain;
    if(debug) {
        cout << "Creating the new bubble chain for " << bubbleChainStringId(e) << endl;
    }

    // Loop over the phased components.
    for(uint64_t i=0; /* Check later */; i++) {

        // Bubbles in-between phased components, or before the first phased component,
        // or after the last phased component.
        {
            const uint64_t beginPositionInBubbleChain =
                (i == 0) ? 0 : phasedComponents[i-1]->maxPositionInBubbleChain + 1;
            const uint64_t endPositionInBubbleChain =
                (i == phasedComponents.size()) ?
                bubbleChain.size() :
                phasedComponents[i]->minPositionInBubbleChain;


            if(debug) {
                cout << "Adding unphased bubbles at positions [" <<
                    beginPositionInBubbleChain << "," << endPositionInBubbleChain << ")" << endl;
            }

            for(uint64_t i=beginPositionInBubbleChain; i<endPositionInBubbleChain; i++) {
                const Bubble& bubble = bubbleChain[i];

                // This unphased bubble will be copied verbatim to the new chain if it is
                // haploid or if it is long.
                bool copyVerbatim = bubble.isHaploid();
                if(not copyVerbatim) {
                    uint64_t averageOffset;
                    uint64_t minOffset;
                    uint64_t maxOffset;
#if 0
                    if(bubbleOffsetNoException(bubble, averageOffset, minOffset, maxOffset)) {
                        copyVerbatim = maxOffset >= longBubbleThreshold;
                    } else {
                        copyVerbatim = false;
                    }
#else
                    bubbleOffset(bubble, averageOffset, minOffset, maxOffset);
                    copyVerbatim = maxOffset >= longBubbleThreshold;
#endif
                }

                if(copyVerbatim) {
                    newBubbleChain.push_back(bubble);
                } else {
                    // Just add a simple haploid bubble with only the source
                    // and target MarkerGraphEdgeIds.
                    Bubble newBubble;
                    newBubble.resize(1);    // Make it haploid
                    Chain& newChain = newBubble.front();    // Its only chain.
                    newChain.push_back(bubble.front().front()); // Source MarkerGraphEdgeId
                    newChain.push_back(bubble.front().back());  // Target MarkerGraphEdgeId
                    newBubbleChain.push_back(newBubble);
                }
            }
        }



        // If we are past the last phased component, we are done.
        if(i == phasedComponents.size()) {
            break;
        }

        // Add a diploid bubble for the i-th phased component.
        const PhasedComponent& phasedComponent = *phasedComponents[i];
        const uint64_t minPositionInBubbleChain = phasedComponent.minPositionInBubbleChain;
        const uint64_t maxPositionInBubbleChain = phasedComponent.maxPositionInBubbleChain;
        if(debug) {
            cout << "Adding phased bubbles at positions " <<
                minPositionInBubbleChain << "-" << maxPositionInBubbleChain << endl;
        }
        newBubbleChain.emplace_back();
        Bubble& newBubble = newBubbleChain.back();
        newBubble.resize(2);    // Make it diploid.
        Chain& newChain0 = newBubble[0];    // The first haplotype after phasing.
        Chain& newChain1 = newBubble[1];    // The second haplotype after phasing.

        // Add the source MarkerGraphEdgeId.
        newChain0.push_back(bubbleChain[minPositionInBubbleChain].front().front());
        newChain1.push_back(bubbleChain[minPositionInBubbleChain].front().front());

        // Add the internal MarkerGraphEdgeIds of all phased diploid bubbles in this PhasedComponent.
        for(const auto& p: phasedComponent) {
            const uint64_t positionInBubbleChain = p.first;
            const int64_t phase = p.second;
            SHASTA_ASSERT(phase==1 or phase==-1);
            const Bubble& bubble = bubbleChain[positionInBubbleChain];
            SHASTA_ASSERT(bubble.isDiploid());
            const Chain& chain0 = (phase==1) ? bubble[0] : bubble[1];
            const Chain& chain1 = (phase==1) ? bubble[1] : bubble[0];
            copy(chain0.begin()+1, chain0.end()-1, back_inserter(newChain0));
            copy(chain1.begin()+1, chain1.end()-1, back_inserter(newChain1));
        }

        // Add the target MarkerGraphEdgeId.
        newChain0.push_back(bubbleChain[maxPositionInBubbleChain].front().back());
        newChain1.push_back(bubbleChain[maxPositionInBubbleChain].front().back());
    }

    // Replace the old BubbleChain with the new one, leaving the id of the edge unchanged.
    newBubbleChain.compress();
    bubbleChain = newBubbleChain;
}



void AssemblyGraph::phaseBubbleChainUsingPhasingTable(
    const string& debugOutputFileNamePrefix,
    edge_descriptor e,
    double phaseErrorThreshold,
    double bubbleErrorThreshold,
    uint64_t longBubbleThreshold)
{
    AssemblyGraph& cGraph = *this;
    BubbleChain& bubbleChain = cGraph[e];

    const bool debug = not debugOutputFileNamePrefix.empty();

    cleanupBubbleChainUsingPhasingTable(
        debug ? (debugOutputFileNamePrefix + "-PreCleanup") : "",
        e,
        phaseErrorThreshold,
        bubbleErrorThreshold,
        longBubbleThreshold);


#if 0
    // If this bubble chain has a single bubble, there is nothing to do.
    // NOT TRUE, WE STILL MAY HAVE TO REMOVE SOME BUBBLES.
    if(bubbleChain.size() == 1) {
        if(debug) {
            cout << "Skipped because it has only one bubble." << endl;
        }
        return;
    }
#endif

    // Create the phasing table for this bubble chain.
    PhasingTable phasingTable(bubbleChain, assembler.markerGraph, phaseErrorThreshold);

    if(phasingTable.empty()) {
        if(debug) {
            cout << "Not phasing because the phasing table is empty." << endl;
        }
        return;
    }
#if 0
    // WE STILL MAY HAVE TO REMOVE SOME BUBBLES.
    if(phasingTable.bubbleCount() < 2) {
        if(debug) {
            cout << "Not phasing because the phasing table has less than 2 bubbles." << endl;
        }
        return;
    }
#endif

    if(debug) {
        const uint64_t totalCount = phasingTable.entryCount();
        const uint64_t ambiguousCount = phasingTable.ambiguousEntryCount();
        const uint64_t unambiguousCount = totalCount - ambiguousCount;
        const uint64_t bubbleCount = phasingTable.bubbleCount();
        const uint64_t orientedReadCount = phasingTable.orientedReadCount();
        const double coverage = double(unambiguousCount) / double(bubbleCount);

        cout << "Phasing table summary for " << bubbleChainStringId(e) << ":" << endl;
        cout << bubbleCount << " diploid bubbles." << endl;
        cout << orientedReadCount << " oriented reads." << endl;
        cout << unambiguousCount << " unambiguous entries." << endl;
        cout << ambiguousCount << " ambiguous entries." << endl;
        cout << "Average coverage " << std::round(coverage) << endl;
        cout << "Average number of diploid bubbles seen by each oriented read " <<
            std::round(double(unambiguousCount)/double(orientedReadCount)) << endl;
    }

    // Phasing of the phasing table.
    phasingTable.greedyPhasing();
    if(debug) {
        uint64_t consistentCount;
        uint64_t inconsistentCount;
        tie(consistentCount, inconsistentCount) = phasingTable.countConsistentEntries();

        cout << "After greedy phasing, the phasing table has " << consistentCount <<
            " consistent entries and " << inconsistentCount <<
            " inconsistent entries (" << consistentCount + inconsistentCount <<
            " total)." << endl;

        phasingTable.writePng(debugOutputFileNamePrefix + "-Consistency.png",
            PhasingTable::ColoringMethod::byConsistency);
        phasingTable.writeCsv(debugOutputFileNamePrefix);
        phasingTable.writePng(debugOutputFileNamePrefix + "-RelativePhase.png",
            PhasingTable::ColoringMethod::byRelativePhase);
        phasingTable.writePng(debugOutputFileNamePrefix + "-DiscreteRelativePhase.png",
            PhasingTable::ColoringMethod::byDiscreteRelativePhase);
    }

    // Create the PhasedComponents.
    phasingTable.constructPhasedComponents(debug);


#if 1
    // Split each PhasedComponent at locations where this is necessary.
    // Check pairs of adjacent consecutive bubbles in the same phased component.
    vector< shared_ptr<PhasedComponent> > splitComponents;
    for(const auto& phasedComponentPointer: phasingTable.phasedComponents) {
        const PhasedComponent& phasedComponent = *phasedComponentPointer;
        if(phasedComponent.size() < 2) {
            break;
        }
        if(debug) {
            cout << "Checking for splitting a PhasedComponent of size " << phasedComponent.size() << endl;
        }
        vector<uint64_t> splitComponentsBegin(1, 0);
        for(uint64_t i=1; i<phasedComponent.size(); i++) {
            const auto& p0 = phasedComponent[i-1];
            const auto& p1 = phasedComponent[i];
            const uint64_t positionInBubbleChain0 = p0.first;
            const uint64_t positionInBubbleChain1 = p1.first;
            const int64_t phase0 = p0.second;
            const int64_t phase1 = p1.second;

            const Bubble& bubble0 = bubbleChain[positionInBubbleChain0];
            const Bubble& bubble1 = bubbleChain[positionInBubbleChain1];
            SHASTA_ASSERT(bubble0.isDiploid());
            SHASTA_ASSERT(bubble1.isDiploid());

            const Chain& chain00 = bubble0[0];
            const Chain& chain01 = bubble0[1];
            const Chain& chain10 = (phase0 == phase1) ? bubble1[0] : bubble1[1];
            const Chain& chain11 = (phase0 == phase1) ? bubble1[1] : bubble1[0];

            MarkerGraphEdgeId e00 = chain00.secondToLast();
            MarkerGraphEdgeId e01 = chain01.secondToLast();
            MarkerGraphEdgeId e10 = chain10.second();
            MarkerGraphEdgeId e11 = chain11.second();

            MarkerGraphEdgePairInfo info;
            SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(e00, e10, info));
            const uint64_t common0 = info.common;
            SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(e01, e11, info));
            const uint64_t common1 = info.common;

            if(debug) {
                cout << "Bubble pair: " <<
                    positionInBubbleChain0 << " " <<
                    positionInBubbleChain1 <<
                    ": side 0 " << e00 << " " << e10 << " " << common0 << " " <<
                    ", side 1 " << e01 << " " << e11 << " " << common1 << endl;
                if(common0 == 0 or common1 == 0) {
                    cout << "No common oriented reads." << endl;
                }
            }

            if(common0 == 0 or common1 == 0) {
                splitComponentsBegin.push_back(i);
            }
        }
        splitComponentsBegin.push_back(phasedComponent.size());


        // Split this phased component, if necessary.
        if(splitComponentsBegin.size() == 2) {
            // No splitting necessary.
            splitComponents.push_back(phasedComponentPointer);
            if(debug) {
                cout << "No splitting was necessary." << endl;
            }
        } else {
            // Split at the split points.
            for(uint64_t i=0; i<splitComponentsBegin.size()-1; i++) {
                const uint64_t begin = splitComponentsBegin[i];
                const uint64_t end = splitComponentsBegin[i+1];
                shared_ptr<PhasedComponent> splitComponentPointer = make_shared<PhasedComponent>();
                copy(phasedComponent.begin() + begin, phasedComponent.begin() + end,
                    back_inserter(*splitComponentPointer));
                splitComponentPointer->computePositionRange();
                splitComponents.push_back(splitComponentPointer);
                if(debug) {
                    cout << "Created a split component at " << begin << " to " << end-1 << " (inclusive)." << endl;
                }
            }
        }
    }
    phasingTable.phasedComponents.swap(splitComponents);
#endif



    // Remove PhasedComponents consisting of only one short bubble.
    {
        vector< shared_ptr<PhasedComponent> > newPhasedComponents;
        for(const auto& phasedComponent: phasingTable.phasedComponents) {
            bool keep = true;
            if(phasedComponent->size() == 1) {
                const uint64_t positionInBubbleChain = phasedComponent->front().first;
                const Bubble& bubble = bubbleChain[positionInBubbleChain];

                uint64_t averageOffset;
                uint64_t minOffset;
                uint64_t maxOffset;
                bubbleOffset(bubble, averageOffset, minOffset, maxOffset);
                if(maxOffset < longBubbleThreshold) {
                    keep = false;
                }
            }
            if(keep) {
                newPhasedComponents.push_back(phasedComponent);
            }
        }
        phasingTable.phasedComponents.swap(newPhasedComponents);
    }



    //  Use the phased components to phase the BubbleChain.
    phaseBubbleChainUsingPhasedComponents(
        debug,
        e,
        phasingTable.phasedComponents,
        longBubbleThreshold);

}



void AssemblyGraph::cleanupBubbleChainUsingPhasingTable(
    const string& debugOutputFileNamePrefix,
    edge_descriptor e,
    double phaseErrorThreshold,
    double bubbleErrorThreshold,
    uint64_t longBubbleThreshold)
{

    AssemblyGraph& cGraph = *this;
    BubbleChain& bubbleChain = cGraph[e];

    const bool debug = not debugOutputFileNamePrefix.empty();
    if(debug) {
        cout << "Before bubble clean up, bubble chain " <<
            bubbleChainStringId(e) << " has " << cGraph[e].size() << " bubbles." << endl;
    }

    // If this bubble chain has a single bubble, there is nothing to do.
    if(bubbleChain.size() == 1) {
        if(debug) {
            cout << "Skipped because it has only one bubble." << endl;
        }
        return;
    }

    // Create the phasing table for this bubble chain.
    PhasingTable phasingTable(bubbleChain, assembler.markerGraph, phaseErrorThreshold);

    if(phasingTable.empty()) {
        return;
    }
    if(phasingTable.bubbleCount() < 2) {
        return;
    }

    if(debug) {
        const uint64_t totalCount = phasingTable.entryCount();
        const uint64_t ambiguousCount = phasingTable.ambiguousEntryCount();
        const uint64_t unambiguousCount = totalCount - ambiguousCount;
        const uint64_t bubbleCount = phasingTable.bubbleCount();
        const uint64_t orientedReadCount = phasingTable.orientedReadCount();
        const double coverage = double(unambiguousCount) / double(bubbleCount);

        cout << "Phasing table summary (for bubble cleanup) " << bubbleChainStringId(e) << ":" << endl;
        cout << bubbleCount << " diploid bubbles." << endl;
        cout << orientedReadCount << " oriented reads." << endl;
        cout << unambiguousCount << " unambiguous entries." << endl;
        cout << ambiguousCount << " ambiguous entries." << endl;
        cout << "Average coverage " << std::round(coverage) << endl;
        cout << "Average number of diploid bubbles seen by each oriented read " <<
            std::round(double(unambiguousCount)/double(orientedReadCount)) << endl;
    }

    // Phasing of the phasing table.
    phasingTable.greedyPhasing();
    if(debug) {
        uint64_t consistentCount;
        uint64_t inconsistentCount;
        tie(consistentCount, inconsistentCount) = phasingTable.countConsistentEntries();

        cout << "After greedy phasing, the phasing table (for bubble cleanup) has " << consistentCount <<
            " consistent entries and " << inconsistentCount <<
            " inconsistent entries (" << consistentCount + inconsistentCount <<
            " total)." << endl;

        phasingTable.writePng(debugOutputFileNamePrefix + "-Consistency.png",
            PhasingTable::ColoringMethod::byConsistency);
        phasingTable.writeCsv(debugOutputFileNamePrefix);
        phasingTable.writePng(debugOutputFileNamePrefix + "-RelativePhase.png",
            PhasingTable::ColoringMethod::byRelativePhase);
        phasingTable.writePng(debugOutputFileNamePrefix + "-DiscreteRelativePhase.png",
            PhasingTable::ColoringMethod::byDiscreteRelativePhase);
    }


    // Use the PhasingTable to create a new BubbleChain that will replace the existing one.
    // In the new bubble chain, we remove:
    // - All diploid bubbles that have a high error rate in the PhasingTable,
    //   unless they are longer than longBubbleThreshold.
    // - All bubbles with ploidy greater than 2,
    //   unless they are longer than longBubbleThreshold.
    // Each bubble that is removed is replaced by a haploid bubble consisting
    // of only the terminal MarkerGraphEdgeIds.
    BubbleChain newBubbleChain;
    for(uint64_t positionInBubbleChain = 0; positionInBubbleChain < bubbleChain.size();
        positionInBubbleChain++) {
        const Bubble& bubble = bubbleChain[positionInBubbleChain];

        // Decide whether this Bubble will be copied verbatim to the new bubble chain.
        bool copyVerbatim = false;
        if(bubble.isHaploid()) {
            copyVerbatim = true;
            if(debug) {
                cout << "Bubble at position in bubble chain " << positionInBubbleChain <<
                    " is haploid and will be kept." << endl;
            }
        } else if(bubble.isDiploid()) {
            const double bubbleErrorRate = phasingTable.bubbleErrorRate(positionInBubbleChain);
            if(debug) {
                cout << "Bubble at phasing table index " << phasingTable.bubblesMap[positionInBubbleChain] <<
                    " position in bubble chain " << positionInBubbleChain <<
                    " has error rate " << bubbleErrorRate;
                if(bubbleErrorRate <= bubbleErrorThreshold) {
                    cout << " and will be kept." << endl;
                } else {
                    cout << " and will be removed." << endl;
                }
            }
            if(bubbleErrorRate <= bubbleErrorThreshold) {
                copyVerbatim = true;
            }
        } else {
            if(debug) {
                cout << "Bubble at position in bubble chain " << positionInBubbleChain <<
                    " has ploidy " << bubble.size() << " and will be removed." << endl;
            }
        }
        if(not copyVerbatim) {
            uint64_t averageOffset;
            uint64_t minOffset;
            uint64_t maxOffset;
            bubbleOffset(bubble, averageOffset, minOffset, maxOffset);
            copyVerbatim = maxOffset >= longBubbleThreshold;
        }

        if(copyVerbatim) {
            newBubbleChain.push_back(bubble);
            if(debug) {
                cout << "Bubble at position in bubble chain " << positionInBubbleChain <<
                    " was copied to the new bubble chain." << endl;
            }
        } else {
            // Just add a simple haploid bubble with only the source
            // and target MarkerGraphEdgeIds.
            Bubble newBubble;
            newBubble.resize(1);    // Make it haploid
            Chain& newChain = newBubble.front();    // Its only chain.
            newChain.push_back(bubble.front().front()); // Source MarkerGraphEdgeId
            newChain.push_back(bubble.front().back());  // Target MarkerGraphEdgeId
            newBubbleChain.push_back(newBubble);
            if(debug) {
                cout << "Bubble at position in bubble chain " << positionInBubbleChain <<
                    " was replaced by a simple haploid bubble in the new bubble chain: " <<
                    bubble.front().front() << " " << bubble.front().back() << endl;
            }
        }
    }
    bubbleChain = newBubbleChain;

    if(debug) {
        cout << "After bubble clean up, bubble chain " <<
            bubbleChainStringId(e) << " has " << newBubbleChain.size() <<
            " bubbles of which " <<
            newBubbleChain.diploidBubbleCount() << " diploid." << endl;
        const string csvFileName = debugOutputFileNamePrefix + "-ChainsDetails-PostBubbleCleanup.csv";
        ofstream csv(csvFileName);
        cout << "For chain details after bubble cleanup, see " << csvFileName << endl;
        writeChainDetailsCsv(csv, e, true);
    }

    // Replace the old BubbleChain with the new one, leaving the id of the edge unchanged.
    bubbleChain.compress();
    if(debug) {
        cout << "After bubble clean up and compression, bubble chain " <<
            bubbleChainStringId(e) << " has " << newBubbleChain.size() <<
            " bubbles of which " <<
            newBubbleChain.diploidBubbleCount() << " diploid." << endl;
        const string csvFileName = debugOutputFileNamePrefix +
            "-ChainsDetails-PostBubbleCleanupSAndCompress.csv";
        ofstream csv(csvFileName);
        cout << "For chain details after bubble cleanup and compress, see " << csvFileName << endl;
        writeChainDetailsCsv(csv, e, true);
    }
}



// Compute the tangle matrix between two incoming chains
// and two outgoing chains, taking into account up to
// n MarkergraphEdgeIds for each Chain.
void AssemblyGraph::computeTangleMatrix(
    const array<const Chain*, 2> inChains,
    const array<const Chain*, 2> outChains,
    uint64_t n,
    TangleMatrix& tangleMatrix) const
{
    // Gather the OrientedReadIds near the end of the inChains.
    array<vector<OrientedReadId>, 2> allOrientedReadIdsIn;
    for(uint64_t i=0; i<2; i++) {
        gatherOrientedReadIdsAtEnd(*inChains[i], n, allOrientedReadIdsIn[i]);

    }

    // Gather the OrientedReadIds near the beginning of the outChains.
    array<vector<OrientedReadId>, 2> allOrientedReadIdsOut;
    for(uint64_t i=0; i<2; i++) {
        gatherOrientedReadIdsAtBeginning(*outChains[i], n, allOrientedReadIdsOut[i]);
    }

    // Discard OrientedReadIds that appear in both inChains.
    array<vector<OrientedReadId>, 2> orientedReadIdsIn;
    for(uint64_t i=0; i<2; i++) {
        std::set_difference(
            allOrientedReadIdsIn[i]  .begin(), allOrientedReadIdsIn[i]  .end(),
            allOrientedReadIdsIn[1-i].begin(), allOrientedReadIdsIn[1-i].end(),
            back_inserter(orientedReadIdsIn[i]));
    }

    // Discard OrientedReadIds that appear in both outChains.
    array<vector<OrientedReadId>, 2> orientedReadIdsOut;
    for(uint64_t i=0; i<2; i++) {
        std::set_difference(
            allOrientedReadIdsOut[i]  .begin(), allOrientedReadIdsOut[i]  .end(),
            allOrientedReadIdsOut[1-i].begin(), allOrientedReadIdsOut[1-i].end(),
            back_inserter(orientedReadIdsOut[i]));
    }

    // Now we can compute the tangle matrix.
    vector<OrientedReadId> commonOrientedReads;
    for(uint64_t i0=0; i0<2; i0++) {
        for(uint64_t i1=0; i1<2; i1++) {
            commonOrientedReads.clear();
            set_intersection(
                orientedReadIdsIn[i0] .begin(), orientedReadIdsIn[i0] .end(),
                orientedReadIdsOut[i1].begin(), orientedReadIdsOut[i1].end(),
                back_inserter(commonOrientedReads));
            tangleMatrix[i0][i1] = commonOrientedReads.size();
        }
    }
}



// Gather OrientedReadIds from up to n MarkergraphEdgeIds
// near the end of a chain.
void AssemblyGraph::gatherOrientedReadIdsAtEnd(
    const Chain& chain,
    uint64_t n,
    vector<OrientedReadId>& orientedReadIds) const
{

    const uint64_t last = chain.size() - 2;                     // Exclude last MarkergraphEdgeId.
    const uint64_t first = (last > (n-1)) ? last + 1 - n : 0;   // Use up to n.

    SHASTA_ASSERT(first < chain.size());
    SHASTA_ASSERT(last < chain.size());

    orientedReadIds.clear();
    for(uint64_t i=first; i<=last; i++) {
        const MarkerGraphEdgeId markerGraphEdgeId = chain[i];
        const auto& markerIntervals =
            assembler.markerGraph.edgeMarkerIntervals[markerGraphEdgeId];
        for(const MarkerInterval& markerInterval: markerIntervals) {
            orientedReadIds.push_back(markerInterval.orientedReadId);
        }
    }
    deduplicate(orientedReadIds);
}



// Gather OrientedReadIds from up to n MarkergraphEdgeIds
// near the beginning of a chain.
void AssemblyGraph::gatherOrientedReadIdsAtBeginning(
    const Chain& chain,
    uint64_t n,
    vector<OrientedReadId>& orientedReadIds) const
{

    const uint64_t first = 1;   // / Exclude first MarkergraphEdgeId.
    const uint64_t last = (chain.size() > (n+1)) ? n : chain.size() - 1;

    SHASTA_ASSERT(first < chain.size());
    SHASTA_ASSERT(last < chain.size());

    orientedReadIds.clear();
    for(uint64_t i=first; i<=last; i++) {
        const MarkerGraphEdgeId markerGraphEdgeId = chain[i];
        const auto& markerIntervals =
            assembler.markerGraph.edgeMarkerIntervals[markerGraphEdgeId];
        for(const MarkerInterval& markerInterval: markerIntervals) {
            orientedReadIds.push_back(markerInterval.orientedReadId);
        }
    }
    deduplicate(orientedReadIds);
}



// To phase the PhasingGraph, we create an optimal spanning tree
// using edges in order of decreasing "significance".
void AssemblyGraph::PhasingGraph::phase(bool debug)
{
    PhasingGraph& phasingGraph = *this;

    // Gather edges by maxDiscordant and minConcordant.
    // edgeTable[maxDiscordant][minConcordant] contains the
    // edges with those values of maxDiscordant and minConcordant.
    // This allows the code later ot process edges in order
    // of increasing maxDiscordant and decreasing minConcordant.
    vector< vector< vector<edge_descriptor> > > edgeTable;
    BGL_FORALL_EDGES(pe, phasingGraph, PhasingGraph) {
        const PhasingGraphEdge& edge = phasingGraph[pe];
        const uint64_t maxDiscordant = edge.maxDiscordant;
        const uint64_t minConcordant = edge.minConcordant;
        if(edgeTable.size() <= maxDiscordant) {
            edgeTable.resize(maxDiscordant + 1);
        }
        vector< vector<edge_descriptor> >& v = edgeTable[maxDiscordant];
        if(v.size() <= minConcordant) {
            v.resize(minConcordant + 1);
        }
        v[minConcordant].push_back(pe);
    }

    // Map vertices to integers.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(pv, phasingGraph, PhasingGraph) {
        vertexIndexMap.insert({pv, vertexIndex++});
    }
    const uint64_t vertexCount = vertexIndexMap.size();



    // Compute optimal spanning tree and connected components.
    vector<uint64_t> rank(vertexCount);
    vector<uint64_t> parent(vertexCount);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t i=0; i<vertexCount; i++) {
        disjointSets.make_set(i);
    }
    uint64_t spanningTreeEdgeCount = 0;
    for(uint64_t maxDiscordant=0; maxDiscordant<edgeTable.size(); maxDiscordant++) {
        const vector< vector<edge_descriptor> >& v = edgeTable[maxDiscordant];
        for(int64_t minConcordant=v.size()-1; minConcordant>=0; minConcordant--) {
            const vector<edge_descriptor>& vv = v[minConcordant];
            if(false) {
                cout << "Processing " << vv.size() << " phasing graph edges with maxDiscordant=" <<
                    maxDiscordant << ", minConcordant=" << minConcordant << endl;
            }
            for(const edge_descriptor e: vv) {
                PhasingGraphEdge& edge = phasingGraph[e];
                const vertex_descriptor pv0 = source(e, phasingGraph);
                const vertex_descriptor pv1 = target(e, phasingGraph);
                const uint64_t vertexIndex0 = vertexIndexMap[pv0];
                const uint64_t vertexIndex1 = vertexIndexMap[pv1];
                const uint64_t componentId0 = disjointSets.find_set(vertexIndex0);
                const uint64_t componentId1 = disjointSets.find_set(vertexIndex1);
                if(componentId0 != componentId1) {
                    disjointSets.union_set(vertexIndex0, vertexIndex1);
                    edge.isSpanningTreeEdge = true;
                    ++spanningTreeEdgeCount;
                }
            }
            if(false) {
                cout << "Found " << spanningTreeEdgeCount << " spanning tree edges so far." << endl;
            }
        }
    }

    // Gather the vertices in each connected component.
    vector< vector<vertex_descriptor> > components(vertexCount);
    BGL_FORALL_VERTICES(pv, phasingGraph, PhasingGraph) {
        const uint64_t componentId = disjointSets.find_set(vertexIndexMap[pv]);
        components[componentId].push_back(pv);
    }

    // Write a histogram of component sizes.
    if(debug) {
        vector<uint64_t> histogram;
        for(const vector<vertex_descriptor>& component: components) {
            const uint64_t componentSize = component.size();
            if(histogram.size() <= componentSize) {
                histogram.resize(componentSize + 1, 0);
            }
            ++histogram[componentSize];
        }

        cout << "Histogram of component sizes:" << endl;
        cout << "Size,Frequency,Vertices" << endl;
        for(uint64_t componentSize=1; componentSize<histogram.size(); componentSize++) {
            const uint64_t frequency = histogram[componentSize];
            if(frequency) {
                cout << componentSize << "," << frequency << ","  << componentSize*frequency << endl;
            }
        }
    }

    // Gather the non-trivial component and sort them by decreasing size.
    vector< pair<uint64_t, uint64_t> > componentTable; // (componentId, componentSize)
    for(uint64_t componentId=0; componentId<vertexCount; componentId++) {
        const vector<vertex_descriptor>& component = components[componentId];
        if(component.size() > 1) {
            componentTable.push_back({componentId, component.size()});
        }
    }
    sort(componentTable.begin(), componentTable.end(),
        OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());



    // Process the non-trivial components in order of decreasing size.
    phasedComponents.clear();
    for(const  pair<uint64_t, uint64_t>& p: componentTable) {
        const uint64_t componentId = p.first;
        const vector<vertex_descriptor>& component = components[componentId];
        SHASTA_ASSERT(component.size() == p.second);
        if(debug) {
            cout << "Processing a phasing component with " << component.size() <<
                " vertices." << endl;
        }

        // Use a BFS on the spanning tree to phase the vertices in this component.
        // Use the spanning tree to phase vertices in the largest component.
        // It does not matter which vertex we start from.
        const vertex_descriptor vFirst = component.front();
        phasingGraph[vFirst].phase = +1;
        std::queue<vertex_descriptor> q;
        q.push(vFirst);
        while(not q.empty()) {
            const vertex_descriptor v0 = q.front();
            q.pop();
            BGL_FORALL_OUTEDGES(v0, e, phasingGraph, PhasingGraph) {
                PhasingGraphEdge& edge = phasingGraph[e];
                if(not edge.isSpanningTreeEdge) {
                    continue;
                }
                const PhasingGraphVertex& vertex0 = phasingGraph[v0];
                const vertex_descriptor v1 = target(e, phasingGraph);
                PhasingGraphVertex& vertex1 = phasingGraph[v1];
                if(vertex1.phase == 0) {
                    vertex1.phase = vertex0.phase;
                    if(edge.phase == -1) {
                        vertex1.phase = - vertex1.phase;
                    }
                    q.push(v1);
                }
            }
        }

        // Count inconsistent edges in this component.
        if(debug) {
            uint64_t inconsistentCount = 0;
            uint64_t totalCount = 0;
            for(const vertex_descriptor v: component) {
                BGL_FORALL_OUTEDGES(v, e, phasingGraph, PhasingGraph) {
                    totalCount++;
                    if(not isConsistent(e)) {
                        ++inconsistentCount;
                    }
                }
            }
            // This counts edges twice.
            inconsistentCount /= 2;
            totalCount /= 2;
            cout << inconsistentCount << " inconsistent edges in this component out of " <<
                totalCount << " total." << endl;
        }


        // Create the PhasedComponent corresponding to this component.
        // Don't include any vertices that overlap previous PhasedComponent.
        shared_ptr<PhasedComponent> phasedComponentPointer = make_shared<PhasedComponent>();
        PhasedComponent& phasedComponent = *phasedComponentPointer;
        for(const vertex_descriptor pv: component) {
            const PhasingGraphVertex& vertex = phasingGraph[pv];
            const uint64_t positionInBubbleChain = vertex.positionInBubbleChain;
            bool overlapsPrevious = false;
            for(const auto& phasedComponent: phasedComponents) {
                if(
                    positionInBubbleChain >= phasedComponent->minPositionInBubbleChain and
                    positionInBubbleChain <= phasedComponent->maxPositionInBubbleChain) {
                    overlapsPrevious = true;
                    break;
                }
            }
            if(not overlapsPrevious) {
                phasedComponent.push_back({vertex.positionInBubbleChain, vertex.phase});
            }
        }
        if(phasedComponent.size() < 2) {
            if(debug) {
                cout << "This component will be discarded due to overlap with previous components." << endl;
            }
            continue;
        }
        phasedComponent.sort();
        if(debug) {
            cout << "Phasing range for this component " << phasedComponent.minPositionInBubbleChain <<
                " " << phasedComponent.maxPositionInBubbleChain << endl;
        }
        phasedComponents.push_back(phasedComponentPointer);
    }

    // Sort the phased components in order of increasing position.
    class SortHelper {
    public:
        bool operator()(
            const shared_ptr<PhasedComponent>& p0,
            const shared_ptr<PhasedComponent>& p1
            ) const
        {
            return p0->minPositionInBubbleChain < p1->minPositionInBubbleChain;
        }
    };
    sort(phasedComponents.begin(), phasedComponents.end(), SortHelper());

    if(debug) {
        cout << "Kept " << phasedComponents.size() << " phased components:" << endl;
        for(const auto& phasedComponent: phasedComponents) {
            cout  << phasedComponent->size() << " diploid bubbles at positions " <<
                phasedComponent->minPositionInBubbleChain << "..." <<
                phasedComponent->maxPositionInBubbleChain << " in bubble chain." << endl;

        }
        phasingGraph.writeGraphviz("PhasingGraph.dot");
    }
}



// Sort edges in order of decreasing significance:
// - If using the Bayesian model, logP.
// - Otherwise, minConcordant/maxDiscordant.
void AssemblyGraph::PhasingGraph::sortEdges(
    vector<edge_descriptor>& sortedEdges,
    bool useBayesianModel) const
{
    const PhasingGraph& phasingGraph = *this;

    if(useBayesianModel) {

        // Gather edges and their logP.
        vector< pair<edge_descriptor, double> > edgeTable;
        BGL_FORALL_EDGES(pe, phasingGraph, PhasingGraph) {
            const PhasingGraphEdge& edge = phasingGraph[pe];
            edgeTable.push_back({pe, edge.logP()});
        }

        // Sort by decreasing logP.
        sort(edgeTable.begin(), edgeTable.end(),
            OrderPairsBySecondOnlyGreater<edge_descriptor, double>());
        sortedEdges.clear();
        for(const auto& p: edgeTable) {
            sortedEdges.push_back(p.first);
        }

    } else {

        // Gather edges by maxDiscordant and minConcordant.
        // edgeTable[maxDiscordant][minConcordant] contains the
        // edges with those values of maxDiscordant and minConcordant.
        vector< vector< vector<edge_descriptor> > > edgeTable;
        BGL_FORALL_EDGES(pe, phasingGraph, PhasingGraph) {
            const PhasingGraphEdge& edge = phasingGraph[pe];
            const uint64_t maxDiscordant = edge.maxDiscordant;
            const uint64_t minConcordant = edge.minConcordant;
            if(edgeTable.size() <= maxDiscordant) {
                edgeTable.resize(maxDiscordant + 1);
            }
            vector< vector<edge_descriptor> >& v = edgeTable[maxDiscordant];
            if(v.size() <= minConcordant) {
                v.resize(minConcordant + 1);
            }
            v[minConcordant].push_back(pe);
        }

        // The sorted edges are in order of increasing maxDiscordant
        // and decreasing minConcordant.
        sortedEdges.clear();
        for(uint64_t maxDiscordant=0; maxDiscordant<edgeTable.size(); maxDiscordant++) {
            const vector< vector<edge_descriptor> >& v = edgeTable[maxDiscordant];
            for(int64_t minConcordant=v.size()-1; minConcordant>=0; minConcordant--) {
                const vector<edge_descriptor>& vv = v[minConcordant];
                for(const edge_descriptor e: vv) {
                    sortedEdges.push_back(e);
                }
            }
        }

    }
}



// To phase the PhasingGraph, we create an optimal spanning tree
// using edges in order of decreasing "significance".
// We do this iteratively. At each iteration we process the largest
// connected component of the surviving PhasingGraph.
void AssemblyGraph::PhasingGraph::phase1(bool debug, bool useBayesianModel)
{
    PhasingGraph& phasingGraph = *this;
    phasedComponents.clear();

    if(debug) {
        cout << "Beginning phasing for a PhasingGraph with " << num_vertices(phasingGraph) <<
            " vertices." << endl;
    }

    // Main iteration loop.
    while(true) {

        // Clear the isSpanningTreeEdge flag of all edges.
        BGL_FORALL_EDGES(pe, phasingGraph, PhasingGraph) {
            phasingGraph[pe].isSpanningTreeEdge = false;
        }

        // Sort edges in order of decreasing significance:
        // - If using the Bayesian model, logP.
        // - Otherwise, minConcordant/maxDiscordant.
        vector<edge_descriptor> sortedEdges;
        sortEdges(sortedEdges, useBayesianModel);

        // Map vertices to integers.
        // This is needed for the computation of the spanning tree and
        // connected components.
        std::map<vertex_descriptor, uint64_t> vertexIndexMap;
        uint64_t vertexIndex = 0;
        BGL_FORALL_VERTICES(pv, phasingGraph, PhasingGraph) {
            vertexIndexMap.insert({pv, vertexIndex++});
        }
        const uint64_t vertexCount = vertexIndexMap.size();

        if(debug) {
            cout << "Beginning a new phasing iteration. The phasing graph has " <<
                vertexCount << " vertices left." << endl;
        }



        // Compute optimal spanning tree and connected components.
        vector<uint64_t> rank(vertexCount);
        vector<uint64_t> parent(vertexCount);
        boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
        for(uint64_t i=0; i<vertexCount; i++) {
            disjointSets.make_set(i);
        }
        uint64_t spanningTreeEdgeCount = 0;

        for(const edge_descriptor e: sortedEdges) {
            PhasingGraphEdge& edge = phasingGraph[e];
            const vertex_descriptor pv0 = source(e, phasingGraph);
            const vertex_descriptor pv1 = target(e, phasingGraph);
            const uint64_t vertexIndex0 = vertexIndexMap[pv0];
            const uint64_t vertexIndex1 = vertexIndexMap[pv1];
            const uint64_t componentId0 = disjointSets.find_set(vertexIndex0);
            const uint64_t componentId1 = disjointSets.find_set(vertexIndex1);
            if(componentId0 != componentId1) {
                disjointSets.union_set(vertexIndex0, vertexIndex1);
                edge.isSpanningTreeEdge = true;
                ++spanningTreeEdgeCount;
            }
        }

        // Gather the vertices in each connected component.
        vector< vector<vertex_descriptor> > components(vertexCount);
        BGL_FORALL_VERTICES(pv, phasingGraph, PhasingGraph) {
            const uint64_t componentId = disjointSets.find_set(vertexIndexMap[pv]);
            components[componentId].push_back(pv);
        }

        // Find the largest connected component.
        uint64_t largestComponentId = invalid<uint64_t>;
        uint64_t largestComponentSize = 0;
        for(uint64_t componentId=0; componentId<vertexCount; componentId++) {
            const uint64_t componentSize = components[componentId].size();
            if(componentSize > largestComponentSize) {
                largestComponentSize = componentSize;
                largestComponentId = componentId;
            }
        }

        // If the largest component has less than two vertices, we are done.
        if(largestComponentSize < 2) {
            if(debug) {
                cout << "Phasing terminates because  only trivial connected components were found." << endl;
            }
            break;
        }

        // Access the largest connected component, which we will be working on
        // for the rest of this iteration.
        const vector<vertex_descriptor>& component = components[largestComponentId];
        SHASTA_ASSERT(component.size() == largestComponentSize);
        if(debug) {
            cout << "The largest component of the current PhasingGraph has " <<
                largestComponentSize << " vertices." << endl;
        }

        // Use a BFS on the spanning tree to phase the vertices in this component.
        // It does not matter which vertex we start from.
        const vertex_descriptor vFirst = component.front();
        phasingGraph[vFirst].phase = +1;
        std::queue<vertex_descriptor> q;
        q.push(vFirst);
        while(not q.empty()) {
            const vertex_descriptor v0 = q.front();
            q.pop();
            BGL_FORALL_OUTEDGES(v0, e, phasingGraph, PhasingGraph) {
                PhasingGraphEdge& edge = phasingGraph[e];
                if(not edge.isSpanningTreeEdge) {
                    continue;
                }
                const PhasingGraphVertex& vertex0 = phasingGraph[v0];
                const vertex_descriptor v1 = target(e, phasingGraph);
                PhasingGraphVertex& vertex1 = phasingGraph[v1];
                if(vertex1.phase == 0) {
                    vertex1.phase = vertex0.phase;
                    if(edge.phase == -1) {
                        vertex1.phase = - vertex1.phase;
                    }
                    q.push(v1);
                }
            }
        }

        // Count inconsistent edges in this component.
        if(debug) {
            uint64_t inconsistentCount = 0;
            uint64_t totalCount = 0;
            for(const vertex_descriptor v: component) {
                BGL_FORALL_OUTEDGES(v, e, phasingGraph, PhasingGraph) {
                    totalCount++;
                    if(not isConsistent(e)) {
                        ++inconsistentCount;
                    }
                }
            }
            // This counts edges twice.
            inconsistentCount /= 2;
            totalCount /= 2;
            cout << inconsistentCount << " inconsistent edges in this component out of " <<
                totalCount << " total." << endl;
        }

        // All vertices in this component have been phased.
        // However, when creating the PhasedComponent, we have to make sure that adjacent
        // phased vertices have common reads.
        // To guarantee this, we find a longest path in this component, in order of increasing
        // positionInBubbleChain. Only vertices in this longest path are then included in the
        // PhasedComponent.

        // To find this longest path, we use an algorithm similar to the one in longestPath.cpp,
        // using the topological ordering induced by positionInBubbleChain.

        // Table of the vertices in order of increasing positionInBubbleChain.
        vector< pair<vertex_descriptor, uint64_t> > vertexTable;
        for(const vertex_descriptor v: component) {
            vertexTable.push_back({v, phasingGraph[v].positionInBubbleChain});
        }
        sort(vertexTable.begin(), vertexTable.end(), OrderPairsBySecondOnly<vertex_descriptor, uint64_t>());

        // The length of the longest path ending at each vertex.
        std::map<vertex_descriptor, uint64_t> lengthMap;
        for(const vertex_descriptor v: component) {
            lengthMap.insert(make_pair(v, 0));
        }

        // Process the vertices in order of increasing positionInBubbleChain.
        for(const auto& p: vertexTable) {
            const vertex_descriptor v0 = p.first;
            const uint64_t positionInBubbleChain0 = phasingGraph[v0].positionInBubbleChain;

            uint64_t maximumLength = 0;
            BGL_FORALL_OUTEDGES_T(v0, e, phasingGraph, PhasingGraph) {
                const vertex_descriptor v1 = target(e, phasingGraph);
                const uint64_t positionInBubbleChain1 = phasingGraph[v1].positionInBubbleChain;

                if(positionInBubbleChain1 < positionInBubbleChain0) {
                    maximumLength = max(maximumLength, lengthMap[v1]);
                }
            }
            lengthMap[v0] = maximumLength + 1;
        }

        // Find the vertex with the longest length.
        // This will be the end of the longest path.
        vertex_descriptor v = PhasingGraph::null_vertex();
        uint64_t maximumLength = 0;
        for(const auto& p: lengthMap) {
            if(p.second > maximumLength) {
                v = p.first;
                maximumLength = p.second;
            }
        }

        // Constuct the path, moving backward from here.
        vector<vertex_descriptor> longestPath;
        longestPath.push_back(v);
        while(true) {
            vertex_descriptor vPrevious = PhasingGraph::null_vertex();
            uint64_t maximumLength = 0;
            BGL_FORALL_OUTEDGES(v, e, phasingGraph, PhasingGraph) {
                const vertex_descriptor v0 = target(e, phasingGraph);
                if(phasingGraph[v0].positionInBubbleChain < phasingGraph[v].positionInBubbleChain) {
                    const uint64_t length = lengthMap[v0];
                    if(length > maximumLength) {
                        vPrevious = v0;
                        maximumLength = length;
                    }
                }
            }
            if(vPrevious == PhasingGraph::null_vertex()) {
                break;
            }
            v = vPrevious;
            longestPath.push_back(v);

        }
        std::reverse(longestPath.begin(), longestPath.end());

        if(debug) {
            cout << "The longest path contains " << longestPath.size() << " vertices." << endl;
        }



        // If the longest path is non-trivial, use it to create a new PhasedComponent.
        if(longestPath.size() > 1) {
            if(debug) {
                cout << "Creating a new PhasedComponent." << endl;
            }
            shared_ptr<PhasedComponent> phasedComponentPointer = make_shared<PhasedComponent>();
            phasedComponents.push_back(phasedComponentPointer);
            PhasedComponent& phasedComponent = *phasedComponentPointer;

            for(const vertex_descriptor v: longestPath) {
                const PhasingGraphVertex& vertex = phasingGraph[v];
                phasedComponent.push_back({vertex.positionInBubbleChain, vertex.phase});
            }
            phasedComponent.minPositionInBubbleChain = phasingGraph[longestPath.front()].positionInBubbleChain;
            phasedComponent.maxPositionInBubbleChain = phasingGraph[longestPath.back()].positionInBubbleChain;
            if(debug) {
                cout << "Phasing range for this component " << phasedComponent.minPositionInBubbleChain <<
                    " " << phasedComponent.maxPositionInBubbleChain << endl;
            }

            // Now remove from the PhasingGraph all vertices of this component
            // plus any vertices with a positionInBubbleChain
            // that overlaps this phased component.
            vector<vertex_descriptor> verticesToBeRemoved = component;
            BGL_FORALL_VERTICES(v, phasingGraph, PhasingGraph) {
                const uint64_t positionInBubbleChain = phasingGraph[v].positionInBubbleChain;
                if( positionInBubbleChain >= phasedComponent.minPositionInBubbleChain and
                    positionInBubbleChain <= phasedComponent.maxPositionInBubbleChain) {
                    verticesToBeRemoved.push_back(v);
                }
            }
            deduplicate(verticesToBeRemoved);
            for(const vertex_descriptor v: verticesToBeRemoved) {
                clear_vertex(v, phasingGraph);
                remove_vertex(v, phasingGraph);
            }
        } else {

            // Now remove from the PhasingGraph all vertices of this component.
            for(const vertex_descriptor v: component) {
                clear_vertex(v, phasingGraph);
                remove_vertex(v, phasingGraph);
            }
        }
    }



    // Sort the phased components in order of increasing position.
    class SortHelper {
    public:
        bool operator()(
            const shared_ptr<PhasedComponent>& p0,
            const shared_ptr<PhasedComponent>& p1
            ) const
        {
            return p0->minPositionInBubbleChain < p1->minPositionInBubbleChain;
        }
    };
    sort(phasedComponents.begin(), phasedComponents.end(), SortHelper());

    if(debug) {
        cout << phasedComponents.size() << " phased components:" << endl;
        for(const auto& phasedComponent: phasedComponents) {
            cout  << phasedComponent->size() << " diploid bubbles at positions " <<
                phasedComponent->minPositionInBubbleChain << "..." <<
                phasedComponent->maxPositionInBubbleChain << " in bubble chain." << endl;

        }
        // phasingGraph.writeGraphviz("PhasingGraph.dot");
    }
}



bool AssemblyGraph::PhasingGraph::isConsistent(edge_descriptor e) const
{
    const PhasingGraph& phasingGraph = *this;
    const vertex_descriptor v0 = source(e, phasingGraph);
    const vertex_descriptor v1 = target(e, phasingGraph);
    const int64_t phase0 = phasingGraph[v0].phase;
    const int64_t phase1 = phasingGraph[v1].phase;
    const int64_t phase = phasingGraph[e].phase;

    SHASTA_ASSERT(phase0==+1 or phase0==-1);
    SHASTA_ASSERT(phase1==+1 or phase1==-1);
    SHASTA_ASSERT(phase==+1 or phase==-1);

    if(phase == +1) {
        return phase0 == phase1;
    } else {
        return phase0 != phase1;
    }
}



void AssemblyGraph::PhasingGraph::writeGraphviz(const string& fileName) const
{
    const PhasingGraph& phasingGraph = *this;

    ofstream dot(fileName);
    dot << "graph PhasingGraph {\n";

    BGL_FORALL_EDGES(e, phasingGraph, PhasingGraph) {
        const vertex_descriptor v0 = source(e, phasingGraph);
        const vertex_descriptor v1 = target(e, phasingGraph);
        dot <<
            phasingGraph[v0].positionInBubbleChain << "--" <<
            phasingGraph[v1].positionInBubbleChain;
        if(phasingGraph[e].isSpanningTreeEdge) {
            dot << " [color=green]";
        } else  if(not isConsistent(e)) {
            dot << " [color=red]";
        }
        dot << ";\n";
    }

    dot << "}\n";
}



void AssemblyGraph::TangleMatrix::analyze(
    uint64_t lowThreshold,
    uint64_t highThreshold,
    int64_t& phase,
    uint64_t& minConcordant,
    uint64_t& maxDiscordant,
    uint64_t& total,
    double epsilon,
    double& logPin, // log[P(in-phase)/P(random)] in decibels
    double& logPout // log[P(out-of-phase)/P(random)] in decibels
    ) const
{
    const TangleMatrix& m = *this;

    // Classify matrix elements:
    // 0 = low (<=lowThreshold)
    // 1 = ambiguous (>lowThreshold, <highThreshold)
    // 2 = high (>=highThreshold)
    array< array<uint64_t, 2>, 2> c;
    total = 0;
    for(uint64_t i=0; i<2; i++) {
        for(uint64_t j=0; j<2; j++) {
            const uint64_t matrixElement = m[i][j];
            total += matrixElement;
            uint64_t& classification = c[i][j];
            if(matrixElement <= lowThreshold) {
                classification = 0;
            } else if(matrixElement >= highThreshold) {
                classification = 2;
            } else {
                classification = 1;
            }
        }
    }

    // Check if this tangle matrix is unambiguously in phase.
    if(c[0][0]==2 and c[1][1]==2 and c[0][1]==0 and c[1][0]==0) {
        phase = +1;
        minConcordant = min(m[0][0], m[1][1]);
        maxDiscordant = max(m[0][1], m[1][0]);
    }

    // Check if this tangle matrix is unambiguously out of phase.
    else if(c[0][1]==2 and c[1][0]==2 and c[0][0]==0 and c[1][1]==0) {
        phase = -1;
        minConcordant = min(m[0][1], m[1][0]);
        maxDiscordant = max(m[0][0], m[0][0]);
    }

    // Otherwise, it is ambiguous.
    else {
        phase = 0;
        minConcordant = 0;
        maxDiscordant = 0;
    }

    tie(logPin, logPout) = diploidBayesianPhase(m, epsilon);
}



// Collapse consecutive haploid bubbles of a BubbleChain.
void BubbleChain::compress()
{
    BubbleChain& bubbleChain = *this;
    BubbleChain newBubbleChain;

    // Find sets of consecutive haploid bubbles.
    for(uint64_t i=0; i<size(); i++) {
        const Bubble& bubble = bubbleChain[i];

        if(bubble.isHaploid()) {

            // This bubble is haploid.
            // If the last bubble of the new bubble is haploid, append it to that.
            // Otherwise apppend it to the last bubble.
            if(not newBubbleChain.empty() and newBubbleChain.back().isHaploid()) {
                const Chain& chain = bubble.front();
                Chain& newChain = newBubbleChain.back().front();
                copy(chain.begin()+1, chain.end(), back_inserter(newChain));
            } else {
                newBubbleChain.push_back(bubble);
            }
        } else {

            // This bubble is not haploid. Just append it to the last bubble.
            newBubbleChain.push_back(bubble);
        }

    }

    // Replace it with the new one.
    bubbleChain = newBubbleChain;
}



void AssemblyGraph::assembleChain(
    Chain& chain,
    uint64_t chainTerminalCommonThreshold)
{
    chain.stepSequences.resize(chain.size() - 1);

    // Do all the assembly steps.
    for(uint64_t positionInChain=0; positionInChain<chain.size()-1; positionInChain++) {
        runAssemblyStep(chain, positionInChain, chainTerminalCommonThreshold);
    }

    combineStepSequences(chain);
    chain.wasAssembled = true;
}



// Multithreaded version of sequence assembly.
// This only assembles the chains that have the shouldBeAssembled flag set.
void AssemblyGraph::assembleChainsMultithreaded(
    uint64_t chainTerminalCommonThreshold,
    uint64_t threadCount)
{
    AssemblyGraph& assemblyGraph = *this;

    // Store the argument so the threads can see it.
    assembleChainsMultithreadedData.chainTerminalCommonThreshold = chainTerminalCommonThreshold;

    // Gather AssemblySteps for all the Chains.
    auto& assemblySteps = assembleChainsMultithreadedData.assemblySteps;
    assemblySteps.clear();

    // Loop over BubbleChains.
    AssemblyStep assemblyStep;
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        assemblyStep.e = e;
        BubbleChain& bubbleChain = assemblyGraph[e];

        // Loop over Bubbles in this BubbleChain.
        for(uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size(); positionInBubbleChain++) {
            assemblyStep.positionInBubbleChain = positionInBubbleChain;
            Bubble& bubble = bubbleChain[positionInBubbleChain];

            // Loop over Chains in this Bubble.
            for(uint64_t indexInBubble=0; indexInBubble<bubble.size(); indexInBubble++) {
                assemblyStep.indexInBubble = indexInBubble;
                Chain& chain = bubble[indexInBubble];
                SHASTA_ASSERT(chain.size() >= 2);

                // If this Chain is not marked to be assembled, skip it.
                if(not chain.shouldBeAssembled) {
                    continue;
                }

                // Prepare the vectors where the threads will store
                // the internal sequence assembled for each AssemblyStep.
                // Each of these vectors will be modified by only one thread.
                chain.stepSequences.resize(chain.size() - 1);

                // Loop over pairs of consecutive vertices in this Chain.
                for(uint64_t positionInChain=0; positionInChain<chain.size()-1; positionInChain++) {
                    assemblyStep.positionInChain = positionInChain;

                    // Compute the offset.
                    const MarkerGraphEdgeId edgeIdA = chain[positionInChain];
                    const MarkerGraphEdgeId edgeIdB = chain[positionInChain + 1];
                    MarkerGraphEdgePairInfo info;
                    SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(
                        edgeIdA, edgeIdB, info));
                    assemblyStep.offsetInBases = info.offsetInBases;

                    // Store this assembly step.
                    assemblySteps.push_back(assemblyStep);
                }
            }
        }
    }

    // For better load balancing, sort them by decreasing offset.
    sort(assemblySteps.begin(), assemblySteps.end());



    // Assemble the steps in parallel.
    setupLoadBalancing(assemblySteps.size(),  1);
    performanceLog << timestamp << "Sequence assembly begins." << endl;
    runThreads(&AssemblyGraph::assembleChainsMultithreadedTheadFunction, threadCount);
    performanceLog << timestamp << "Sequence assembly ends." << endl;



    // Now that all the AssemblySteps have been computed, the stepSequences
    // of each Chain have been filled in.
    // Combine those with the marker graph edge sequences to obtain the
    // complete sequence of each chain.
    // This can be parallelized.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        assemblyStep.e = e;
        BubbleChain& bubbleChain = assemblyGraph[e];

        // Loop over Bubbles in this BubbleChain.
        for(uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size(); positionInBubbleChain++) {
            assemblyStep.positionInBubbleChain = positionInBubbleChain;
            Bubble& bubble = bubbleChain[positionInBubbleChain];

            // Loop over Chains in this Bubble.
            for(uint64_t indexInBubble=0; indexInBubble<bubble.size(); indexInBubble++) {
                assemblyStep.indexInBubble = indexInBubble;
                Chain& chain = bubble[indexInBubble];
                if(chain.shouldBeAssembled) {
                    combineStepSequences(chain);
                    chain.wasAssembled = true;
                }
            }
        }
    }
}



// This sets the shouldBeAssembled flag for all chains, then
// calls assembleChainsMultithreaded.
void AssemblyGraph::assembleAllChainsMultithreaded(
    uint64_t chainTerminalCommonThreshold,
    uint64_t threadCount)
{
    AssemblyGraph& assemblyGraph = *this;

    // Loop over all bubble chains.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        BubbleChain& bubbleChain = assemblyGraph[e];

        // Loop over Bubbles in this BubbleChain.
        for(uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size(); positionInBubbleChain++) {
            Bubble& bubble = bubbleChain[positionInBubbleChain];

            // Loop over Chains in this Bubble.
            for(uint64_t indexInBubble=0; indexInBubble<bubble.size(); indexInBubble++) {
                Chain& chain = bubble[indexInBubble];
                chain.shouldBeAssembled = true;
            }
        }
    }

    assembleChainsMultithreaded(chainTerminalCommonThreshold, threadCount);
}



// This clears the shouldBeAssembled flag from all Chains.
void AssemblyGraph::clearAllShouldBeAssembledFlags()
{
    AssemblyGraph& assemblyGraph = *this;

    // Loop over all bubble chains.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        BubbleChain& bubbleChain = assemblyGraph[e];

        // Loop over Bubbles in this BubbleChain.
        for(uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size(); positionInBubbleChain++) {
            Bubble& bubble = bubbleChain[positionInBubbleChain];

            // Loop over Chains in this Bubble.
            for(uint64_t indexInBubble=0; indexInBubble<bubble.size(); indexInBubble++) {
                Chain& chain = bubble[indexInBubble];
                chain.shouldBeAssembled = false;
            }
        }
    }

}



// Combine stepSequences of a Chain with the marker graph edge sequences to obtain the
// complete sequence of the chain.
void AssemblyGraph::combineStepSequences(Chain& chain)
{
    chain.sequence.clear();
    for(uint64_t positionInChain=0; /* Check later */ ; positionInChain++) {

        // Add the sequence for the marker graph primary edge.
        const MarkerGraphEdgeId edgeId = chain[positionInChain];
        const auto edgeSequence = assembler.markerGraph.edgeSequence[edgeId];
        copy(edgeSequence.begin(), edgeSequence.end(), back_inserter(chain.sequence));

        // If this was the last primary edge for the chain, we are done.
        if(positionInChain == chain.size() - 1) {
            break;
        }

        // Add assembled sequence between this marker graph primayr edge and the next in the chain.
        const vector<Base>& stepSequence = chain.stepSequences[positionInChain];
        copy(stepSequence.begin(), stepSequence.end(), back_inserter(chain.sequence));

    }
}



// This writes the details of sequence assembly for all Chains in the AssemblyGraph.
void AssemblyGraph::writeAssemblyDetails() const
{
    const AssemblyGraph& assemblyGraph = *this;

    // Opeb the csv file and write the header.
    ofstream csv("AssemblyDetails-" + to_string(componentId) + ".csv");
    csv << "Chain,Component,Bubble chain,Position in bubble chain,Index in bubble,"
        "Position in chain,Type,Marker graph edge id,"
        "Length,Sequence begin,Sequence end,Coverage,Common\n";

    // Loop over all bubble chains.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const BubbleChain& bubbleChain = assemblyGraph[e];

        // Loop over Bubbles in this BubbleChain.
        for(uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size(); positionInBubbleChain++) {
            const Bubble& bubble = bubbleChain[positionInBubbleChain];

            // Loop over Chains in this Bubble.
            for(uint64_t indexInBubble=0; indexInBubble<bubble.size(); indexInBubble++) {
                const Chain& chain = bubble[indexInBubble];
                SHASTA_ASSERT(chain.wasAssembled);
                SHASTA_ASSERT(chain.stepSequences.size() == chain.size() - 1);
                const string chainString = chainStringId(e, positionInBubbleChain, indexInBubble);

                // Loop over positions in this Chain.
                uint64_t positionInSequence = 0;
                for(uint64_t positionInChain=0; /* Check later */ ; positionInChain++) {

                    // Write one line to csv with information about the sequence
                    // contributed by this the marker graph primary edge.
                    {
                        const MarkerGraphEdgeId edgeId = chain[positionInChain];
                        const uint64_t coverage = assembler.markerGraph.edgeMarkerIntervals[edgeId].size();
                        const uint64_t edgeSequenceLength = assembler.markerGraph.edgeSequence[edgeId].size();
                        const uint64_t beginInSequence = positionInSequence;
                        const uint64_t endInSequence = positionInSequence + edgeSequenceLength;
                        csv << chainString << ",";
                        csv << componentId << ",";
                        csv << assemblyGraph[e].id << ",";
                        csv << positionInBubbleChain << ",";
                        csv << indexInBubble << ",";
                        csv << positionInChain << ",";
                        csv << "E,";
                        csv << edgeId << ",";
                        csv << edgeSequenceLength << ",";
                        csv << beginInSequence << ",";
                        csv << endInSequence << ",";
                        csv << coverage << ",";
                        csv << ",";
                        csv << "\n";
                        positionInSequence = endInSequence;
                    }


                    // If this was the last primary edge for the chain, we are done.
                    if(positionInChain == chain.size() - 1) {
                        SHASTA_ASSERT(positionInSequence == chain.sequence.size());
                        break;
                    }

                    // Write one line to csv with information about the sequence
                    // contributed by the assemby step between this marker graph primary edge
                    // and the next in the chain.
                    {
                        const MarkerGraphEdgeId edgeId = chain[positionInChain];
                        const MarkerGraphEdgeId nextEdgeId = chain[positionInChain + 1];
                        const uint64_t commonCount = assembler.countCommonOrientedReadsUnsafe(
                            edgeId, nextEdgeId);
                        const uint64_t stepSequenceLength = chain.stepSequences[positionInChain].size();
                        const uint64_t beginInSequence = positionInSequence;
                        const uint64_t endInSequence = positionInSequence + stepSequenceLength;
                        csv << chainString << ",";
                        csv << componentId << ",";
                        csv << assemblyGraph[e].id << ",";
                        csv << positionInBubbleChain << ",";
                        csv << indexInBubble << ",";
                        csv << ",";
                        csv << "S,";
                        csv << ",";
                        csv << stepSequenceLength << ",";
                        csv << beginInSequence << ",";
                        csv << endInSequence << ",";
                        csv << ",";
                        csv << commonCount << ",";
                        csv << "\n";
                        positionInSequence = endInSequence;
                    }

                }
            }
        }
    }
}



void AssemblyGraph::assembleChainsMultithreadedTheadFunction(uint64_t threadId)
{
    const uint64_t chainTerminalCommonThreshold = assembleChainsMultithreadedData.chainTerminalCommonThreshold;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all assembly steps assigned to this batch.
        for(uint64_t i=begin; i!=end; ++i) {
            const auto& assemblyStep = assembleChainsMultithreadedData.assemblySteps[i];
            runAssemblyStep(chainTerminalCommonThreshold, assemblyStep);
        }
    }
}



void AssemblyGraph::runAssemblyStep(
    uint64_t chainTerminalCommonThreshold,
    const AssemblyStep& assemblyStep)
{
    AssemblyGraph& assemblyGraph = *this;

    // Get the BubbleChain.
    BubbleChain& bubbleChain = assemblyGraph[assemblyStep.e];

    // Get the Bubble.
    Bubble& bubble = bubbleChain[assemblyStep.positionInBubbleChain];

    // Get the Chain.
    Chain& chain = bubble[assemblyStep.indexInBubble];
    SHASTA_ASSERT(chain.size() >= 2);

    // Do it.
    runAssemblyStep(chain, assemblyStep.positionInChain, chainTerminalCommonThreshold);
}



void AssemblyGraph::runAssemblyStep(
    Chain& chain,
    uint64_t positionInChain,
    uint64_t chainTerminalCommonThreshold)
{

    // Find the MarkerGraphEdgeIds for this local assembly.
    const MarkerGraphEdgeId edgeIdA = chain[positionInChain];
    const MarkerGraphEdgeId edgeIdB = chain[positionInChain + 1];

    // Suppress html output from LocalAssembly.
    ostream html(0);



    // Figure out if we should use the oriented reads on edgeIdA and edgeIdB.
    bool useA = true;
    bool useB = true;
    // For chains of length 2, we leave useA and useB set to true.
    // For the usual case of longer chains, there is more checking.
    if(chain.size() != 2) {

        // If we are at the beginning or end of the chain, we need to check
        // the number of common oriented reads.
        MarkerGraphEdgePairInfo info;
        if((positionInChain == 0) or (positionInChain == chain.size() - 2)) {
            SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeIdA, edgeIdB, info));
        }

        // If this is the first step of the Chain, we want to set useA to false
        // to avoid using reads that don't belong. But we only do it
        // if this leaves us with enough reads to assemble.
        if(positionInChain == 0) {
            if(info.common >= chainTerminalCommonThreshold) {
                useA = false;
            }
        }

        // If this is the last step of the Chain, we want to set useB to false
        // to avoid using reads that don't belong. But we only do it
        // if this leaves us with enough reads to assemble.
        else if(positionInChain == chain.size() - 2) {
            if(info.common >= chainTerminalCommonThreshold) {
                useB = false;
            }
        }
    }



    // Do the local assembly between these two MarkerGraphEdgeIds.
    try {
        LocalAssembly localAssembly(assembler, edgeIdA, edgeIdB, 0, html, useA, useB);
        localAssembly.getSecondarySequence(chain.stepSequences[positionInChain]);
    } catch (...) {
        std::lock_guard<std::mutex> lock(mutex);
        cout << "Error occurred in local assembly between marker graph edges " <<
            edgeIdA << " and " << edgeIdB << endl;
        throw;
    }
}



// Make a copy of an edge, truncating it at its end by removing the last MarkerGraphEdgeId.
// Return the target vertex of the newly created edge.
// The last bubble of the bubble chain of the given edge must be haploid.
// If the bubble chain consists of just a single haploid bubble with a chain of length 2,
// no new edge is created, and this simply returns the source vertex of the given edge.
AssemblyGraph::vertex_descriptor
    AssemblyGraph::cloneAndTruncateAtEnd(edge_descriptor ce)
{
    AssemblyGraph& cGraph = *this;
    const AssemblyGraphEdge& edge = cGraph[ce];
    const vertex_descriptor cv0 = source(ce, cGraph);
    const BubbleChain& bubbleChain = cGraph[ce];

    // Sanity checks.
    SHASTA_ASSERT(not bubbleChain.empty());
    SHASTA_ASSERT(bubbleChain.lastBubble().isHaploid());



    // Case where the bubble chain consists of a single bubble, which must be haploid,
    // that is, consist of a single chain.
    if(bubbleChain.size() == 1) {
        const Bubble& bubble = bubbleChain.lastBubble();
        SHASTA_ASSERT(bubble.isHaploid());
        const Chain& chain = bubble.front();
        SHASTA_ASSERT(chain.size() > 1);

        // If the Chain has length 2, we can't truncate it.
        // So we don't create a new edge, and instead just return cv0.
        // Detangling code will connect there, as prescribed by the tangle matrix.
        if(chain.size() == 2) {
            return cv0;
        }

        // Create the new edge, without adding it to the graph for now.
        AssemblyGraphEdge newEdge = edge;
        newEdge.id = nextEdgeId++;
        BubbleChain& newBubbleChain = newEdge;
        SHASTA_ASSERT(newBubbleChain.size() == 1);
        Bubble& newBubble = newBubbleChain.lastBubble();
        SHASTA_ASSERT(newBubble.isHaploid());
        Chain& newChain = newBubble.front();
        SHASTA_ASSERT(chain.size() > 2);
        newChain.pop_back();    // Remove the last MarkerGraphEdgeId.

        // Add it to the graph.
        // It will be dangling at its end.
        // Detangling code will later connect it s prescribed by the tangle matrix.
        const vertex_descriptor cv2 = createVertex(newBubbleChain.lastMarkerGraphEdgeId());
        add_edge(cv0, cv2, newEdge, cGraph);
        return cv2;
    }



    // Case where the bubble chain consists of more than one bubble.
    else {
        const Bubble& lastBubble = bubbleChain.lastBubble();
        SHASTA_ASSERT(lastBubble.isHaploid());
        const Chain& lastChain = lastBubble.front();
        SHASTA_ASSERT(lastChain.size() > 1);

        // Create the new edge, without adding it to the graph for now.
        AssemblyGraphEdge newEdge = edge;
        newEdge.id = nextEdgeId++;
        BubbleChain& newBubbleChain = newEdge;
        SHASTA_ASSERT(newBubbleChain.size() > 1);
        Bubble& newLastBubble = newBubbleChain.lastBubble();
        SHASTA_ASSERT(newLastBubble.isHaploid());
        Chain& newLastChain = newLastBubble.front();

        // If the last chain has length 2, just remove the last bubble from newBubbleChain.
        // Otherwise, remove the last MarkerGraphEdgeId from the lastChain.
        if(newLastChain.size() == 2) {
            newBubbleChain.pop_back();
        } else {
            newLastChain.pop_back();
        }

        // Add it to the graph.
        // It will be dangling at its end.
        // Detangling code will later connect it s prescribed by the tangle matrix.
        const vertex_descriptor cv2 = createVertex(newBubbleChain.lastMarkerGraphEdgeId());
        add_edge(cv0, cv2, newEdge, cGraph);
        return cv2;
    }

}





// Make a copy of an edge, truncating it at its beginning by removing the first MarkerGraphEdgeId.
// Return the source vertex of the newly created edge.
// The first bubble of the bubble chain of the given edge must be haploid.
// If the bubble chain consists of just a single haploid bubble with a chain of length 2,
// no new edge is created, and this simply returns the target vertex of the given edge.
AssemblyGraph::vertex_descriptor
    AssemblyGraph::cloneAndTruncateAtBeginning(edge_descriptor ce)
{
    AssemblyGraph& cGraph = *this;
    const AssemblyGraphEdge& edge = cGraph[ce];
    const vertex_descriptor cv1 = target(ce, cGraph);
    const BubbleChain& bubbleChain = cGraph[ce];

    // Sanity checks.
    SHASTA_ASSERT(not bubbleChain.empty());
    SHASTA_ASSERT(bubbleChain.firstBubble().isHaploid());



    // Case where the bubble chain consists of a single bubble, which must be haploid,
    // that is, consist of a single chain.
    if(bubbleChain.size() == 1) {
        const Bubble& bubble = bubbleChain.firstBubble();
        SHASTA_ASSERT(bubble.isHaploid());
        const Chain& chain = bubble.front();
        SHASTA_ASSERT(chain.size() > 1);

        // If the Chain has length 2, we can't truncate it.
        // So we don't create a new edge, and instead just return cv1.
        // Detangling code will connect there, as prescribed by the tangle matrix.
        if(chain.size() == 2) {
            return cv1;
        }

        // Create the new edge, without adding it to the graph for now.
        AssemblyGraphEdge newEdge = edge;
        newEdge.id = nextEdgeId++;
        BubbleChain& newBubbleChain = newEdge;
        SHASTA_ASSERT(newBubbleChain.size() == 1);
        Bubble& newBubble = newBubbleChain.firstBubble();
        SHASTA_ASSERT(newBubble.isHaploid());
        Chain& newChain = newBubble.front();
        SHASTA_ASSERT(chain.size() > 2);
        newChain.erase(newChain.begin());    // Remove the first MarkerGraphEdgeId.

        // Add it to the graph.
        // It will be dangling at its beginning.
        // Detangling code will later connect it s prescribed by the tangle matrix.
        const vertex_descriptor cv2 = createVertex(newBubbleChain.firstMarkerGraphEdgeId());
        add_edge(cv2, cv1, newEdge, cGraph);
        return cv2;
    }



    // Case where the bubble chain consists of more than one bubble.
    else {
        const Bubble& firstBubble = bubbleChain.firstBubble();
        SHASTA_ASSERT(firstBubble.isHaploid());
        const Chain& firstChain = firstBubble.front();
        SHASTA_ASSERT(firstChain.size() > 1);

        // Create the new edge, without adding it to the graph for now.
        AssemblyGraphEdge newEdge = edge;
        newEdge.id = nextEdgeId++;
        BubbleChain& newBubbleChain = newEdge;
        SHASTA_ASSERT(newBubbleChain.size() > 1);
        Bubble& newFirstBubble = newBubbleChain.firstBubble();
        SHASTA_ASSERT(newFirstBubble.isHaploid());
        Chain& newFirstChain = newFirstBubble.front();

        // If the last chain has length 2, just remove the first bubble from newBubbleChain.
        // Otherwise, remove the first MarkerGraphEdgeId from the lastChain.
        if(newFirstChain.size() == 2) {
            newBubbleChain.erase(newBubbleChain.begin());
        } else {
            newFirstChain.erase(newFirstChain.begin());
        }

        // Add it to the graph.
        // It will be dangling at its end.
        // Detangling code will later connect it s prescribed by the tangle matrix.
        const vertex_descriptor cv2 = createVertex(newBubbleChain.firstMarkerGraphEdgeId());
        add_edge(cv2, cv1, newEdge, cGraph);
        return cv2;
    }

}


// Create a new edge connecting the cv0 and cv1.
// The new edge will consist of a simple BubbleChain with a single
// haploid Bubble with a Chain of length 2.
AssemblyGraph::edge_descriptor AssemblyGraph::connect(vertex_descriptor cv0, vertex_descriptor cv1)
{
    AssemblyGraph& cGraph = *this;

    edge_descriptor ceNew;
    tie(ceNew, ignore) = add_edge(cv0, cv1, cGraph);
    AssemblyGraphEdge& newEdge = cGraph[ceNew];
    newEdge.id = nextEdgeId++;
    BubbleChain& newBubbleChain = newEdge;

    // The new BubbleChain consists of a single Bubble.
    newBubbleChain.resize(1);
    Bubble& bubble = newBubbleChain.front();

    // The new Bubble is haploid, that is, consists of a single Chain.
    bubble.resize(1);

    // The new Bubble consists of just the two MarkerGraphEdgeIds
    // corresponding to cv0 and cv1.
    Chain& chain = bubble.front();
    chain.push_back(cGraph[cv0].edgeId);
    chain.push_back(cGraph[cv1].edgeId);

    return ceNew;

}



void AssemblyGraph::save(const string& fileName) const
{
    ofstream file(fileName);
    boost::archive::binary_oarchive archive(file);
    archive << *this;
}



void AssemblyGraph::load(const string& fileName)
{
    ifstream file(fileName);
    boost::archive::binary_iarchive archive(file);
    archive >> *this;
}



// Optimize chains before assembly, to remove assembly steps with
// less that minCommon reads.
void AssemblyGraph::optimizeChains(
    bool debug,
    uint64_t minCommon,
    uint64_t k)
{
    AssemblyGraph& cGraph = *this;

    BGL_FORALL_EDGES(ce, cGraph, AssemblyGraph) {
        BubbleChain& bubbleChain = cGraph[ce];

        for(uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size(); positionInBubbleChain++) {
            Bubble& bubble = bubbleChain[positionInBubbleChain];
            const uint64_t ploidy = bubble.size();

            for(uint64_t indexInBubble=0; indexInBubble<ploidy; indexInBubble++) {
                Chain& chain = bubble[indexInBubble];
                SHASTA_ASSERT(chain.size() >= 2);

                if(debug) {
                    cout << "Optimizing chain " << chainStringId(ce, positionInBubbleChain, indexInBubble) << endl;
                }
                optimizeChain(debug, chain, minCommon, k);
            }
        }
    }

}



// Optimize a chain before assembly, to remove assembly steps with
// less that minCommon reads.
void AssemblyGraph::optimizeChain(
    bool debug,
    Chain& chain,
    uint64_t minCommon,
    uint64_t k)
{
    if(debug) {
        cout << "Optimizing a chain of length " << chain.size() << endl;
    }
    SHASTA_ASSERT(chain.size() >= 2);


    // A directed graph describing the initial and final chains.
    // Each vertex stores a MarkerGraphEdgeId.
    // Each edge stores the number of common oriented reads.
    class ChainGraphVertex {
    public:
        MarkerGraphEdgeId edgeId;
        uint64_t immediateDominator = invalid<uint64_t>;
    };
    class ChainGraphEdge {
    public:
        uint64_t commonCount;
        bool keep = false;
    };
    using ChainGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::bidirectionalS,
        ChainGraphVertex,
        ChainGraphEdge>;
    class ChainGraph : public ChainGraphBaseClass {
    public:
    };
    ChainGraph chainGraph;

    class PathInspector {
    public:
        PathInspector(ChainGraph& chainGraph, bool debug) : chainGraph(chainGraph), debug(debug) {}
        ChainGraph& chainGraph;
        bool debug;
        using Path = vector<ChainGraph::edge_descriptor>;
        Path bestPath;
        uint64_t bestPathMinCommonCount = 0;
        void operator()(const Path& path)
        {
            // Compute the minimum number of common oriented reads over edges of this path.
            uint64_t minCommonCount = invalid<uint64_t>;
            for(const ChainGraph::edge_descriptor e: path) {
                minCommonCount = min(minCommonCount, chainGraph[e].commonCount);
            }

            if(debug) {
                cout << "Path with minCommonCount " << minCommonCount << ":";
                for(const ChainGraph::edge_descriptor e: path) {
                    cout << " " << source(e, chainGraph);
                }
                cout << " " << target(path.back(), chainGraph) << "\n";
            }

            // A Path is better if it has a higher minCommonCount or
            // it has the same minCommonCount and is longer.
            //
            if( (minCommonCount >  bestPathMinCommonCount) or
                (minCommonCount == bestPathMinCommonCount and path.size() > bestPath.size())) {
                bestPath = path;
                bestPathMinCommonCount = minCommonCount;
            }
        }

    };

    // Construct the initial ChainGraph.

    // Add the vertices.
    // We are using vecS as the second template argument for ChainGraph,
    // so positions in the chain are also vertex descriptors in the ChainGraph.
    for(const MarkerGraphEdgeId edgeId: chain) {
        add_vertex({edgeId}, chainGraph);
    }

    // Add the edges that correspond to the initial Chain.
    for(uint64_t i1=1; i1<chain.size(); i1++) {
        const uint64_t i0 = i1 - 1;
        const MarkerGraphEdgeId edgeId0 = chainGraph[i0].edgeId;
        const MarkerGraphEdgeId edgeId1 = chainGraph[i1].edgeId;
        MarkerGraphEdgePairInfo info;
        SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(edgeId0, edgeId1, info));
        add_edge(i0, i1, {info.common}, chainGraph);
    }



    // Add edges that skip around any edges with less than minCommon common oriented reads.
    uint64_t totalAddedEdgesCount = 0;
    uint64_t totalRemovedEdgesCount = 0;
    for(uint64_t i1=1; i1<chain.size(); i1++) {
        const uint64_t i0 = i1 - 1;
        ChainGraph::edge_descriptor e;
        bool edgeWasFound = false;
        tie(e, edgeWasFound) = edge(i0, i1, chainGraph);
        SHASTA_ASSERT(edgeWasFound);

        // If this edge has enough common reads, don't do anything.
        if(chainGraph[e].commonCount >= minCommon) {
            continue;
        }

        if(debug) {
            cout << i0 << "->" << i1 << " " << chainGraph[i0].edgeId << "->" << chainGraph[i1].edgeId <<
                " has " << chainGraph[e].commonCount << " common oriented reads, adding edges to skip it." << endl;
        }

        // Loop over pairs of predecessors of v0 and successors of v1.
        uint64_t addedEdgesCount = 0;
        const uint64_t j0First = (k < i0) ? (i0 - k) : 0;
        const uint64_t j0Last = i0;
        const uint64_t j1First = i1;
        const uint64_t j1Last = min(i1 + k, chain.size() - 1);
        for(uint64_t j0=j0First; j0<=j0Last; j0++) {
            for(uint64_t j1=j1First; j1<=j1Last; j1++) {
                if(j0==i0 and j1 == i1) {
                    // We already have the edge between v0 and v1.
                    continue;
                }
                MarkerGraphEdgePairInfo info;
                SHASTA_ASSERT(assembler.analyzeMarkerGraphEdgePair(chainGraph[j0].edgeId, chainGraph[j1].edgeId, info));

                // If the number of common reads is better than for e, add this edge.
                if(info.common > chainGraph[e].commonCount) {
                    add_edge(j0, j1, {info.common}, chainGraph);
                    ++addedEdgesCount;
                    if(debug) {
                    cout << " Added " << j0 << "->" << j1 << " " << chainGraph[j0].edgeId << "->" << chainGraph[j1].edgeId <<
                        " with " << info.common << " common oriented reads." << endl;
                    }
                } else {
                    if(debug) {
                        cout << "Found " << j0 << "->" << j1 << " " << chainGraph[j0].edgeId << "->" << chainGraph[j1].edgeId <<
                            " with " << info.common << " common oriented reads." << endl;

                    }
                }
            }
        }
        totalAddedEdgesCount += addedEdgesCount;

        // If we added any edges skipping e, we can remove e.
        if(addedEdgesCount > 0) {
            if(debug) {
                cout << "Removed " << i0 << "->" << i1 << " " << chainGraph[i0].edgeId << "->" << chainGraph[i1].edgeId <<
                    " with " << chainGraph[e].commonCount << " common oriented reads." << endl;
            }
            // DON'T REMOVE THE EDGE. THIS IS NECESSARY TO MAKE SURE WE
            // STILL HAVE A PATH FROM THE ENTRANCE TO THE EXIT.
            // boost::remove_edge(e, chainGraph);
            // ++totalRemovedEdgesCount;
        } else {
            if(debug) {
                cout << "Did not find any suitable replacement edges." << endl;
            }
        }
    }


    // If we did not add or remove any edges, leave this Chain alone.
    if(totalAddedEdgesCount == 0) {
        SHASTA_ASSERT(totalRemovedEdgesCount == 0);
        if(debug) {
            cout << "No edges were added or removed, so this Chain will be left unchanged." << endl;
        }
        return;
    }

    if(debug) {
        cout << "This chain will be optimized." << endl;
    }



    // To find the optimized chain, we want to do path enumeration on the ChainGraph,
    // looking for paths that only use edges with large numbers of common oriented reads.
    // Specifically, we use as the new chain the path that maximizes the minimum
    // number of common oriented reads encountered on edges along the path.
    // For efficiency of the path enumeration, we compute a dominator tree
    // for the ChainGraph, with entrance at the beginning of the chain.
    // The unique path on that tree from the entrance to the exit
    // divides the graph in segments, and we can do path enumeration on one segment at a time.
    shasta::lengauer_tarjan_dominator_tree(chainGraph, 0,
        boost::get(&ChainGraphVertex::immediateDominator, chainGraph));

    // The unique path on the dominator tree from the entrance to the exit.
    vector<ChainGraph::vertex_descriptor> dominatorTreePath;
    ChainGraph::vertex_descriptor v = chain.size() - 1;
    while(true) {
        dominatorTreePath.push_back(v);
        if(v == 0) {
            break;
        }
        v = chainGraph[v].immediateDominator;
        if(v == invalid<uint64_t>) {
            cout << "Assertion failure at " << v << endl;
        }
        SHASTA_ASSERT(v != invalid<uint64_t>);
    }
    if(debug) {
        cout << "Dominator tree path length " << dominatorTreePath.size() << endl;
    }
    reverse(dominatorTreePath.begin(), dominatorTreePath.end());

    if(false) {
        cout << "Dominator tree path:" << endl;
        for(uint64_t i=0; i<dominatorTreePath.size(); i++) {
            const uint64_t v = dominatorTreePath[i];
            cout << i << "," << v << "," << chainGraph[v].edgeId << "\n";
        }
    }



    // The dominator tree path divides the graph in segments,
    // and we can do path enumeration on one segment at a time.
    // For each segment we find the best path and mark the edges on that
    // best path as to be kept in the final chain.
    for(uint64_t i1=1; i1<dominatorTreePath.size(); i1++) {
        const uint64_t i0 = i1 - 1;
        const ChainGraph::vertex_descriptor v0 = dominatorTreePath[i0];
        const ChainGraph::vertex_descriptor v1 = dominatorTreePath[i1];

        // Fast handling of the most common case.
        if(v1 == v0+1 and out_degree(v0, chainGraph)==1 and in_degree(v1, chainGraph)==1) {
            ChainGraph::edge_descriptor e;
            bool edgeWasFound = true;
            tie(e, edgeWasFound) = edge(v0, v1, chainGraph);
            if(edgeWasFound) {
                chainGraph[e].keep = true;
                continue;
            }
        }

        // If getting here, we have to do path enumeration.
        if(debug) {
            cout << "Starting path enumeration between " << v0 << " " << v1 << endl;
        }

        // Enumerate paths starting at v0 and ending at v1.
        PathInspector pathInspector(chainGraph, debug);
        enumeratePathsBetween(chainGraph, v0, v1, pathInspector);

        if(debug) {
            if(debug) {
                cout << "The best path has minCommonCount " << pathInspector.bestPathMinCommonCount << ":";
                for(const ChainGraph::edge_descriptor e: pathInspector.bestPath) {
                    cout << " " << source(e, chainGraph);
                }
                cout << " " << target(pathInspector.bestPath.back(), chainGraph) << "\n";
            }
        }

        // Mark as to be kept all edges on the best path.
        for(const ChainGraph::edge_descriptor e: pathInspector.bestPath) {
            chainGraph[e].keep = true;
        }
    }


    // Remove all edges not marked to be kept.
    vector<ChainGraph::edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, chainGraph, ChainGraph) {
        if(not chainGraph[e].keep) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const ChainGraph::edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, chainGraph);
    }

    // The remaining edges should form a path in the ChainGraph
    // which defines the optimized Chain.
    SHASTA_ASSERT(in_degree(0, chainGraph) == 0);
    SHASTA_ASSERT(out_degree(0, chainGraph) == 1);
    SHASTA_ASSERT(in_degree(chain.size()-1, chainGraph) == 1);
    SHASTA_ASSERT(out_degree(chain.size()-1, chainGraph) == 0);
    for(uint64_t i=1; i<chain.size()-1; i++) {
        const uint64_t inDegree = in_degree(i, chainGraph);
        const uint64_t outDegree = out_degree(i, chainGraph);
        SHASTA_ASSERT(
            (inDegree==1 and outDegree==1) or   // In the new chain.
            (inDegree==0 and outDegree==0)      // Now isolated.
            );
    }

    // Find the path from the entrance to the exit in the update ChainGraph.
    vector<uint64_t> newPath;
    v = 0;
    while(true) {
        newPath.push_back(v);
        if(v == chain.size()-1) {
            break;
        }

        // Move forward.
        SHASTA_ASSERT(out_degree(v, chainGraph) == 1);
        ChainGraph::out_edge_iterator it;
        tie(it, ignore) = out_edges(v, chainGraph);
        const ChainGraph::edge_descriptor e = *it;
        v = target(e, chainGraph);
    }

    // Sanity check that the path is moving forward.
    for(uint64_t i=1; i<newPath.size(); i++) {
        SHASTA_ASSERT(newPath[i] > newPath[i-1]);
    }

    // Construct the new Chain.
    chain.clear();
    chain.sequence.clear();
    for(const uint64_t v: newPath) {
        chain.push_back(chainGraph[v].edgeId);
    }

}



bool AssemblyGraph::removeSelfComplementaryEdges()
{
    AssemblyGraph& cGraph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(ce, cGraph, AssemblyGraph) {
        const vertex_descriptor v0 = source(ce, cGraph);
        const vertex_descriptor v1 = target(ce, cGraph);
        const MarkerGraphEdgeId edgeId0 = cGraph[v0].edgeId;
        const MarkerGraphEdgeId edgeId1 = cGraph[v1].edgeId;

        if(assembler.markerGraph.reverseComplementEdge[edgeId0] == edgeId1) {
            SHASTA_ASSERT(assembler.markerGraph.reverseComplementEdge[edgeId1] == edgeId0);
            edgesToBeRemoved.push_back(ce);
        }
    }

    for(const edge_descriptor ce: edgesToBeRemoved) {
        boost::remove_edge(ce, cGraph);
    }

    return not edgesToBeRemoved.empty();
}



// Split terminal haploid bubbles out of bubble chains, to facilitate detangling.
void AssemblyGraph::splitTerminalHaploidBubbles()
{
    AssemblyGraph& cGraph = *this;

    vector<edge_descriptor> allEdges;
    BGL_FORALL_EDGES(ce, cGraph, AssemblyGraph) {
        allEdges.push_back(ce);
    }

    for(const edge_descriptor e: allEdges) {
        splitTerminalHaploidBubbles(e);
    }
}



void AssemblyGraph::splitTerminalHaploidBubbles(edge_descriptor ce)
{
    AssemblyGraph& cGraph = *this;
    BubbleChain& bubbleChain = cGraph[ce];

    // Skip trivial bubble chains consisting of a single bubble.
    if(bubbleChain.size() < 2) {
        return;
    }

    // Access the first and last bubble in the bubble chain.
    // We already checked that the bubble chain has at least two bubbles,
    // so these two are distinct.
    const Bubble& firstBubble = bubbleChain.front();
    const Bubble& lastBubble = bubbleChain.back();

    // Skip bubble chains consisting of two haploid bubbles.
    // After compress() is called, there should be none of these.
    if(bubbleChain.size() == 2 and firstBubble.isHaploid() and lastBubble.isHaploid()) {
        return;
    }

    // Figure out if we need to split the first or last bubble, or both.
    bool splitFirstBubble = false;
    bool splitLastBubble = false;
    if(firstBubble.isHaploid()) {
        splitFirstBubble = true;
    }
    if(lastBubble.isHaploid()) {
        splitLastBubble = true;
    }
    if(splitFirstBubble and splitLastBubble) {
        SHASTA_ASSERT(bubbleChain.size() > 2);
    }

    // If there is nothing to do, we are done.
    if(not (splitFirstBubble or splitLastBubble)) {
        return;
    }

    // The source and target vertices of the edge we are splitting.
    const vertex_descriptor cv0 = source(ce, cGraph);
    const vertex_descriptor cv1 = target(ce, cGraph);
    vertex_descriptor cv2 = null_vertex();
    vertex_descriptor cv3 = null_vertex();



    // Create a new edge with just the first bubble, if necessary.
    if(splitFirstBubble) {

        // Get the target vertex for the new edge.
        const Chain& firstChain = firstBubble.front();
        const MarkerGraphEdgeId markerGraphEdgeId2 = firstChain.back();
        cv2 = createVertex(markerGraphEdgeId2);

        // Add the new edge.
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(cv0, cv2, cGraph);
        AssemblyGraphEdge& newEdge = cGraph[eNew];
        newEdge.id = nextEdgeId++;

        // Copy the first bubble to the new edge.
        newEdge.push_back(firstBubble);

    }



    // Create a new edge with just the last bubble, if necessary.
    if(splitLastBubble) {

        // Get the source vertex for the new edge.
        const Chain& lastChain = lastBubble.front();
        const MarkerGraphEdgeId markerGraphEdgeId3 = lastChain.front();
        cv3 = createVertex(markerGraphEdgeId3);

        // Add the new edge.
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(cv3, cv1, cGraph);
        AssemblyGraphEdge& newEdge = cGraph[eNew];
        newEdge.id = nextEdgeId++;

        // Copy the last bubble to the new edge.
        newEdge.push_back(lastBubble);

    }



    // Create a new edge for the rest of the bubble chain.
    edge_descriptor eNew;
    tie(eNew, ignore) = add_edge(
        splitFirstBubble ? cv2 : cv0,
        splitLastBubble ? cv3 : cv1,
        cGraph);
        AssemblyGraphEdge& newEdge = cGraph[eNew];
        newEdge.id = nextEdgeId++;

    // Copy the rest of the bubble chain to the new edge.
    auto it0 = bubbleChain.begin();
    auto it1 = bubbleChain.end();
    if(splitFirstBubble) {
        ++it0;
    }
    if(splitLastBubble) {
        --it1;
    }
    copy(it0, it1, back_inserter(newEdge));


    // Now we can remove the old BubbleChain we just split.
    boost::remove_edge(ce, cGraph);

}



// Bubble cleanup (all bubbles), with the purpose of eliminating most bubbles caused by errors.
uint64_t AssemblyGraph::cleanupBubbles(
    bool debug,
    uint64_t maxOffset,
    uint64_t chainTerminalCommonThreshold,
    uint64_t threadCount)
{
    AssemblyGraph& graph = *this;
    performanceLog << timestamp << "AssemblyGraph::cleanupBubbles begins." << endl;



    // First, assemble sequence for all the chains of diploid bubbles with a small offset.
    clearAllShouldBeAssembledFlags();
    BGL_FORALL_EDGES(e, graph, AssemblyGraph) {
        BubbleChain& bubbleChain = graph[e];
        for(Bubble& bubble: bubbleChain) {

            // If this bubble is not diploid, skip it.
            if(bubble.size() != 2) {
                continue;
            }

            // The bubble is diploid. Compute its maxOffset.
            uint64_t averageOffset;
            uint64_t minOffset;
            uint64_t bubbleMaxOffset;
            const uint64_t offsetWasComputed = bubbleOffsetNoException(
                bubble, averageOffset, minOffset, bubbleMaxOffset);

            // If the offset is large or could not be computed, we don't need to
            // assemble this bubble.
            if((not offsetWasComputed) or bubbleMaxOffset>maxOffset) {
                continue;
            }

            // We need to assemble the Chains of this bubble.
            for(Chain& chain: bubble) {
                chain.shouldBeAssembled = true;
            }
        }
    }
    assembleChainsMultithreaded(chainTerminalCommonThreshold, threadCount);
    performanceLog << timestamp << "Sequence assembly for AssemblyGraph::cleanupBubbles ends." << endl;



    uint64_t removedCount = 0;
    BGL_FORALL_EDGES(ce, graph, AssemblyGraph) {
        removedCount += cleanupBubbles(debug, ce, maxOffset, chainTerminalCommonThreshold);
    }

    performanceLog << timestamp << "AssemblyGraph::cleanupBubbles ends." << endl;
    return removedCount;
}



// Bubble cleanup for a bubble chain, with the purpose of eliminating most bubbles caused by errors.
uint64_t AssemblyGraph::cleanupBubbles(bool debug, edge_descriptor ce,
    uint64_t maxOffset, uint64_t chainTerminalCommonThreshold)
{
    AssemblyGraph& cGraph = *this;
    BubbleChain& bubbleChain = cGraph[ce];
    BubbleChain newBubbleChain;

    uint64_t removedCount = 0;
    for(uint64_t positionInBubbleChain=0; positionInBubbleChain<bubbleChain.size(); positionInBubbleChain++) {
        Bubble& bubble = bubbleChain[positionInBubbleChain];

        if(debug) {
            cout << "cleanupBubbles working on Bubble " << bubbleStringId(ce, positionInBubbleChain) <<
                " ploidy " << bubble.size() << endl;
            cout << "Entrance " << bubble.front().front() << ", exit " << bubble.front().back() << endl;
        }

        bool keepBubble = false;

        if(bubble.isHaploid()) {

            // The bubble is haploid. Keep it.
            keepBubble = true;

            if(debug) {
                cout << "Keeping this bubble because it is haploid." << endl;
            }

        } else {

            // The bubble is not haploid. Compute its maxOffset.
            uint64_t averageOffset;
            uint64_t minOffset;
            uint64_t bubbleMaxOffset;
            const bool offsetWasComputed = bubbleOffsetNoException(bubble, averageOffset, minOffset, bubbleMaxOffset);

            if((not offsetWasComputed) or bubbleMaxOffset>maxOffset) {

                // The bubble is not haploid but the offset is large. Keep it.
                keepBubble = true;

                if(debug) {
                    cout << "Keeping this bubble because it is not haploid but its offset is large." << endl;
                }

            } else {

                // The bubble is not haploid and has a small offset.

                if(bubble.size() > 2) {

                    // The bubble has a small offset and ploidy greater than 2. Remove it.
                    keepBubble = false;

                    if(debug) {
                        cout << "Removing this bubble because it has a small offset and ploidy greater than 2." << endl;
                    }

                } else {

                    // The bubble has a small offset and ploidy 2.
                    // Check that we assembled the sequence of its two sides.
                    for(Chain& chain: bubble) {
                        SHASTA_ASSERT(chain.wasAssembled);
                    }

                    if(debug) {
                        for(uint64_t indexInBubble=0; indexInBubble<2; indexInBubble++) {
                            const auto& sequence = bubble[indexInBubble].sequence;
                            cout << ">" << chainStringId(ce, positionInBubbleChain, indexInBubble) <<
                                " " << sequence.size() << "\n";
                            copy(sequence.begin(), sequence.end(), ostream_iterator<shasta::Base>(cout));
                            cout << "\n";
                        }
                    }
                    if(bubble[0].sequence == bubble[1].sequence) {
                        keepBubble = false;
                        if(debug) {
                            cout << "The two sides have identical sequence." << endl;
                        }
                    } else {

                        // Figure out if they differ by a copy number of short periodicity.
                        const uint64_t period = isCopyNumberDifference(bubble[0].sequence, bubble[1].sequence, 4);
                        if(debug) {
                            cout << "Period " << period << "\n";
                        }
                        keepBubble = (period == 0);
                    }
                }
            }


        }

        if(keepBubble) {
            newBubbleChain.push_back(bubble);
            if(debug) {
                cout << "Kept this bubble." << endl;
            }
        } else {
            // Remove the bubble and replace it with a haploid bubble
            // consisting of only the terminal MarkerGraphEdgeIds.
            Chain newChain;
            newChain.push_back(bubble.front().front());
            newChain.push_back(bubble.front().back());
            Bubble newBubble;
            newBubble.push_back(newChain);
            newBubbleChain.push_back(newBubble);
            ++removedCount;
            if(debug) {
                cout << "Removed this bubble." << endl;
            }
        }
    }

    bubbleChain.swap(newBubbleChain);
    return removedCount;
}



// This finds squares of the form:
// A->B
// A->B'
// B->A'
// B'->A'
// where a prime sign indicates reverse complementing.
// It then one of two pairs of self-complementary edges:
//     A->B  and B'->A'
// or
//     A->B' and B->A'
// The pair to be removed is selected in such a way that its removal
// does not introduce any dead ends.
// The code uses the following names:
// A0 = A
// A1 = A'
// B0 = B
// B1 = B'
void AssemblyGraph::removeSelfComplementarySquares()
{
    AssemblyGraph& cGraph = *this;
    const bool debug = true;

    vector< pair<edge_descriptor, vertex_descriptor> > outEdgesA0;


    // Do this iteratively.
    while(true) {


        // Loop over all possible choices for A0.
        bool done = false;
        BGL_FORALL_VERTICES(A0, cGraph, AssemblyGraph) {

            // Gather the children of A.
            outEdgesA0.clear();
            BGL_FORALL_OUTEDGES(A0, ce, cGraph, AssemblyGraph) {
                outEdgesA0.push_back({ce, target(ce, cGraph)});
            }

            // Look for a reverse complementary pair (B0, B1)
            // with edges B0->A1 and B1->A1.
            for(uint64_t i1=0; i1<outEdgesA0.size(); i1++) {
                const vertex_descriptor B1 = outEdgesA0[i1].second;
                const uint64_t edgeIdB1 = cGraph[B1].edgeId;
                const uint64_t edgeIdB0 = assembler.markerGraph.reverseComplementEdge[edgeIdB1];
                for(uint64_t i0=0; i0<i1; i0++) {
                    const vertex_descriptor B0 = outEdgesA0[i0].second;
                    if(cGraph[B0].edgeId == edgeIdB0) {

                        // We found it.

                        // Look for the edges B0->A1 and B1->A1.
                        const uint64_t edgeIdA0 = cGraph[A0].edgeId;
                        const uint64_t edgeIdA1 = assembler.markerGraph.reverseComplementEdge[edgeIdA0];

                        edge_descriptor B0A1;
                        vertex_descriptor A10 = null_vertex();
                        BGL_FORALL_OUTEDGES(B0, ce, cGraph, AssemblyGraph) {
                            const vertex_descriptor v = target(ce, cGraph);
                            if(cGraph[v].edgeId == edgeIdA1) {
                                B0A1 = ce;
                                A10 = v;
                                break;
                            }
                        }
                        if(A10 == null_vertex()) {
                            continue;
                        }

                        edge_descriptor B1A1;
                        vertex_descriptor A11 = null_vertex();
                        BGL_FORALL_OUTEDGES(B1, ce, cGraph, AssemblyGraph) {
                            const vertex_descriptor v = target(ce, cGraph);
                            if(cGraph[v].edgeId == edgeIdA1) {
                                B1A1 = ce;
                                A11 = v;
                                break;
                            }
                        }
                        if(A11 == null_vertex()) {
                            continue;
                        }

                        if(A10 != A11) {
                            continue;
                        }
                        const vertex_descriptor A1 = A10;

                        // We found a self-complementary square.
                        const edge_descriptor A0B0 = outEdgesA0[i0].first;
                        const edge_descriptor A0B1 = outEdgesA0[i1].first;

                        if(debug) {
                            cout << "Found a self-complementary square:\n" <<
                                cGraph[A0].edgeId << " " <<
                                cGraph[B0].edgeId << " " <<
                                cGraph[B1].edgeId << " " <<
                                cGraph[A1].edgeId << "\n" <<
                                bubbleChainStringId(A0B0) << " " <<
                                bubbleChainStringId(A0B1) << " " <<
                                bubbleChainStringId(B0A1) << " " <<
                                bubbleChainStringId(B1A1) << "\n";
                        }

                        // Remove two of the edges in the square,
                        // making sure to not introduce dead ends.
                        if(out_degree(A0, cGraph) > 1 and in_degree(A1, cGraph) > 1) {
                            if(in_degree (B0, cGraph) > 1 and out_degree(B1, cGraph) > 1) {
                                boost::remove_edge(A0B0, cGraph);
                                boost::remove_edge(B1A1, cGraph);
                                done = true;
                            } else if(in_degree(B1, cGraph) > 1 and out_degree(B0, cGraph) > 1) {
                                boost::remove_edge(A0B1, cGraph);
                                boost::remove_edge(B0A1, cGraph);
                                done = true;
                            }
                        }

                        if(done) {
                            break;
                        }
                    }

                    if(done) {
                        break;
                    }

                }
                if(done) {
                    break;
                }
            }
            if(done) {
                break;
            }
        }

        // If nothing happened, stop the outer iteration.
        if(not done) {
            break;
        }
    }
}