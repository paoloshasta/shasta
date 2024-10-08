#ifndef SHASTA_LOCAL_MARKER_GRAPH1_HPP
#define SHASTA_LOCAL_MARKER_GRAPH1_HPP

// Shasta.
#include "shastaTypes.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "iosfwd.hpp"
#include <map>
#include "string.hpp"
#include "vector.hpp"

namespace shasta {

    class LocalMarkerGraph1Vertex;
    class LocalMarkerGraph1Edge;
    class LocalMarkerGraph1;
    using LocalMarkerGraph1BaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        LocalMarkerGraph1Vertex,
        LocalMarkerGraph1Edge
        >;

    class CompressedMarker;
    class MarkerGraph;
    namespace MemoryMapped {
        template<class T, class Int> class VectorOfVectors;
    }
}


class shasta::LocalMarkerGraph1Vertex {
public:

    // The id of the corresponding marker graph vertex.
    MarkerGraphVertexId vertexId;

    // The distance from the start vertex.
    uint64_t distance;

    LocalMarkerGraph1Vertex(
        MarkerGraphVertexId vertexId,
        uint64_t distance) :
        vertexId(vertexId),
        distance(distance)
    {
    }

};



class shasta::LocalMarkerGraph1Edge {
public:

    // The id of the corresponding marker graph edge.
    MarkerGraphEdgeId edgeId;

    LocalMarkerGraph1Edge(MarkerGraphEdgeId edgeId) :
        edgeId(edgeId)
    {
    }

};



class shasta::LocalMarkerGraph1 :
    public LocalMarkerGraph1BaseClass {
public:

    LocalMarkerGraph1(
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph&,
        MarkerGraphVertexId,
        uint64_t maxDistance,
        uint64_t minVertexCoverage,
        uint64_t minEdgeCoverage
    );

    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    const MarkerGraph& markerGraph;
    uint64_t maxDistance;

    std::map<MarkerGraphVertexId, vertex_descriptor> vertexMap;
    std::map<MarkerGraphEdgeId, edge_descriptor> edgeMap;
    vertex_descriptor addVertex(MarkerGraphVertexId, uint64_t distance);

    void writeGfa(const string& fileName) const;
    void writeHtml0(
        ostream&,
        uint64_t sizePixels,
        uint64_t quality,
        double timeout,
        bool useSvg) const;
    void writeHtml1(
        ostream&,
        uint64_t sizePixels,
        double thicknessScaling,
        uint64_t quality,
        double edgeResolution,
        const string& coloring,
        uint64_t redCoverage,
        uint64_t greenCoverage,
        MarkerGraphEdgeId readFollowingStartEdgeId,
        int64_t firstMarkerOffset,
        int64_t lastMarkerOffset,
        bool showLabels,
        double timeout) const;

    void pruneLowCoverageLeaves(uint64_t maxPruneCoverage);
private:
    void pruneLowCoverageForwardLeaves(uint64_t maxPruneCoverage);
    void pruneLowCoverageBackwardLeaves(uint64_t maxPruneCoverage);

public:

    void removeLongLowCoverageChains(
        uint64_t maxChainCoverage,
        uint64_t minLength);
private:
    void findLowCoverageChains(
        uint64_t maxChainCoverage,
        vector< vector<vertex_descriptor> >&
        ) const;

};

#endif
