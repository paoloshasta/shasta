#pragma once

// Shasta.
#include "Kmer.hpp"
#include "invalid.hpp"
#include "MarkerInterval.hpp"
#include "MarkerKmers.hpp"
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"

// Boost libraries.
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/property_tree/ptree_fwd.hpp>

// Standard library.
#include "cstdint.hpp"
#include "memory.hpp"
#include "span.hpp"



namespace shasta {

    class Base;
    class CompressedMarker;
    class MarkerGraph;
    class MarkerInterval;
    class Reads;
    struct VariantPositionContext;


    // The main input to mode 3 assembly is a set of anchors.
    // Each anchor consists of a span of AnchorMarkerInterval, with the following requirements:
    // - All AnchorMarkerInterval correspond to exactly the same sequence in the corresponding oriented reads, and:
    //      * Those portions of the oriented reads are believed to be aligned.
    //      * They appear in a low number of copies in the genome being sequenced.
    // - There are no duplicate oriented reads in an anchor.
    // - The anchor coverage (number of oriented reads) is in [minPrimaryCoverage, maxPrimaryCoverage].

    namespace mode3 {

        using AnchorId = uint64_t;
        class Anchor;
        class AnchorMarkerInterval;
        class Anchors;
        class AnchorInfo;
        class AnchorPairInfo;

        using AnchorBaseClass = span<const AnchorMarkerInterval>;

        string anchorIdToString(AnchorId);
        AnchorId anchorIdFromString(const string&);
    }
}



// The second ordinal of the marker interval is not stored.
// It can be obtained by adding to ordinal0 the value returned
// by Anchors::ordinalOffset(AnchorId).
// This value is the same for all marker intervals of an Anchor, by construction.
// Currently this is also the same value for all Anchors and equal to 1,
// but this could change.
class shasta::mode3::AnchorMarkerInterval {
public:
    OrientedReadId orientedReadId;
    uint32_t ordinal0;
    uint32_t positionInJourney = invalid<uint32_t>;

    AnchorMarkerInterval() {}

    AnchorMarkerInterval(
        OrientedReadId orientedReadId,
        uint32_t ordinal0) :
        orientedReadId(orientedReadId),
        ordinal0(ordinal0)
    {}
};



class shasta::mode3::AnchorInfo {
public:
    uint32_t ordinalOffset = invalid<uint32_t>;
    uint32_t componentId = invalid<uint32_t>;
    uint64_t localAnchorIdInComponent = invalid<uint64_t>;
};



// An Anchor is a set of AnchorMarkerIntervals.
class shasta::mode3::Anchor : public AnchorBaseClass {
public:

    Anchor(const AnchorBaseClass& s) : AnchorBaseClass(s) {}

    void check() const;

    uint64_t coverage() const
    {
        return size();
    }

    // Return the number of common oriented reads with another Anchor.
    uint64_t countCommon(const Anchor& that, bool ignoreNegativeOffsets = false) const;
};



class shasta::mode3::Anchors :
    public MultithreadedObject<Anchors>,
    public MappedMemoryOwner {
public:

    // This constructor creates the Anchors from marker graph edges.
    Anchors(
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph&,
        uint64_t minPrimaryCoverage,
        uint64_t maxPrimaryCoverage,
        uint64_t threadCount);

    // This constructor creates the Anchors from marker k-mers.
    Anchors(
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        shared_ptr<MarkerKmers>,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage,
        uint64_t threadCount);

    // This constructor creates the Anchors from a json file.
    Anchors(
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const vector<string>& jsonFileNames,
        uint64_t minPrimaryCoverage,
        uint64_t maxPrimaryCoverage,
        uint64_t threadCount);

    // This constructor creates the Anchors from heterozygous sites (variant clustering).
    Anchors(
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const std::vector<uint64_t>& clusterRepresentatives,
        DisjointSets& disjointSets,
        const MemoryMapped::Vector<std::pair<OrientedReadId, uint32_t>>& positionPairs,
        const MemoryMapped::Vector<uint8_t>& positionPairAlleles,
        const MemoryMapped::Vector<VariantPositionContext>& positionPairContexts,
        uint64_t minClusterCoverage,
        uint64_t minAlleleCoverage,
        uint64_t threadCount);

    // This constructor access existing Anchors.
    Anchors(
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers);

    Anchor operator[](AnchorId) const;
    uint64_t size() const;


    // This returns the sequence between the midpoint of the first marker of the
    // anchor and the midpoint of the second marker.
    // When anchorCreationMethod is "FromMarkerKmers", this sequence is empty.
    span<const Base> anchorSequence(AnchorId) const;

    // This returns the sequence between the beginning of the first marker of the
    // anchor and the end of the second marker.
    vector<Base> anchorExtendedSequence(AnchorId) const;

    // Return the number of common oriented reads between two Anchors.
    uint64_t countCommon(AnchorId, AnchorId, bool ignoreNegativeOffsets = false) const;

    // Analyze the oriented read composition of two anchors.
    void analyzeAnchorPair(AnchorId, AnchorId, AnchorPairInfo&) const;
    void writeHtml(AnchorId, AnchorId, AnchorPairInfo&, ostream&) const;

    // Return true if the second Anchor is adjacent to the first one.
    // For precise definition see the code.
    bool areAdjacentAnchors(AnchorId, AnchorId) const;

    // The offset to be added to ordinal0 of an Anchor to obtain ordinal1.
    // * When constructing anchors from the marker graph, this is the same for all Anchors
    //   and always equal to 1.
    // * When reading the anchors from a json file, each anchor can have a different value.
    uint32_t ordinalOffset(AnchorId anchorId) const
    {
        return anchorInfos[anchorId].ordinalOffset;
    }

    void writeCoverageHistogram() const;

private:
    MemoryMapped::VectorOfVectors<AnchorMarkerInterval, uint64_t> anchorMarkerIntervals;

public:

    // Get the first ordinal for the AnchorMarkerInterval corresponding to a
    // given AnchorId and OrientedReadId.
    // This asserts if the given AnchorId does not contain an AnchorMarkerInterval
    // for the requested OrientedReadId.
    uint32_t getFirstOrdinal(AnchorId, OrientedReadId) const;

    const Reads& reads;
    uint64_t k;
    uint64_t kHalf;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;

    // The sequences of the anchors.
    // We assume that marker length k is even, and anchor sequence
    // starts at the midpoint of the first marker of the anchor
    // and ends at the midpoint of the second marker of the anchor.
    // This means:
    // - If the Anchor ordinal are consecutive (as it happens
    //   when getting Anchors from marker graph edges),
    //   its sequence is guaranteed to have at least one base.
    // - If the Anchor ordinal are identical (as it may happen
    //   in a future alignment free formulation), its sequence
    //   is empty.
    MemoryMapped::VectorOfVectors<Base, uint64_t> anchorSequences;


    // The journey of each oriented read is the sequence of AnchorIds
    // encountered by the oriented read.
    MemoryMapped::VectorOfVectors<AnchorId, uint64_t> journeys;
    void computeJourneys(uint64_t threadCount);
    void writeJourneys() const;
    void writeAnchorGapsByRead() const;
private:
    void computeJourneysThreadFunction1(uint64_t threadId);
    void computeJourneysThreadFunction2(uint64_t threadId);
    void computeJourneysThreadFunction12(uint64_t pass);
    void computeJourneysThreadFunction3(uint64_t threadId);
    void computeJourneysThreadFunction4(uint64_t threadId);

    // Temporary storage of journeys with ordinals.
    MemoryMapped::VectorOfVectors<pair<uint64_t, uint32_t>, uint64_t> journeysWithOrdinals;

    void check() const;

public:

    // For a given AnchorId, follow the read journeys forward/backward by one step.
    // Return a vector of the AnchorIds reached in this way.
    // The count vector is the number of oriented reads each of the AnchorIds.
    void findChildren(
        AnchorId,
        vector<AnchorId>&,
        vector<uint64_t>& count,
        uint64_t minCoverage = 0) const;
    void findParents(
        AnchorId,
        vector<AnchorId>&,
        vector<uint64_t>& count,
        uint64_t minCoverage = 0) const;


    // In addition to the marker intervals, we also store an AnchorInfo for each Anchor.
    MemoryMapped::Vector<AnchorInfo> anchorInfos;
public:
        void storeAnchorInfo(
            AnchorId anchorId,
            uint32_t componentId,
            uint64_t localAnchorIdInComponent)
        {
            AnchorInfo& anchorInfo = anchorInfos[anchorId];
            anchorInfo.componentId =  componentId;
            anchorInfo.localAnchorIdInComponent =  localAnchorIdInComponent;
        }
        uint64_t getComponent(AnchorId anchorId) const
        {
            return anchorInfos[anchorId].componentId;
        }
        uint64_t getLocalAnchorIdInComponent(AnchorId anchorId) const
        {
            return anchorInfos[anchorId].localAnchorIdInComponent;
        }


        // Read following.
        void followOrientedReads(
            AnchorId,
            uint64_t direction,                         // 0 = forward, 1 = backward
            uint64_t minCommonCount,
            double minJaccard,
            double minCorrectedJaccard,
            vector< pair<AnchorId, AnchorPairInfo> >&
            ) const;

private:



    // Data and functions used when constructing the Anchors from the MarkerGraph.
    class ConstructFromMarkerGraphData {
    public:
        uint64_t minPrimaryCoverage;
        uint64_t maxPrimaryCoverage;

        const MarkerGraph* markerGraphPointer;

        // The marker intervals of the anchors found by each thread.
        class ThreadMarkerInterval {
        public:
            OrientedReadId orientedReadId;
            uint32_t ordinal0;
        };
        vector< shared_ptr< MemoryMapped::VectorOfVectors<ThreadMarkerInterval, uint64_t> > > threadMarkerIntervals;

        // The corresponding sequences
        vector< shared_ptr< MemoryMapped::VectorOfVectors<Base, uint64_t> > > threadSequences;
    };
    ConstructFromMarkerGraphData constructFromMarkerGraphData;
    void constructFromMarkerGraphThreadFunction(uint64_t threadId);


    // Data and functions used when constructing the Anchors from marker k-mers.
    // New code that gets the Kmers from MarkerKmers.
    class ConstructFromMarkerKmersData {
    public:
        uint64_t minAnchorCoverage;
        uint64_t maxAnchorCoverage;

        shared_ptr<MarkerKmers> markerKmers;

        // The MarkerInfo objects for the candidate anchors found by each thread.
        using MarkerInfo = MarkerKmers::MarkerInfo;
        vector< shared_ptr<MemoryMapped::VectorOfVectors<MarkerInfo, uint64_t> > > threadAnchors;
    };
    ConstructFromMarkerKmersData constructFromMarkerKmersData;
    void constructFromMarkerKmersThreadFunction(uint64_t threadId);

    // Data and functions used when constructing the Anchors from het sites using variant clustering.
    // New code that gets the Kmers from VariantClustering.
    class ConstructFromHetSitesData {
    public:
        uint64_t minClusterCoverage;
        uint64_t minAlleleCoverage;
        

        const std::vector<uint64_t>* clusterRepresentatives;
        DisjointSets* disjointSets;
        std::vector< std::vector<uint64_t> > membersByRepIdx;

        const MemoryMapped::Vector< std::pair<OrientedReadId, uint32_t> >* positionPairs;
        const MemoryMapped::Vector<uint8_t>* positionPairAlleles;
        const MemoryMapped::Vector<VariantPositionContext>* positionPairContexts;

        // The MarkerInfo objects for the candidate anchors found by each thread.
        using MarkerInfo = MarkerKmers::MarkerInfo;

        // The anchors found by each thread.
        vector< shared_ptr<MemoryMapped::VectorOfVectors<MarkerInfo, uint64_t> > > threadAnchors;
    };
    ConstructFromHetSitesData constructFromHetSitesData;
    void constructFromHetSitesThreadFunction(uint64_t threadId);


#if 0
    // Data and functions used when constructing the Anchors from marker k-mers (old code).
    class ConstructFromMarkerKmersData {
    public:
        uint64_t minPrimaryCoverage;
        uint64_t maxPrimaryCoverage;

        class MarkerInfo {
        public:
            OrientedReadId orientedReadId;
            uint32_t ordinal;
        };

        // A hash table that will contain a MarkerInfo object
        // for each marker in all oriented reads.
        // Indexed by bucketId.
        // The bucket is computed by hashing the k-mer of each marker,
        // so all markers with the same k-mer end up in the same bucket.
        MemoryMapped::VectorOfVectors<MarkerInfo, uint64_t> buckets;

        // The number of buckets is chosen equal to a power of 2,
        // so the bucketId can be obtained with a siple bitwise and
        // with a mask equal to the number of buckets minus 1.
        uint64_t mask;

        // The MarkerInfo objects for the candidate anchors found by each thread.
        vector< shared_ptr<MemoryMapped::VectorOfVectors<MarkerInfo, uint64_t> > > threadAnchors;

        // Another hash table where we information for each Kmer.
        class KmerInfo {
        public:
            Kmer kmer;
            uint64_t isForbidden:1;
            uint64_t frequency:63;
            KmerInfo(const Kmer kmer, uint64_t frequencyArgument) : kmer(kmer)
            {
                isForbidden = 0;
                frequency = frequencyArgument & 0x7FFFFFFFFFFFFFFFUL;
            }
            KmerInfo()
            {
                isForbidden = 0;
                frequency = 0;
            }
        };
        MemoryMapped::VectorOfVectors<KmerInfo, uint64_t > kmerInfo;
    };
    ConstructFromMarkerKmersData constructFromMarkerKmersData;
    void constructFromMarkerKmersGatherMarkersPass1(uint64_t threadId);
    void constructFromMarkerKmersGatherMarkersPass2(uint64_t threadId);
    void constructFromMarkerKmersGatherMarkersPass12(uint64_t pass);
    void constructFromMarkerKmersComputeKmerFrequencyPass1(uint64_t threadId);
    void constructFromMarkerKmersComputeKmerFrequencyPass2(uint64_t threadId);
    void constructFromMarkerKmersComputeKmerFrequencyPass12(uint64_t pass);
    void constructFromMarkerKmersFlagForbiddenKmers(uint64_t threadId);
    void constructFromMarkerKmersCreateAnchors(uint64_t threadId);
#endif


    // Process a candidate anchor from json input.
    using Ptree = boost::property_tree::basic_ptree<string, string>;
    bool processCandidateAnchor(const Ptree&, const string& name);
};



// Information about the read composition similarity of two anchors A and B.
class shasta::mode3::AnchorPairInfo {
public:

    // The total number of OrientedReadIds in each of the anchors A and B.
    uint64_t totalA = 0;
    uint64_t totalB = 0;

    // The number of common oriented reads.
    uint64_t common = 0;

    // The number of oriented reads present in A but not in B.
    uint64_t onlyA = 0;

    // The number of oriented reads present in B but not in A.
    uint64_t onlyB = 0;

    // The rest of the statistics are only valid if the number
    // of common oriented reads is not 0.

    // The estimated offset between the two Anchors.
    // The estimate is done using the common oriented reads.
    int64_t offsetInMarkers = invalid<int64_t>;
    int64_t offsetInBases = invalid<int64_t>;

    // The number of onlyA reads which are too short to be on edge B,
    // based on the above estimated offset.
    uint64_t onlyAShort = invalid<uint64_t>;

    // The number of onlyB reads which are too short to be on edge A,
    // based on the above estimated offset.
    uint64_t onlyBShort = invalid<uint64_t>;

    uint64_t intersectionCount() const
    {
        return common;
    }
    uint64_t unionCount() const {
        return totalA + totalB - common;
    }
    uint64_t correctedUnionCount() const
    {
        return unionCount() - onlyAShort - onlyBShort;
    }
    double jaccard() const
    {
        return double(intersectionCount()) / double(unionCount());
    }
    double correctedJaccard() const
    {
        return double(intersectionCount()) / double(correctedUnionCount());
    }

    void reverse()
    {
        swap(totalA, totalB);
        swap(onlyA, onlyB);
        swap(onlyAShort, onlyBShort);
        offsetInMarkers = - offsetInMarkers;
        offsetInBases = - offsetInBases;
    }

};
