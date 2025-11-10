/********************************************************************************

Creation of Anchors from json input.

This expects on input a json file containing a list of candidate anchors,
where each candidate anchor is a list of base intervals, and each base interval
is of the form:

["readName", strand, begin, end]
where:
- readName identifies the read using its name as in the input fasta or fastq.
- strand can be 1 or 0 to indicate whether the read is reverse complemented or not.
- begin is the first base position in the read that belongs to the anchor.
- end if the first base position in the read after the anchor.
If strand is "-", begin and end refer to base positions after reverse complement\ing,
so end is always > begin.

All the anchor read sequences in each candidate anchor must be exactly identical.

Candidate anchors are clipped to the first and last markers entirely contained in the candidate anchor.
If no markers are entirely contained in the candidate anchor, the candidate anchor is discarded.
Otherwise, the clipped candidate anchor is used to generate a pair of referce complemented anchors.

So in summary the input looks something like this (spacing only used for readability):
[
    [
        ["read1", "+". 400, 450],
        ["read2", "-", 2000, 2050],
        ["read3", "+", 6100 ,6150]
    [
        ["read4", "-". 12300, 12400],
        ["read5", "-", 4300, 4400],
        ["read6", "+", 8470 ,8570]
    ]
]

********************************************************************************/

// Shasta.
#include "mode3-Anchor.hpp"
#include "Base.hpp"
#include "Marker.hpp"
#include "Reads.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

// Standatd library.
#include "fstream.hpp"



Anchors::Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const vector<string>& jsonFileNames,
    uint64_t minPrimaryCoverage,
    uint64_t maxPrimaryCoverage,
    uint64_t /* threadCount */) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    k(k),
    markers(markers)
{
    kHalf = k / 2;

    // Initialize anchor data structures.
    anchorMarkerIntervals.createNew(
            largeDataName("AnchorMarkerIntervals"),
            largeDataPageSize);
    anchorSequences.createNew(
        largeDataName("AnchorSequences"), largeDataPageSize);
    anchorInfos.createNew(largeDataName("AnchorInfos"), largeDataPageSize);

    // Open a csv file that will contain the anchors generated for each candidate anchor.
    ofstream csv("AnchorsFromJson.csv");


    // Loop over the json input files.
    uint64_t candidateAnchorCount = 0;
    uint64_t candidateAnchorKeptCount = 0;
    uint64_t candidateAnchorDiscardedDueToCoverageCount = 0;
    uint64_t candidateAnchorDiscardedDueToLengthCount = 0;
    for(const string& jsonFileName: jsonFileNames) {

        // Read the json file into a boost property tree.
        Ptree json;
        boost::property_tree::read_json(jsonFileName, json);

        cout << "File " << jsonFileName << " contains " << json.size() << " candidate anchors." << endl;



        // Loop over the candidate anchors.
        for(const auto& p: json) {
            const Ptree& candidateAnchorWithName = p.second;
            ++candidateAnchorCount;

            try {

                // Sanity checks, then get the name.
                if(candidateAnchorWithName.size() != 2) {
                    cout << "Invalid format for the following candidate anchor. "
                        "Must have two elements." << endl;
                    boost::property_tree::write_json(cout, candidateAnchorWithName);
                    cout << endl;
                    throw runtime_error("Invalid format for the following candidate anchor.");
                }

                auto it = candidateAnchorWithName.begin();
                const Ptree& namePtree = (*it).second;

                const string name = namePtree.get<string>("");

                if(name.empty()) {
                    cout << "Empty name for the following candidate anchor." << endl;
                    boost::property_tree::write_json(cout, candidateAnchorWithName);
                    cout << endl;
                    throw runtime_error("Empty name for candidate anchor");
                }

                // Get the list of the base intervals.
                ++it;
                const Ptree& candidateAnchor = (*it).second;

                // If not in the desired coverage range, skip it.
                const uint64_t coverage = candidateAnchor.size();
                if(
                    (coverage < minPrimaryCoverage) or
                    (coverage > maxPrimaryCoverage)) {
                    ++candidateAnchorDiscardedDueToCoverageCount;
                    csv << name << "," << "Discarded due to coverage," << coverage << "\n";
                    continue;
                }

                if(processCandidateAnchor(candidateAnchor, name)) {
                    ++candidateAnchorKeptCount;
                    csv << name << "," << anchorIdToString(anchorInfos.size() - 2) << "," <<
                        anchorIdToString(anchorInfos.size() - 1) << "\n";
                } else {
                    ++candidateAnchorDiscardedDueToLengthCount;
                    csv << name << "," << "Discarded when clipping to markers\n";
                }
            } catch(...) {
                cout << "An error occurred while processing the following candidate anchor:" << endl;
                boost::property_tree::write_json(cout, candidateAnchorWithName);
                cout << endl;
                throw;
            }
        }
    }
    const uint64_t anchorCount = anchorInfos.size();

    cout << "Of " << candidateAnchorCount << " candidate anchors on input, " << candidateAnchorKeptCount <<
        " were kept, " << candidateAnchorDiscardedDueToCoverageCount << " were discarded due to coverage, and " <<
        candidateAnchorDiscardedDueToLengthCount << " were discarded due to length." << endl;
    cout << "Number of anchors created is " << anchorCount << endl;
    cout << "Anchor coverage in use is " << minPrimaryCoverage << " through " <<
        maxPrimaryCoverage << " included." << endl;

    // Sanity checks.
    SHASTA_ASSERT(
        candidateAnchorKeptCount +
        candidateAnchorDiscardedDueToCoverageCount +
        candidateAnchorDiscardedDueToLengthCount ==
        candidateAnchorCount);
    SHASTA_ASSERT(anchorCount == 2 * candidateAnchorKeptCount);
    SHASTA_ASSERT(anchorMarkerIntervals.size() == anchorCount);
    SHASTA_ASSERT(anchorSequences.size() == anchorCount);
}



// Process a candidate anchor from json input.
bool Anchors::processCandidateAnchor(
    const Ptree& candidateAnchor,
    const string& name)
{
    const bool debug = false;
    if(debug) {
        cout << "Processing a candidate anchor with coverage " << candidateAnchor.size() << endl;
    }
    using Ptree = boost::property_tree::basic_ptree<string, string>;

    class Interval {
    public:
        OrientedReadId orientedReadId;
        uint64_t begin;
        uint64_t end;
        uint64_t clippedBegin;
        uint64_t clippedEnd;
        uint32_t ordinal0;
        uint32_t ordinal1;
        uint64_t length() const
        {
            return end - begin;
        }
        uint64_t leftClip() const
        {
            return clippedBegin - begin;
        }
        uint64_t rightClip() const
        {
            return end - clippedEnd;
        }
        bool operator<(const Interval& that) const
        {
            return orientedReadId < that.orientedReadId;
        }
    };

    // Gather the intervals of this candidate anchor.
    vector<Interval> intervals;
    for(const auto& p: candidateAnchor) {
        const Ptree& candidateAnchorInterval = p.second;

        // The interval must have 4 entries (read name strand, begin, end).
        if(candidateAnchorInterval.size() != 4) {
            std::ostringstream s;
            boost::property_tree::write_json(s, candidateAnchorInterval);
            throw runtime_error("Anchor interval has size " + to_string(candidateAnchorInterval.size()) +
                ", must be 4: " + s.str());
        }

        // Parse the interval.
        try {
            // Get the read name.
            auto it = candidateAnchorInterval.begin();
            const string readName = ((it++)->second).get<string>("");

            // Get the ReadId.
            const ReadId readId = reads.getReadId(readName);
            if(readId == invalidReadId) {
                cout << "Read " << readName << " does not exist." << endl;
                throw runtime_error("Read does not exist.");
            }

            // Get the Strand
            const Strand strand = ((it++)->second).get<Strand>("");
            if((strand != 0) and (strand != 1)) {
                cout << "Invalid strand." << endl;
                throw runtime_error("Invalid strand.");
            }
            const OrientedReadId orientedReadId(readId, strand);

            // Get the begin, end of the interval.
            const uint64_t readLength = reads.getReadRawSequenceLength(readId);
            const uint64_t begin = ((it++)->second).get<uint64_t>("");
            const uint64_t end = ((it++)->second).get<uint64_t>("");
            if((begin >= readLength) or (end > readLength)) {
                cout << "Invalid begin/end for " << readName << endl;
                cout << "Read length is " << readLength << endl;
                throw runtime_error("Invalid begin/end.");
            }
            intervals.push_back({orientedReadId, begin, end});

            if(debug) {
                cout << readName << " " << orientedReadId << " " << begin << " " << end << endl;
            }

        } catch (...) {
            std::ostringstream s;
            boost::property_tree::write_json(s, candidateAnchorInterval);
            throw runtime_error("Invalid anchor interval : " + s.str());
        }
    }


    // Check that the sequences on these intervals are identical.
    const Interval& interval0 = intervals[0];
    const uint64_t length0 = interval0.length();
    bool foundDiscrepancy = false;
    if(length0 < k) {
    	if(debug) {
    		cout << "Anchor is too short." << endl;
    	}
        return false;
    }
    for(uint64_t i=1; i<intervals.size(); i++) {
        const Interval& interval = intervals[i];

        // Check the lengths.
        const uint64_t length = interval.length();
        if(length != length0) {
            std::ostringstream s;
            cout << "Length for interval " << i << " is " << length << endl;
            cout << "Length for interval 0 is " << length0 << endl;
            throw runtime_error("Interval lengths must all be identical.");
        }

        // Check the bases.
        for(uint64_t j=0; j<length; j++) {
            const Base b0 = reads.getOrientedReadBase(interval0.orientedReadId, uint32_t(interval0.begin + j));
            const Base b = reads.getOrientedReadBase(interval.orientedReadId, uint32_t(interval.begin + j));
            if(debug) {
            	cout << "Checking sequence " <<
					interval0.orientedReadId << " " << interval.orientedReadId << " "  <<
					j << " " << interval0.begin + j << " " << interval.begin + j << " " <<
					b0 << " " << b << endl;
            }
            if(b != b0) {
            	if(debug) {
            		cout << "Found sequence discrepancy." << endl;
            	}
                foundDiscrepancy = true;
                break;
            }
        }
        if(foundDiscrepancy) {
            break;
        }
    }



    // If a discrepancy in the sequence was found, write all the sequences.
    if(foundDiscrepancy) {
        cout << "The sequences of all intervals for anchor " << name << " are not identical:" << endl;
        for(const Interval& interval: intervals) {
            const uint64_t length = interval.length();
            for(uint64_t j=0; j<length; j++) {
                const Base b = reads.getOrientedReadBase(interval.orientedReadId, uint32_t(interval.begin + j));
                cout << b;
            }
            cout << " ";
            const span<const char> readName = reads.getReadName(interval.orientedReadId.getReadId());
            copy(readName.begin(), readName.end(), ostream_iterator<char>(cout));
            cout <<
                " " << interval.orientedReadId.getStrand() <<
                " " << interval.orientedReadId <<
                " " << interval.begin <<
                " " << interval.end << endl;
        }
        throw runtime_error("The sequences of all intervals for an anchor are not identical.");
    }



    // Object to compare CompressedMarkers by position.
    class Comparator {
    public:
        bool operator()(
            const CompressedMarker& x,
            const CompressedMarker& y) const
        {
            return x.position < y.position;
        }
    };
    const Comparator comparator;

    // Clip to the first/last marker entirely contained in each interval.
    for(Interval& interval: intervals) {
        if(debug) {
            cout << "Clipping " << interval.orientedReadId << " " << interval.begin << " " << interval.end << endl;
        }
        const span<const CompressedMarker> orientedReadMarkers = markers[interval.orientedReadId.getValue()];

        // Find the ordinal of the first marker entirely contained in this interval.
        auto it = std::lower_bound(
            orientedReadMarkers.begin(), orientedReadMarkers.end(),
            CompressedMarker({Uint24(uint32_t(interval.begin))}), comparator);
        if(it == orientedReadMarkers.end()) {
            if(debug) {
                cout << "The interval begins after the last marker." << endl;
            }
            return false;
        }

        const uint32_t ordinal0 = uint32_t(it - orientedReadMarkers.begin());
        const uint32_t position0 = it->position;
        if(position0 + k > interval.end) {
            if(debug) {
                cout << "No marker is entirely inside this interval." << endl;
            }
            return false;
        }

        // Ok, the marker at ordinal0 is entirely inside the interval.
        // Find the ordinal of the last marker entirely contained in this interval.
        uint32_t ordinal1 = ordinal0;
        uint32_t position1 = position0;
        for(; ordinal1<orientedReadMarkers.size(); ++ordinal1) {
            const uint32_t newPosition1 = orientedReadMarkers[ordinal1].position;
            if(newPosition1 + k > interval.end) {
                break;
            }
            position1 = newPosition1;
        }
        ordinal1 -= 1;
        position1 += uint32_t(k);
        if(debug) {
            cout << "Clipped to ordinals " << ordinal0 << " " << ordinal1 <<
                ", positions " << position0 << " " << position1 << endl;
        }

        interval.ordinal0 = ordinal0;
        interval.ordinal1 = ordinal1;
        interval.clippedBegin = position0;
        interval.clippedEnd = position1;
    }

    if(debug) {
        cout << "Clipped bases: " << intervals.front().leftClip() << " on left, " <<
            intervals.front().rightClip() << " on right." <<endl;
    }

    // Check that the number of bases clipped is the same for all the intervals.
    const uint64_t leftClip0 = intervals[0].leftClip();
    const uint64_t rightClip0 = intervals[0].rightClip();
    for(uint64_t i=1; i<intervals.size(); i++) {
        const Interval& interval = intervals[i];
        if((interval.leftClip() != leftClip0) or (interval.rightClip() != rightClip0)) {
            cout << "i,OrientedReadId,Begin,End,Length,ClipperBegin,ClippedEnd,Ordinal0,Ordinal1,LeftClip,RightClip," << endl;
            for(uint64_t i=0; i<intervals.size(); i++) {
                const Interval& interval = intervals[i];
                cout << i << ",";
                cout << interval.orientedReadId << ",";
                cout << interval.begin << ",";
                cout << interval.end << ",";
                cout << interval.length() << ",";
                cout << interval.clippedBegin << ",";
                cout << interval.clippedEnd << ",";
                cout << interval.ordinal0 << ",";
                cout << interval.ordinal1 << ",";
                cout << interval.leftClip() << ",";
                cout << interval.rightClip() << ",";
                cout << endl;
            }
            throw runtime_error("Clip inconsistency for interval " + to_string(i) +
                ". See above message for details.");
        }
    }


    // Check that the same ReadId does not appear twice.
    sort(intervals.begin(), intervals.end());
    for(uint64_t i=1; i<intervals.size(); i++) {
        const ReadId readId = intervals[i].orientedReadId.getReadId();
        const ReadId previousReadId = intervals[i-1].orientedReadId.getReadId();
        if(readId == previousReadId) {
            const auto readName = reads.getReadName(readId);
            string readNameString;
            copy(readName.begin(), readName.end(), back_inserter(readNameString));
            cout << "Duplicate read " << readNameString << " on anchor " << name << ". This anchor was suppressed." << endl;
            return false;
        }
    }

    // All good. Generate an anchor.
    {

        if(debug) {
            const AnchorId anchorId = anchorSequences.size();
            cout << "Generating anchor " << anchorIdToString(anchorId) << endl;
        }

        // Generate the sequence for this anchor.
        anchorSequences.appendVector();
        const Interval& interval0 = intervals[0];
        const OrientedReadId orientedReadId0 = interval0.orientedReadId;
        const uint64_t begin0 = interval0.clippedBegin + kHalf;
        const uint64_t end0 = interval0.clippedEnd - kHalf;
        for(uint64_t position=begin0; position<end0; ++position) {
            const Base b = reads.getOrientedReadBase(orientedReadId0, uint32_t(position));
            anchorSequences.append(b);
        }
        if(debug) {
            cout << "Anchor sequence:" << endl;
            const auto sequence = anchorSequences[anchorSequences.size() - 1];
            copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(cout));
            cout << endl;
        }

        // Generate the marker intervals for this anchor.
        anchorMarkerIntervals.appendVector();
        for(const Interval& interval: intervals) {
            anchorMarkerIntervals.append(AnchorMarkerInterval(interval.orientedReadId, interval.ordinal0));
            if(debug) {
                cout << interval.orientedReadId << " " << interval.ordinal0 << " " << interval.ordinal1 << endl;
            }
        }

        // Generate the AnchorInfo and store the ordinal offset.
        anchorInfos.resize(anchorInfos.size() + 1);
        AnchorInfo& anchorInfo = anchorInfos.back();
        anchorInfo.ordinalOffset = interval0.ordinal1 - interval0.ordinal0;



        // Also generate the reverse complemented anchor.
        anchorSequences.appendVector();
        if(end0 > begin0) {
            for(uint64_t position=end0-1; /* Check later */ ; --position) {
                const Base b = reads.getOrientedReadBase(orientedReadId0, uint32_t(position));
                anchorSequences.append(b.complement());
                if(position == begin0) {
                    break;
                }
            }
        }

        anchorMarkerIntervals.appendVector();
        for(const Interval& interval: intervals) {
            const uint32_t markerCount = uint32_t(markers[interval.orientedReadId.getValue()].size());
            AnchorMarkerInterval markerInterval(interval.orientedReadId, markerCount - 1 - interval.ordinal1);
            markerInterval.orientedReadId.flipStrand();
            anchorMarkerIntervals.append(markerInterval);
        }

        anchorInfos.resize(anchorInfos.size() + 1);
        anchorInfos.back().ordinalOffset = interval0.ordinal1 - interval0.ordinal0;
    }


    return true;
}

